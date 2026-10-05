#!/usr/bin/env python3
"""Freeze an airway-prioritized RefSeq panel and verify it for paired simulation.

Downloads are resumable; accepted accessions and checksums are frozen in a lock
file. Rehydration uses only that lock, never a fresh catalog query.
"""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

import pandas as pd
import yaml
from Bio import SeqIO
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from genome_reference import extract_markers, RANKS


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def download(url, path):
    path = Path(path)
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part")
    subprocess.run(["curl", "--fail", "--location", "--silent", "--show-error", "--retry", "3",
                    "--connect-timeout", "30", "--max-time", "180", url, "-o", str(temporary)], check=True)
    temporary.replace(path)


def ranked_taxonomy(archive, needed):
    result = {}
    with tarfile.open(archive, "r:gz") as tar:
        name = next(m for m in tar.getmembers() if m.name.endswith("rankedlineage.dmp"))
        with tar.extractfile(name) as stream:
            for raw in stream:
                parts = [p.strip() for p in raw.decode().split("|")]
                if parts[0] not in needed:
                    continue
                # NCBI rankedlineage's established lower-rank order is fixed;
                # Domain is identified by name to accommodate the rank migration.
                result[parts[0]] = dict(zip(RANKS,
                    [next((x for x in parts[8:] if x in {"Bacteria", "Archaea"}), ""),
                     parts[7], parts[6], parts[5], parts[4], parts[3], parts[2] or parts[1]]))
    return result


def candidate_plan(policy, summary, taxonomy_archive):
    eligible = []
    with open(summary) as stream:
        header = None
        for line in stream:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").strip().split("\t")
                continue
            if header is None:
                raise ValueError("NCBI assembly summary header missing")
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            if (row.get("group") not in {"bacteria", "archaea"}
                or row.get("assembly_level") not in policy["assembly_levels"]
                or row.get("version_status") != "latest"
                or row.get("excluded_from_refseq") not in {"na", ""}
                or not row.get("ftp_path", "").startswith("https://ftp.ncbi.nlm.nih.gov/")
                or not 300000 <= int(row.get("genome_size", "0")) <= policy["max_genome_bases"]):
                continue
            eligible.append(row)
    taxonomy = ranked_taxonomy(taxonomy_archive, {r["species_taxid"] for r in eligible})
    eligible.sort(key=lambda r: (r["refseq_category"] not in {"reference genome", "representative genome"}, r["assembly_accession"]))
    species = Counter()
    groups = defaultdict(list)
    for row in eligible:
        sid = row["species_taxid"]
        tax = taxonomy.get(sid, {})
        if tax.get("Domain") not in policy["domains"] or sum(bool(tax.get(r)) for r in RANKS) < 5:
            continue
        if species[sid] >= policy["max_assemblies_per_species"]:
            continue
        attempt = species[sid]; species[sid] += 1
        genus = tax["Genus"]
        priority = genus in policy["priority_genera"]
        contaminant = genus in policy.get("contaminant_genera", [])
        url = row["ftp_path"].rstrip("/")
        item = {"accession": row["assembly_accession"], "taxid": row["taxid"], "species_taxid": sid,
                "organism_name": row["organism_name"], "refseq_category": row["refseq_category"],
                "assembly_level": row["assembly_level"], "genome_size": int(row["genome_size"]),
                "priority": "reagent_contaminant" if contaminant else ("airway_prioritized" if priority else "background"), "attempt": attempt,
                "reference_role": "contaminant" if contaminant else "biological",
                "url": url + "/" + url.split("/")[-1] + "_genomic.fna.gz",
                "md5_url": url + "/md5checksums.txt", **tax}
        groups[(0 if contaminant else (1 if priority else 2), genus)].append(item)
    def seeded(value):
        return hashlib.sha256(f"{policy['seed']}:{value}".encode()).hexdigest()
    ordered = []
    for (tier, genus), rows in groups.items():
        rows.sort(key=lambda r: (r["attempt"], seeded(r["species_taxid"]), r["accession"]))
        for index, row in enumerate(rows):
            ordered.append(((tier, index, seeded(genus)), row))
    return [row for _, row in sorted(ordered, key=lambda x:x[0])]


def fetch_and_screen(row, cache, reference):
    row = dict(row)
    accession = row["accession"]
    path = cache / f"{accession}.fna.gz"
    md5path = cache / f"{accession}.md5checksums.txt"
    try:
        download(row["md5_url"], md5path)
        remote_name = row["url"].split("/")[-1]
        checks = [line.split()[0] for line in md5path.read_text().splitlines()
                  if len(line.split()) == 2 and line.split()[1].removeprefix("./") == remote_name]
        if len(checks) != 1:
            raise ValueError("Unique NCBI genome MD5 entry not found")
        download(row["url"], path)
        if digest(path, "md5") != checks[0]:
            raise ValueError("Genome MD5 mismatch")
        row.update(ncbi_md5=checks[0], sha256=digest(path), bytes=path.stat().st_size)
        with gzip.open(path, "rt") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
        if not records or any(len(r) < 1000 or set(str(r.seq).upper()) - set("ACGT") for r in records):
            row["status"] = "short_or_ambiguous_contigs"
            return row
        hits = [h for r in records for h in extract_markers(str(r.seq), reference["forward_primer"],
                    reference["reverse_primer"], reference["min_length"], reference["max_length"])]
        alleles = {h[0] for h in hits}
        if not alleles:
            row["status"] = "no_exact_primer_bounded_V4"
        elif len(alleles) > 1:
            row["status"] = "multiple_V4_alleles"
        else:
            marker = hits[0][0]; reverse = str(Seq(marker).reverse_complement())
            copies = sum(str(r.seq).upper().count(marker) + (str(r.seq).upper().count(reverse) if marker != reverse else 0) for r in records)
            if copies != len(hits):
                row["status"] = "marker_copy_mismatch"
            else:
                row.update(status="eligible", ASV_ID="ASV_"+hashlib.sha256(marker.encode()).hexdigest()[:16],
                           marker_copies=copies, genome_length=sum(len(r) for r in records))
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        row.update(status="download_or_parse_error", error=str(exc))
    return row


def build(args):
    policy = yaml.safe_load(args.policy.read_text())
    cfg = yaml.safe_load(args.config.read_text())
    root = args.output.resolve(); root.mkdir(parents=True, exist_ok=True)
    metadata = root / "metadata"; metadata.mkdir(exist_ok=True)
    cache = root / "cache"; cache.mkdir(exist_ok=True)
    lock_path = root / "panel.lock.json"
    if lock_path.exists():
        raise ValueError("Panel already frozen; use verify or rehydrate, or choose a new panel version")
    download(policy["assembly_source"], metadata / "assembly_summary_refseq.txt")
    download(policy["taxonomy_source"], metadata / "new_taxdump.tar.gz")
    plan_path = metadata / "candidates.json"
    provenance = {"policy_sha256": digest(args.policy), "config_sha256": digest(args.config),
                  "assembly_summary_sha256": digest(metadata / "assembly_summary_refseq.txt"),
                  "taxonomy_sha256": digest(metadata / "new_taxdump.tar.gz")}
    if plan_path.exists():
        cached = json.loads(plan_path.read_text())
        if cached["provenance"] != provenance:
            raise ValueError("Inputs changed since candidate planning; use a new panel directory")
        plan = cached["candidates"]
    else:
        plan = candidate_plan(policy, metadata / "assembly_summary_refseq.txt", metadata / "new_taxdump.tar.gz")
        plan_path.write_text(json.dumps({"provenance": provenance, "candidates": plan}, indent=2)+"\n")
    print(f"Planned {len(plan)} candidates; {sum(r['priority']=='airway_prioritized' for r in plan)} in priority genera", flush=True)
    if args.plan_only:
        return
    selected, genus_counts, seen_species, seen_markers = [], Counter(), set(), set()
    audit_path = root / "screening.jsonl"
    previous = {}
    if audit_path.exists():
        for line in audit_path.read_text().splitlines():
            row = json.loads(line); previous[row["accession"]] = row
    def process(row):
        old = previous.get(row["accession"])
        if old and old["status"] != "download_or_parse_error":
            return old
        return fetch_and_screen(row, cache, cfg["reference"])
    screened = 0
    def cap_for(row):
        if row["reference_role"] == "contaminant":
            return 1
        return policy["max_retained_per_priority_genus"] if row["priority"]=="airway_prioritized" else policy["max_retained_per_background_genus"]
    def contaminants_full():
        return sum(r["reference_role"]=="contaminant" for r in selected) >= policy.get("contaminant_count", 0)
    with ThreadPoolExecutor(max_workers=args.workers) as pool, audit_path.open("a") as audit:
        for start in range(0, len(plan), 16):
            batch = []
            for row in plan[start:start+16]:
                cap = cap_for(row)
                if row["reference_role"]=="contaminant" and contaminants_full():
                    continue
                if row["species_taxid"] not in seen_species and genus_counts[row["Genus"]] < cap:
                    batch.append(row)
            for result in pool.map(process, batch):
                audit.write(json.dumps(result)+"\n"); audit.flush(); screened += 1
                if result["status"] != "eligible" or result["ASV_ID"] in seen_markers or result["species_taxid"] in seen_species:
                    continue
                cap = cap_for(result)
                if result["reference_role"]=="contaminant" and contaminants_full():
                    continue
                if genus_counts[result["Genus"]] >= cap:
                    continue
                selected.append(result); seen_species.add(result["species_taxid"]); seen_markers.add(result["ASV_ID"]); genus_counts[result["Genus"]]+=1
                if len(selected) == policy["target_genomes"]:
                    break
            print(f"Screened {screened}; retained {len(selected)}/{policy['target_genomes']}; priority {sum(r['priority']=='airway_prioritized' for r in selected)}", flush=True)
            if len(selected) == policy["target_genomes"] or screened >= policy["max_candidates"]:
                break
    if len(selected) < policy["target_genomes"]:
        raise RuntimeError(f"Only {len(selected)} eligible unique genomes; screening audit retained, panel not frozen")
    if not contaminants_full():
        raise RuntimeError("Not enough distinct reagent-contaminant genomes; panel not frozen")
    raw = root / "raw"; raw.mkdir(exist_ok=True)
    for row in selected:
        filename = row["accession"] + ".fna.gz"
        source = cache / filename
        if digest(source) != row["sha256"]:
            raise ValueError(f"Cached genome changed: {filename}")
        target = raw / filename
        if not target.exists():
            shutil.copyfile(source, target)
        row["fasta"] = filename
    pd.DataFrame(selected)[["fasta", *RANKS, "accession", "taxid", "species_taxid", "priority", "reference_role"]].to_csv(raw / "taxonomy.tsv", sep="\t", index=False)
    pd.DataFrame(selected).to_csv(root / "panel.tsv", sep="\t", index=False)
    lock = {"panel_id": policy["panel_id"], "created_utc": datetime.now(timezone.utc).isoformat(),
            "policy": policy, "provenance": provenance, "reference_settings": cfg["reference"],
            "taxonomy_table_sha256": digest(raw / "taxonomy.tsv"), "genomes": selected}
    lock_path.write_text(json.dumps(lock, indent=2)+"\n")
    print(f"Frozen {len(selected)} genomes in {lock_path}", flush=True)


def replay(args):
    lock = json.loads(args.lock.read_text())
    raw = args.output.resolve(); raw.mkdir(parents=True, exist_ok=True)
    for row in lock["genomes"]:
        path = raw / row["fasta"]
        if args.command == "rehydrate":
            download(row["url"], path)
        if digest(path) != row["sha256"]:
            raise ValueError(f"Checksum mismatch: {path}")
    if args.command == "rehydrate":
        pd.DataFrame(lock["genomes"])[["fasta", *RANKS, "accession", "taxid", "species_taxid", "priority", "reference_role"]].to_csv(raw / "taxonomy.tsv", sep="\t", index=False)
    if digest(raw / "taxonomy.tsv") != lock["taxonomy_table_sha256"]:
        raise ValueError("Taxonomy table checksum mismatch")
    print(f"Verified {len(lock['genomes'])} pinned genomes and taxonomy", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build")
    b.add_argument("--policy", type=Path, default=Path("references/airway_v1/policy.yaml"))
    b.add_argument("--config", type=Path, default=Path("config/airway_paired.yaml"))
    b.add_argument("--output", type=Path, default=Path("data/reference_panels/airway_v1"))
    b.add_argument("--workers", type=int, default=4)
    b.add_argument("--plan-only", action="store_true")
    for command in ("verify", "rehydrate"):
        p = sub.add_parser(command)
        p.add_argument("--lock", type=Path, required=True)
        p.add_argument("--output", type=Path, required=True, help="Raw genome directory")
    args = parser.parse_args()
    if args.command == "build": build(args)
    else: replay(args)
