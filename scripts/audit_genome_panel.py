#!/usr/bin/env python3
"""Audit a frozen panel against its prepared DECOI reference and write a report."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.build_genome_panel import digest
from wgs import load_genomes
from Bio import SeqIO


def audit(lock_path, raw, prepared):
    lock = json.loads(lock_path.read_text())
    rows = lock["genomes"]
    assert len(rows) == lock["policy"]["target_genomes"]
    assert len({r["accession"] for r in rows}) == len(rows)
    assert len({r["ASV_ID"] for r in rows}) == len(rows)
    assert len({r["species_taxid"] for r in rows}) == len(rows)
    for row in rows:
        assert digest(raw / row["fasta"]) == row["sha256"], row["accession"]
        assert digest(raw / row["fasta"], "md5") == row["ncbi_md5"], row["accession"]
    assert digest(raw / "taxonomy.tsv") == lock["taxonomy_table_sha256"]
    markers = {r.id: str(r.seq) for r in SeqIO.parse(prepared / "v4_asv_pool.fasta", "fasta")}
    # Validates every marker against the real genome sequences and copy numbers.
    registry, genomes = load_genomes(prepared / "genomes", markers)
    assert set(registry.index) == {r["ASV_ID"] for r in rows}
    tax = pd.read_csv(prepared / "v4_asv_taxonomy.tsv", sep="\t").set_index("ASV_ID")
    for row in rows:
        r = registry.loc[row["ASV_ID"]]
        assert r.accession == row["accession"]
        assert r.source_sha256 == row["sha256"]
        assert r.reference_role == row["reference_role"]
        assert r.marker_copies == row["marker_copies"]
        assert r.genome_length == row["genome_length"]
        assert tax.loc[row["ASV_ID"], "reference_role"] == row["reference_role"]
    del genomes
    contaminants = [r for r in rows if r["reference_role"] == "contaminant"]
    assert len(contaminants) == lock["policy"]["contaminant_count"]
    assert len({r["Genus"] for r in contaminants}) == len(contaminants)
    report = {"panel_id":lock["panel_id"], "lock_sha256":digest(lock_path), "status":"passed",
              "genomes":len(rows), "roles":dict(Counter(r["reference_role"] for r in rows)),
              "priority_counts":dict(Counter(r["priority"] for r in rows)),
              "genus_counts":dict(sorted(Counter(r["Genus"] for r in rows).items())),
              "domain_counts":dict(Counter(r["Domain"] for r in rows)),
              "genome_bases":sum(r["genome_length"] for r in rows),
              "compressed_genome_bytes":sum(r["bytes"] for r in rows),
              "contaminants":[{k:r[k] for k in ("accession","organism_name","Genus","ASV_ID")} for r in contaminants],
              "checks":["NCBI MD5", "frozen SHA-256", "taxonomy table checksum", "unique accessions/species/markers",
                        "all prepared V4 sequences present in their genomes", "marker copy counts", "genome lengths",
                        "provenance and contaminant roles preserved"]}
    return report


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--lock",type=Path,required=True)
    p.add_argument("--raw",type=Path,required=True)
    p.add_argument("--prepared",type=Path,required=True)
    p.add_argument("--report",type=Path,required=True)
    args=p.parse_args()
    report=audit(args.lock,args.raw,args.prepared)
    args.report.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
