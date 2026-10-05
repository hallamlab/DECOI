"""Build linked WGS/V4 reference bundles directly from genome FASTAs."""
import gzip
import hashlib
import json
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

RANKS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
IUPAC = dict(zip("ACGTRYSWKMBDHVN", ["A", "C", "G", "T", "AG", "CT", "CG", "AT", "GT", "AC", "CGT", "AGT", "ACT", "ACG", "ACGT"]))


def primer_pattern(primer):
    if not primer or set(primer.upper()) - set(IUPAC):
        raise ValueError("Primer must be a nonempty IUPAC DNA sequence")
    return re.compile("(?=(" + "".join("[" + IUPAC[c] + "]" for c in primer.upper()) + "))")


def extract_markers(sequence, forward, reverse, minimum, maximum):
    """Return primer-bounded inserts on both strands, with zero-based coordinates."""
    sequence = sequence.upper()
    fwd = primer_pattern(forward)
    rev = primer_pattern(str(Seq(reverse).reverse_complement()))
    hits = []
    for strand, oriented in (("+", sequence), ("-", str(Seq(sequence).reverse_complement()))):
        for match in fwd.finditer(oriented):
            start = match.start() + len(forward)
            # Search only the valid interval; avoid a genome-wide all-pairs scan.
            for end_match in rev.finditer(oriented, start + minimum, min(len(oriented), start + maximum + len(reverse))):
                end = end_match.start()
                marker = oriented[start:end]
                if minimum <= len(marker) <= maximum and not set(marker) - set("ACGT"):
                    left, right = (start, end) if strand == "+" else (len(sequence)-end, len(sequence)-start)
                    hits.append((marker, strand, left, right))
    return hits


def prepare_genome_reference(cfg, source, output):
    """Discover genomes, extract markers, and write the existing DECOI reference schema.

    Initial selection retains one unique V4 allele per genome and one genome per
    allele. Every excluded genome is recorded with its reason. No genome sequence
    is fabricated, padded, or altered to make a marker match.
    """
    source, output = Path(source).resolve(), Path(output).resolve()
    if not source.is_dir():
        raise ValueError(f"Genome source directory does not exist: {source}")
    if output.exists():
        raise ValueError(f"Use a fresh reference output directory: {output}")
    rcfg = cfg["reference"]
    minimum, maximum = int(rcfg.get("min_length", 230)), int(rcfg.get("max_length", 310))
    if minimum < 1 or maximum < minimum:
        raise ValueError("Invalid marker length interval")
    primer_pattern(rcfg["forward_primer"]); primer_pattern(rcfg["reverse_primer"])
    files = sorted(p for p in source.rglob("*") if p.is_file() and any(str(p).endswith(s) for s in (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")))
    if not files:
        raise ValueError("No genome FASTAs found (one assembly per .fa/.fasta/.fna file, optionally gzipped)")
    taxonomy = {}
    taxfile = source / "taxonomy.tsv"
    if taxfile.exists():
        tax = pd.read_csv(taxfile, sep="\t", dtype=str).fillna("")
        if "fasta" not in tax or tax.fasta.duplicated().any():
            raise ValueError("taxonomy.tsv requires unique fasta paths relative to the source directory")
        taxonomy = tax.set_index("fasta").to_dict("index")
    silva = {}
    if rcfg.get("taxonomy_reference_dir"):
        tax = pd.read_csv(Path(rcfg["taxonomy_reference_dir"]) / "v4_asv_taxonomy.tsv", sep="\t", dtype=str).fillna("")
        silva = tax.set_index("ASV_ID").to_dict("index")
    output.mkdir(parents=True)
    bundle = output / "genomes"
    bundle.mkdir()
    audit, registry, taxrows, markers, coordinates = [], [], [], {}, []
    filt = cfg["simulation"].get("selection", {})
    for path in files:
        relative = path.relative_to(source).as_posix()
        checksum = hashlib.sha256(path.read_bytes()).hexdigest()
        gid = "GENOME_" + checksum[:20]
        with (gzip.open(path, "rt") if path.suffix == ".gz" else path.open()) as handle:
            records = list(SeqIO.parse(handle, "fasta"))
        reason = ""
        hits = []
        if not records or any(len(r) < 1000 or set(str(r.seq).upper()) - set("ACGT") for r in records):
            reason = "requires_contigs_at_least_1000bp_and_unambiguous_ACGT"
        else:
            for i, record in enumerate(records):
                hits.extend((i, *h) for h in extract_markers(str(record.seq), rcfg["forward_primer"], rcfg["reverse_primer"], minimum, maximum))
            alleles = {h[1] for h in hits}
            if not alleles:
                reason = "no_exact_primer_bounded_V4"
            elif len(alleles) > 1:
                reason = "multiple_V4_alleles_not_yet_supported"
        aid, marker, ranks = "", "", {}
        if not reason:
            marker = hits[0][1]
            aid = "ASV_" + hashlib.sha256(marker.encode()).hexdigest()[:16]
            ranks = {r: taxonomy.get(relative, {}).get(r) or silva.get(aid, {}).get(r, "") for r in RANKS}
            copies = sum(str(r.seq).upper().count(marker) + (str(r.seq).upper().count(str(Seq(marker).reverse_complement())) if str(Seq(marker).reverse_complement()) != marker else 0) for r in records)
            if copies != len(hits):
                reason = "V4_occurrences_disagree_with_primer_bounded_loci"
            elif filt.get("domains") and ranks["Domain"] not in filt["domains"]:
                reason = "domain_filter_or_missing_taxonomy"
            elif sum(bool(v) for v in ranks.values()) < int(filt.get("min_taxonomy_ranks", 0)):
                reason = "insufficient_taxonomy"
            elif aid in markers:
                reason = "V4_shared_with_previously_retained_genome"
        audit.append({"source_fasta": relative, "source_sha256": checksum, "genome_id": gid,
                      "ASV_ID": aid, "primer_bounded_loci": len(hits), "status": reason or "retained"})
        if reason:
            continue
        markers[aid] = marker
        role=taxonomy.get(relative,{}).get("reference_role", "")
        if role and role not in {"biological", "contaminant"}:
            raise ValueError(f"Invalid reference_role for {relative}: {role}")
        role_fields={"reference_role":role} if role else {}
        provenance_fields={key:taxonomy.get(relative,{}).get(key,"") for key in ("accession","taxid","species_taxid","priority") if key in taxonomy.get(relative,{})}
        filename = gid + ".fasta"
        SeqIO.write(records, bundle / filename, "fasta")
        registry.append({"ASV_ID": aid, "genome_id": gid, "fasta": filename, "marker_copies": copies,
                         "source_fasta": relative, "source_sha256": checksum, **ranks, **role_fields, **provenance_fields})
        taxrows.append({"ASV_ID": aid, **ranks, "representative_reference_id": gid, "n_identical_references": copies, **role_fields, **provenance_fields})
        coordinates.extend({"ASV_ID": aid, "genome_id": gid, "contig_index": i+1, "source_contig_id": records[i].id,
                            "strand": strand, "start_0based": start, "end_exclusive": end} for i, _, strand, start, end in hits)
    pd.DataFrame(audit).to_csv(output / "genome_selection.tsv", sep="\t", index=False)
    if not registry:
        raise ValueError(f"No eligible genomes; see {output / 'genome_selection.tsv'} for exclusion reasons")
    pd.DataFrame(registry).to_csv(bundle / "genomes.tsv", sep="\t", index=False)
    pd.DataFrame(taxrows).to_csv(output / "v4_asv_taxonomy.tsv", sep="\t", index=False)
    pd.DataFrame(coordinates).to_csv(output / "genome_marker_loci.tsv", sep="\t", index=False)
    pd.DataFrame([{"ASV_ID": r["ASV_ID"], "reference_id": r["genome_id"]} for r in registry]).to_csv(output / "v4_asv_source_map.tsv", sep="\t", index=False)
    SeqIO.write([SeqRecord(Seq(s), id=a, description="") for a, s in markers.items()], output / "v4_asv_pool.fasta", "fasta")
    metadata = {"source": str(source), "source_type": "genomes", "n_source_genomes": len(files), "n_unique_v4": len(markers),
                "genome_bundle": "genomes", "extraction": "exact IUPAC primer matches, both strands, linear contigs",
                "selection_policy": "one V4 allele per genome; one genome per allele; configured taxonomy filters",
                "primers": {"forward": rcfg["forward_primer"], "reverse": rcfg["reverse_primer"]},
                "taxonomy_source_sha256": hashlib.sha256(taxfile.read_bytes()).hexdigest() if taxfile.exists() else None}
    (output / "reference_manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata
