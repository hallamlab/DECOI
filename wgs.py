"""Matched short-read metagenomics from DECOI's shared marker abundance state."""
from pathlib import Path
import hashlib
import re

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def load_genomes(directory, marker_sequences, selected_ids=None):
    """Load a one-genome/one-V4-allele reference bundle; verify sequence linkage."""
    directory = Path(directory).resolve()
    table = pd.read_csv(directory / "genomes.tsv", sep="\t", dtype=str).fillna("")
    required = {"ASV_ID", "genome_id", "fasta", "marker_copies"}
    if not required.issubset(table.columns) or table.empty:
        raise ValueError("genomes.tsv requires ASV_ID, genome_id, fasta, marker_copies and at least one row")
    if table.ASV_ID.duplicated().any() or table.genome_id.duplicated().any():
        raise ValueError("Initial WGS support requires a one-to-one ASV/genome mapping")
    if selected_ids is not None:
        missing=set(selected_ids)-set(table.ASV_ID)
        if missing:
            raise ValueError(f"Genome mappings missing for {len(missing)} existing features; study was not modified")
        table=table[table.ASV_ID.isin(selected_ids)]
    genomes, rows = {}, []
    for row in table.to_dict("records"):
        aid, gid = row["ASV_ID"], row["genome_id"]
        if not re.fullmatch(r"[A-Za-z0-9_-]+", gid):
            raise ValueError(f"Unsafe genome identifier: {gid!r}")
        if aid not in marker_sequences:
            raise ValueError(f"Mapped ASV absent from prepared reference: {aid}")
        path = (directory / row["fasta"]).resolve()
        if not path.is_relative_to(directory):
            raise ValueError("Genome FASTAs must reside inside the reference bundle")
        copies = int(row["marker_copies"])
        records = list(SeqIO.parse(path, "fasta"))
        if copies < 1 or not records:
            raise ValueError(f"Invalid marker copy number or empty FASTA: {gid}")
        marker = marker_sequences[aid]
        reverse = str(Seq(marker).reverse_complement())
        observed = 0
        for i, rec in enumerate(records):
            sequence = str(rec.seq).upper()
            if len(sequence) < 1000 or set(sequence) - set("ACGT"):
                raise ValueError(f"{gid}: each contig must contain >=1000 unambiguous ACGT bases")
            observed += sequence.count(marker)
            if reverse != marker:
                observed += sequence.count(reverse)
            rec.seq = Seq(sequence)
            rec.id = f"{gid}__contig_{i+1}"
            rec.description = ""
        if observed != copies:
            raise ValueError(f"{gid}: marker_copies={copies}, but found {observed} exact V4 occurrences")
        genomes[aid] = records
        rows.append({**row, "marker_copies": copies,
                     "genome_length": sum(len(r) for r in records),
                     "fasta_sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    return pd.DataFrame(rows).set_index("ASV_ID"), genomes


def allocate_reads(marker_counts, registry, read_pairs, rng):
    """Marker equivalents / copies = cell equivalents; DNA mass scales with length."""
    if isinstance(read_pairs, bool) or int(read_pairs) != read_pairs or read_pairs < 0:
        raise ValueError("wgs.read_pairs must be a nonnegative integer")
    if not marker_counts.index.isin(registry.index).all():
        raise ValueError("Missing genome mappings for shared community features")
    reg = registry.loc[marker_counts.index]
    cells = marker_counts.div(reg.marker_copies, axis=0)
    weights = cells.mul(reg.genome_length, axis=0)
    probabilities = weights.div(weights.sum(axis=0), axis=1).fillna(0)
    pairs = pd.DataFrame(0, index=marker_counts.index, columns=marker_counts.columns, dtype=int)
    for sample in pairs:
        if weights[sample].sum() > 0:
            pairs[sample] = rng.multinomial(int(read_pairs), probabilities[sample].to_numpy())
    return cells, probabilities, pairs


def simulate_wgs(cfg, biological, final, contaminants, registry, genomes, output,
                 run, count_fastq_records, stable_uuid):
    """Render WGS using an independent RNG, retaining shared samples and blanks."""
    settings = cfg["wgs"]
    directory = output / "wgs"
    directory.mkdir()
    inputs = directory / "iss_inputs"
    fastq = directory / "fastq"
    inputs.mkdir(); fastq.mkdir()
    ids = list(biological.index) + list(contaminants)
    shared = pd.DataFrame(0, index=ids, columns=final.columns, dtype=int)
    shared.loc[biological.index, biological.columns] = biological
    if contaminants:
        shared.loc[contaminants] = final.loc[contaminants]
    controls = [s for s in final.columns if s not in biological.columns]
    if controls:
        shared.loc[:, controls] = final.loc[ids, controls]
    seed = np.random.SeedSequence([int(cfg["seed"]), 574753])
    rng = np.random.default_rng(seed)
    reg = registry.loc[ids].copy()
    reg["feature_uuid"] = [stable_uuid("mock16s-chem-feature", a) for a in ids]
    reg["feature_type"] = ["contaminant" if a in contaminants else "biological" for a in ids]
    reg.to_csv(directory / "genome_registry.tsv", sep="\t")
    cells, probabilities, pairs = allocate_reads(shared, reg, settings.get("read_pairs", 100000), rng)
    for name, frame in [("shared_marker_counts", shared), ("cell_equivalents", cells),
                        ("expected_dna_fractions", probabilities), ("genome_read_pairs", pairs)]:
        frame.rename(index=reg.genome_id).to_csv(directory / f"{name}.tsv", sep="\t", index_label="genome_id")
    rows, contig_truth = [], []
    for sample in pairs.columns:
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", str(sample)) or sample in {".", ".."}:
            raise ValueError(f"Unsafe sample identifier: {sample!r}")
        fa, rc = inputs / f"{sample}.fasta", inputs / f"{sample}.readcounts.tsv"
        with fa.open("w") as fh, rc.open("w") as rh:
            for aid, n in pairs[sample].items():
                records = genomes[aid]
                lengths = np.array([len(r) for r in records], dtype=float)
                allocations = rng.multinomial(int(n), lengths / lengths.sum())
                for rec, count in zip(records, allocations):
                    contig_truth.append((sample, reg.loc[aid, "genome_id"], rec.id, int(count)))
                    if count:
                        SeqIO.write(rec, fh, "fasta")
                        rh.write(f"{rec.id}\t{2 * int(count)}\n")
        expected = int(pairs[sample].sum())
        prefix = fastq / sample
        compressed = bool(settings.get("gzip", True))
        suffix = ".fastq.gz" if compressed else ".fastq"
        r1, r2 = Path(f"{prefix}_R1{suffix}"), Path(f"{prefix}_R2{suffix}")
        sample_seed = int(rng.integers(1, 2**31-1))
        if expected:
            cmd = ["iss", "generate", "--genomes", str(fa), "--readcount_file", str(rc),
                   "--sequence_type", "metagenomics", "--model", str(settings.get("model", "miseq")),
                   "--cpus", str(settings.get("cpus", 1)), "--seed", str(sample_seed), "--output", str(prefix)]
            if compressed:
                cmd.append("--compress")
            run(cmd)
        else:
            import gzip
            for path in (r1, r2):
                with (gzip.open(path, "wt") if compressed else path.open("w")):
                    pass
        n1, n2 = count_fastq_records(r1), count_fastq_records(r2)
        if n1 != expected or n2 != expected:
            raise RuntimeError(f"WGS read-count mismatch for {sample}: expected {expected}, got {n1}/{n2}")
        rows.append((sample, expected, n1, n2, sample_seed, f"fastq/{r1.name}", f"fastq/{r2.name}"))
    validation = pd.DataFrame(rows, columns=["sample_id", "expected_pairs", "observed_R1", "observed_R2", "iss_seed", "fastq_r1", "fastq_r2"])
    validation.to_csv(directory / "fastq_validation.tsv", sep="\t", index=False)
    validation[["sample_id", "fastq_r1", "fastq_r2"]].to_csv(directory / "fastq_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(contig_truth, columns=["sample_id", "genome_id", "contig_id", "read_pairs"]).to_csv(directory / "contig_read_pairs.tsv", sep="\t", index=False)
    return {"samples": len(rows), "genomes": len(ids), "total_read_pairs": int(pairs.to_numpy().sum()),
            "reference_registry": "wgs/genome_registry.tsv", "abundance_basis": "pre-PCR marker equivalents / marker copies * genome length",
            "excluded_amplicon_features": ["chimeras", "mitochondrial marker-only fixtures"]}
