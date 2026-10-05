#!/usr/bin/env python3
"""Create small synthetic references for a local paired-assay integration test."""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mock16s_chem import stable_asv_id


def create(destination, genome_first=False):
    destination = destination.resolve()
    destination.mkdir(parents=True, exist_ok=False)
    bundle = destination / "genomes"
    bundle.mkdir()
    rng = np.random.default_rng(123)
    rows, silva = [], []
    for i in range(8):
        marker = "".join(rng.choice(list("ACGT"), 250))
        aid = stable_asv_id(marker)
        genome = "".join(rng.choice(list("ACGT"), 10000 + i*1000))
        molecule = "GTGCCAGCAGCCGCGGTAA" + marker + str(Seq("GGACTACAGGGGTATCTAAT").reverse_complement())
        genome = genome[:2000] + molecule + genome[2000:]
        (bundle / f"genome_{i}.fa").write_text(f">synthetic_{i}\n{genome}\n")
        rows.append((aid, f"synthetic_{i}", f"genome_{i}.fa", 1))
        silva.append(f">Bacteria;Synthetic;Synthetic;Synthetic;Synthetic;Synthetic;species_{i};\n{molecule}\n")
    if genome_first:
        pd.DataFrame([{"fasta":row[2], "Domain":"Bacteria", "Phylum":"Synthetic", "Class":"Synthetic", "Order":"Synthetic", "Family":"Synthetic", "Genus":"Synthetic", "Species":f"synthetic_{i}"} for i,row in enumerate(rows)]).to_csv(bundle / "taxonomy.tsv",sep="\t",index=False)
    else:
        pd.DataFrame(rows, columns=["ASV_ID", "genome_id", "fasta", "marker_copies"]).to_csv(bundle / "genomes.tsv", sep="\t", index=False)
    (destination / "markers.fa").write_text("".join(silva))
    (destination / "study.yaml").write_text("study_name: paired_smoke\ngroups:\n  - name: smoke\n    n_samples: 2\n")
    root = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load((root / "config/example.yaml").read_text())
    cfg["study_design_file"] = str(destination / "study.yaml")
    cfg["reference"].update(dada2_silva_fasta=str(destination / "markers.fa"), output_dir=str(destination / "reference"))
    cfg["simulation"].update(n_asvs=4, median_read_depth=100)
    cfg["reference_fixtures"]["mitochondria"]["source_fasta"] = str(root / cfg["reference_fixtures"]["mitochondria"]["source_fasta"])
    cfg["chemistry"].update(compounds=["compound_1"], drivers_per_compound=1, driver_min_prevalence=0)
    for key in ("group_differential_abundance", "microbial_network", "microbiome_batch_effect", "chemistry_batch_effect"):
        cfg["artifacts"][key]["enabled"] = False
    cfg["artifacts"]["contaminants"].update(n_asvs=1, prevalence=1, mean_reads=10)
    cfg["artifacts"]["chimeras"].update(n_chimeras=1)
    cfg["artifacts"]["extraction_controls"].update(labels=["blank"], contaminant_mean_reads=20)
    cfg["fastq"]["cpus"] = 1
    cfg["wgs"] = dict(enabled=True, reference_dir=str(bundle), read_pairs=100, model="miseq", cpus=1, gzip=True)
    if genome_first:
        cfg["reference"].update(source="genomes", genome_dir=str(bundle))
        cfg["wgs"].pop("reference_dir")
    (destination / "config.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
    print(destination / "config.yaml")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--genome-first", action="store_true", help="Derive V4 references from raw genomes without a prelinked ASV list")
    args=parser.parse_args()
    create(args.destination, args.genome_first)
