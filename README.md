# DECOI: **D**ata **E**mulator for **C**ommunity **O**mics tool Benchmark**I**ng

DECOI generates truth-aware V4 16S amplicons, correlated chemistry and optional
matched short-read metagenomics from one mock microbial study. It records imposed
artifacts and sequence provenance for workflow and statistical-method testing.
The chemistry is correlated by construction, not a mechanistic metabolic model.

[Full user guide](https://hallamlab-decoi.readthedocs.io/en/latest/) ·
[Reviewer test](https://hallamlab-decoi.readthedocs.io/en/latest/reviewer-test.html) ·
[Issues and feature requests](https://github.com/hallamlab/DECOI/issues)

![DECOI workflow](docs/assets/workflow-brief.svg)

## Quick start

```bash
git clone https://github.com/hallamlab/DECOI.git
cd DECOI
mamba env create -f environment.yml
conda activate decoi
Rscript scripts/install_sparsedossa2.R
```

## Small reviewer run

```bash
python scripts/create_paired_smoke.py reviewer --genome-first
nextflow run main.nf \
  --config reviewer/config.yaml --study reviewer/study.yaml \
  --genome_dir reviewer/genomes --threads 1 --outdir reviewer/results
python -m http.server 8765 --directory reviewer/results/dataset/mock_dataset
```

Open `http://localhost:8765/report.html`. The small synthetic input exercises
amplicon reads, WGS reads, chemistry, controls and their ground-truth tables.
For a full study, optional DADA2 validation, scientific limits and citations,
see the [user guide](docs/index.md). Matched WGS support is experimental.

## Citation

Cite the tools and references used in your run, listed in
[software and reference citations](docs/citations.md). A DECOI publication
citation has not yet been assigned.
