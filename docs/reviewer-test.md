# Small reviewer test

This test uses eight tiny synthetic genomes, two biological samples and one
extraction blank. It generates matched amplicon and WGS read pairs, chemistry,
metadata and ground-truth tables. It exercises the real tools without downloading
SILVA or the full airway panel. It is an installation test, not a biological benchmark.

After [installation](installation.md), from the repository root:

```bash
python scripts/create_paired_smoke.py reviewer --genome-first
nextflow run main.nf \
  --config reviewer/config.yaml --study reviewer/study.yaml \
  --genome_dir reviewer/genomes --threads 1 --outdir reviewer/results
```

The generator requires a new destination directory. The default Nextflow memory
reservations are 8 GB for preparation and 32 GB for simulation; these are scheduling
requests, not measurements of this tiny test's memory use.

Inspect:

```bash
ls reviewer/results/dataset/mock_dataset
python -m http.server 8765 --directory reviewer/results/dataset/mock_dataset
```

Open `http://localhost:8765/report.html`. Stop the server with Ctrl-C.
`fastq_validation.tsv` and `wgs/fastq_validation.tsv` should have identical expected
and observed R1/R2 counts. `assay_manifest.tsv` links both assays to the same sample
IDs. Chemistry is generated for biological samples, not extraction blanks.

Repeat the same Nextflow command with `-resume` to check cache reuse. Keep the
work directory when testing resume. Optional `--run_dada2 true` exercises amplicon
validation; the tiny read depth may leave few or no retained ASVs, which is not
an installation failure. Validation does not analyze the WGS reads.

For developer regression tests:

```bash
python -m pytest -q
```

The opt-in real-tool regression suite is described in the full paired-study guide.
