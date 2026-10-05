# Small installation test

This test uses eight tiny synthetic genomes, two biological samples and one
extraction blank. It generates matched amplicon and WGS reads, chemistry,
metadata, ground-truth tables and a report. It runs the real tools without
requiring SILVA or the full airway genome panel.

## Install

On Linux with mamba and Git available, start from a shell without an active Python
virtual environment (`.venv`). If DECOI is already installed, activate its mamba
environment and skip to the run commands.

```bash
git clone --branch refactor/user-guide-installation https://github.com/hallamlab/DECOI.git
cd DECOI
mamba env create -f environment.yml
conda activate decoi
Rscript scripts/install_sparsedossa2.R
python -m pip install .
```

Installation downloads software dependencies. The tiny reference inputs are
generated locally; the supporting Python/R environment is larger than the demo.
See [installation](installation.md) for details.

## Run

These commands work from any directory after installation:

```bash
decoi check
decoi test -o test-output
```

`check` verifies the required Python modules and external programs, including the
R packages. `test` creates its inputs under `test-output/inputs/` and requests
**one CPU and 4 GB RAM per task**. Reference preparation and study simulation
must both complete successfully. The controller prints the output and log locations.

Expected outputs:

- Three libraries per assay: two biological samples and one extraction blank.
- Twelve compressed FASTQs: R1 and R2 for both assays in each library.
- 100 WGS pairs per library (300 total); amplicon depths follow the simulation.
- One chemical measurement per biological sample, with no chemistry for the blank.
- Shared metadata, assay manifests, ground-truth tables and `report.html`.

## Review and resume

```bash
decoi report -o test-output --serve
```

Open `http://localhost:8765/report.html`. Stop the server with Ctrl-C. For a remote
server, use `--no-browser` and [SSH port forwarding](troubleshooting.md).
The dataset is in `test-output/dataset/mock_dataset/`. Both `fastq_validation.tsv`
and `wgs/fastq_validation.tsv` should have matching expected/R1/R2 counts.

```bash
decoi test -o test-output --resume
```

Completed matching Nextflow tasks should be cached. Keep `.decoi/`, including its
work directory, when testing resume. Resume operates at workflow-stage level,
not within a partially generated sample. To start independently, choose a new
output directory.

DADA2 is not part of this tiny test: its error-model training needs more reads.
Use a larger study for [optional amplicon validation](validation.md).
