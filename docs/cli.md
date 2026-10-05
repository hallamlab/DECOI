# Command-line guide

DECOI uses Nextflow to run the workflow, while the `decoi` controller handles
input paths, resources, logs and output locations. Install it in the supporting
mamba environment, then use it from any working directory.

```bash
decoi --help
decoi run --help
decoi check
decoi test -o test-output
```

## Run a study

```bash
decoi run --config study-config.yaml -o study-output --threads 8
```

The configuration defines the scientific study and which assays are enabled.
`study_design_file` selects the study design; `reference.source: genomes` and
`reference.genome_dir` select genome-first preparation. A SILVA-based configuration
can supply `reference.dada2_silva_fasta`, or let the workflow download its pinned
SILVA reference. Paired WGS requires genome-first inputs or a prelinked
`wgs.reference_dir` bundle with `wgs.enabled: true`.

Relative paths inside a custom YAML are resolved against that YAML's directory.
The bundled example YAMLs retain their resource-root-relative layout. Paths
passed as CLI flags are resolved against your current directory. You can override
inputs with `--study`, `--genome_dir`, `--silva_fasta` and `--wgs_reference`.
Raw genomes and prelinked WGS references cannot be selected together.

Use `--run_dada2` or `validation.run_dada2: true` for downstream amplicon validation.
The controller forwards the configured primers, truncation lengths and error
thresholds. This is separate from the minimal installation test.

## Resources and Slurm

`run` defaults to eight simulation threads and 32 GB per task. `test` defaults
to one thread and 4 GB. Reference preparation uses one CPU; simulation and optional
DADA2 receive `--threads`. On a local machine, threads are capped by `--max_cpus`
(or available CPUs when omitted). `--max_memory` optionally limits the aggregate
local memory budget; it must accommodate the per-task `--memory` request.

```bash
decoi run --config study-config.yaml -o study-output \
  --threads 8 --memory '32 GB' --max_cpus 16 --max_memory '64 GB'
```

From a cluster login node or suitable controller allocation:

```bash
decoi run --config /shared/study-config.yaml -o /shared/study-output \
  --executor slurm --account YOUR_ACCOUNT --partition YOUR_PARTITION \
  --threads 8 --memory '32 GB' --max_tasks 2 --time_limit 24h
```

Partition can be omitted to use the cluster default. Slurm jobs request one node.
`--max_tasks` limits submitted/running jobs; it defaults to one. Local aggregate
CPU/memory flags are not used for Slurm: requests are per job. The environment,
inputs, installed DECOI resources, outputs and work directory must be visible on
compute nodes. Keep the controller running for the duration of the workflow.

A study simulation is currently **one task**, with read generation inside it.
Raising `--max_tasks` does not distribute samples across cluster nodes. Use
`--threads` to control the simulator workers. Full-scale Slurm execution remains
unvalidated; the controller generates Slurm configuration, but local tests cannot
verify a particular cluster's scheduler policies.

Local resource requests guide scheduling; they do not impose a hard operating-system
CPU limit. The tiny test trace showed usage above its one-CPU request despite the
simulator receiving one thread. Additional library threading/accounting remains
under investigation. Do not interpret the requested CPUs as a verified usage cap.

## Logs, outputs and resume

A run produces:

```text
study-output/
├── dataset/mock_dataset/    # FASTQs, tables, manifest and report.html
├── reference/               # prepared reference outputs
├── validation/              # optional DADA2 results
├── logs/<run-id>/           # console, Nextflow log, trace, timeline and execution report
└── .decoi/                  # resolved settings, Nextflow cache and work directory
```

Each launch saves its resolved settings and resource configuration with its logs.
Use a fresh output directory for a new run. To reuse matching completed stages:

```bash
decoi run --config study-config.yaml -o study-output --resume
```

Keep `.decoi/` and the work directory. `--work_dir` can place work elsewhere, but
resume must use the same location. Changing inputs/settings can invalidate cached
stages. Resume does not checkpoint within simulation of an individual sample.
Only one controller can use an output directory at a time.

## Standalone operations

```bash
decoi prepare-reference --config study-config.yaml -o prepared-reference
decoi add-wgs --config linked-wgs.yaml -o existing-dataset --threads 8
```

Reference preparation requires a fresh destination and runs directly in the active
environment. `add-wgs` requires an existing dataset containing `manifest.json`,
and exact sequence-linked genome coverage of its biological and contaminant ASVs.
For a controller run, pass `study-output/dataset/mock_dataset` as that dataset.
These operations save logs but do not use Nextflow checkpoints or Slurm flags.
They must not run concurrently with another operation on the same dataset.

## Open results

```bash
decoi report -o study-output --serve
```

Open `http://localhost:8765/report.html`. `report` also accepts a direct dataset
path. Without `--serve` it prints the report location. Use `--no-browser` remotely,
`--port` to select another port, and SSH forwarding to view it locally. The server
binds to localhost only; Ctrl-C stops it.
