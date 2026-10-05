# Paired airway study: setup and full test on a new machine

This guide applies to `feature/wgs-short-read-metagenomics`. It generates V4
amplicons, short-read WGS, and chemistry for **one shared mock study**, not two
independently sampled communities. Run commands from the repository root in Bash.
The full real-genome study has not yet been completed; this is the procedure for
that test, not a claim that it has passed.

## 1. Obtain the branch and install the environment

After the branch has been pushed:

```bash
git clone --branch feature/wgs-short-read-metagenomics https://github.com/hallamlab/DECOI.git
cd DECOI
mamba env create -f environment.yml
conda activate mock16s-chem-v1
Rscript scripts/install_sparsedossa2.R
python -m pip install pytest
```

The R installer installs SparseDOSSA2 from `biobakery/SparseDOSSA2@26a998a`
when it is absent. It does not replace an existing installation, so use a fresh
environment when reproducing that dependency. `environment.yml` is an environment
specification, **not an exact package lock**: a new solve can install different
package versions. Record the environment used for the larger run.

```bash
python --version
iss --version
cutadapt --version
Rscript -e 'stopifnot(requireNamespace("SparseDOSSA2", quietly=TRUE)); print(packageVersion("SparseDOSSA2"))'
python -c 'import Bio, numpy, pandas, yaml; print("Python imports OK")'
```

Use an adequately provisioned Linux workstation or a scheduler allocation. The
configuration requests eight InSilicoSeq workers. As conservative planning
figures, allow **32 GB RAM and at least 200 GB free disk for the direct Python
route**; these are not measured peak requirements. Allow more storage for
Nextflow, which retains work files and publishes copies of outputs.

The selected 503 genomes contain about 1.76 billion bases. WGS retains a separate
input FASTA per sample containing contigs allocated reads. If all selected
genomes were represented in every library, those inputs alone would approach
113 GB before FASTA formatting, FASTQs, and simulator temporary files. Actual
space depends on sparsity. `fastq.keep_iss_inputs` does **not** control WGS input
retention; WGS currently has no automatic cleanup setting. Monitor free space.

## 2. Run the small integration test first

```bash
DECOI_RUN_INTEGRATION=1 python -m pytest -q
```

The development environment passed all 29 tests, including real-tool paired
smokes. Those tests use small synthetic genomes, not the full frozen panel.
Without `DECOI_RUN_INTEGRATION=1`, the real-tool integration tests are skipped.

## 3. Download the exact frozen panel

Git contains the accession/checksum lock and small metadata files, **not the
genome sequences**. Rehydrate v1; do not run `build` to obtain a different live
selection.

```bash
python scripts/build_genome_panel.py rehydrate \
  --lock references/airway_v1/panel.lock.json \
  --output data/reference_panels/airway_v1/raw
python scripts/build_genome_panel.py verify \
  --lock references/airway_v1/panel.lock.json \
  --output data/reference_panels/airway_v1/raw
```

Rehydration needs access to NCBI HTTPS download endpoints. It checks the pinned
SHA-256 values and the regenerated taxonomy table. An interrupted download can be
retried with the same command; a completed file with a wrong checksum is an
error, not silently accepted. If an upstream pinned file becomes unavailable,
transfer the original `raw/` directory from the development machine and run
`verify`—do not substitute a newer assembly.

The panel contains 550 unique bacterial genome/V4 pairs: 547 biological references
and three reserved representatives of reagent-contaminant genera. No
background-tier genomes were needed. See the [panel policy and provenance](references/airway_v1/README.md)
and [supporting publications/DOIs](references/airway_v1/CITATIONS.md). This is an
airway-prioritized benchmark pool, not 550 proven human lung isolates.

## 4. Prepare and audit the paired reference

```bash
python mock16s_chem.py --config config/airway_paired.yaml prepare-reference
python scripts/audit_genome_panel.py \
  --lock references/airway_v1/panel.lock.json \
  --raw data/reference_panels/airway_v1/raw \
  --prepared data/reference_panels/airway_v1/prepared \
  --report data/reference_panels/airway_v1/validation.json
```

Preparation requires a fresh output directory. If this prepared directory already
exists and was successfully built, run the audit instead of preparing it again.
The expected result is 550 retained genome/marker pairs and `status: passed`.
The audit checks source hashes, exact V4 matches, copy numbers, genome lengths,
unique mappings, and preservation of biological/contaminant roles.

## 5. Run the full shared study

`config/airway_paired.yaml` points at the prepared panel. Its
`study_design_file: study/example_study.yaml` expands to 60 biological samples
from 20 participants; the YAML's fallback `simulation.n_samples: 30` is overridden
by this study design. Four extraction controls bring the total to 64.

| Setting | Full-test configuration |
| --- | --- |
| Biological community | 500 features sampled from the 547 biological references |
| Contaminants | Three reserved Bradyrhizobium, Ralstonia, and Sphingomonas representatives |
| Other amplicon features | Five chimeras and one mitochondrial marker fixture |
| Chemistry | Six compounds for the 60 biological samples; no invented chemistry for blanks |
| Amplicon depth | SparseDOSSA2 median 30,000 pairs, modified by configured artifacts/controls |
| WGS depth | 100,000 pairs per nonempty library; 6.4 million pairs expected across 64 libraries |
| Sequencing | InSilicoSeq MiSeq model, gzip output, eight workers per assay |
| Seed / abundance template | 42 / SparseDOSSA2 `Stool` (not calibrated to lung ecology) |

Choose a **fresh** output destination. Do not point this command at a historical
SILVA study: genome-first simulation is a new realization of the same study
design, not an in-place preservation of its old ASVs or chemistry.

```bash
mkdir -p results
set -o pipefail
test ! -e results/airway_paired_full && \
python mock16s_chem.py --config config/airway_paired.yaml \
  simulate --output results/airway_paired_full \
  2>&1 | tee results/airway_paired_full.log
```

The existence check prevents this example from reusing an output directory.
Use a persistent terminal or
batch allocation; this is not a background service. The direct Python command has
no resume/checkpoint mode. On failure, retain the log and partial output for
diagnosis; do not blindly rerun into it or delete unrelated reference/results
directories. Use a fresh output path for a new attempt.

Save environment provenance alongside this run:

```bash
conda list --explicit > results/airway_paired_full.conda-explicit.txt
python -m pip freeze > results/airway_paired_full.pip-freeze.txt
Rscript -e 'sessionInfo(); print(packageDescription("SparseDOSSA2")[c("Version","RemoteSha")])' \
  > results/airway_paired_full.R-session.txt
git rev-parse HEAD > results/airway_paired_full.git-commit.txt
git diff > results/airway_paired_full.source.patch
sha256sum references/airway_v1/panel.lock.json > results/airway_paired_full.panel-sha256.txt
```

Use a clean, committed checkout for the test: `git diff` does not capture untracked
source files. `manifest.json` records resolved settings and selected tool versions,
but does not replace these full environment records.

## 6. Inspect completion

For direct Python execution, the output root is `results/airway_paired_full/`:

```text
airway_paired_full/
├── sample_metadata.tsv          # 64 rows, shared sample IDs
├── chemistry.tsv                # 60 biological sample rows
├── asv_counts_biological.tsv    # shared pre-PCR biological community
├── asv_counts_final.tsv         # amplicon requests including artifacts/controls
├── ground_truth_*.tsv
├── fastq/                      # amplicon R1/R2
├── fastq_validation.tsv
├── fastq_manifest.tsv
├── assay_manifest.tsv           # 128 sample/assay rows
├── manifest.json
├── report.html
└── wgs/
    ├── fastq/                  # WGS R1/R2
    ├── iss_inputs/             # retained per-sample genomes and read requests
    ├── genome_registry.tsv     # 503 genomes linked to ASV IDs and feature UUIDs
    ├── shared_marker_counts.tsv
    ├── cell_equivalents.tsv
    ├── expected_dna_fractions.tsv
    ├── genome_read_pairs.tsv
    ├── contig_read_pairs.tsv
    ├── fastq_validation.tsv
    └── fastq_manifest.tsv
```

Check that the process exits successfully, `manifest.json` and `report.html`
exist, both `fastq_validation.tsv` tables have matching expected/R1/R2 counts,
and each sample occurs once per assay in `assay_manifest.tsv`. DECOI counts the
actual FASTQ records during generation and fails on discrepancies. Amplicon and
WGS abundances need not be equal: WGS converts shared marker counts to cell
equivalents using copy number, then weights by genome length. Chimeras and the
short mitochondrial fixture are amplicon-only. DADA2 validation is amplicon-only
and is not enabled in this full-test configuration.

## Optional: Nextflow on the new machine

The Python route above avoids the machine-specific Conda directives currently
present in `main.nf`. To use Nextflow, supply a local override file outside the
repository (for example `/scratch/decoi-machine.config`) with the following
content, replacing the environment path with the new machine's actual path:

```groovy
process {
    withName: '.*' {
        conda = '/absolute/path/to/envs/mock16s-chem-v1'
    }
    withName: SIMULATE_STUDY {
        cpus = 8
        memory = '32 GB'
    }
    withName: PREPARE_GENOMES {
        memory = '8 GB'
    }
}
```

Named process selectors override inline directives; see the official
[Nextflow configuration precedence](https://docs.seqera.io/nextflow/config#selector-priority).
The environment must already contain the installed SparseDOSSA2 package.

```bash
nextflow -c /scratch/decoi-machine.config run main.nf \
  -profile standard \
  --config config/airway_paired.yaml \
  --study study/example_study.yaml \
  --genome_dir data/reference_panels/airway_v1/raw \
  --outdir results/airway_paired_nextflow
```

Pass `--genome_dir` explicitly: the Nextflow route is selected by this parameter,
not by the YAML reference-source setting alone. It bypasses SILVA and prepares
the raw genomes inside its work directory, so step 4's prepared reference is not
reused by this route. The dataset is published under
`results/airway_paired_nextflow/dataset/mock_dataset/`. Preserve `work/` and
`.nextflow/` to allow `-resume` to reuse completed tasks; the simulation is one
task, so a failure partway through it does not resume at a particular sample.

For Slurm, use `-profile slurm` inside the site's supported submission workflow
and add its queue/time/resource settings to the override. Keep allocated CPUs
consistent with both `fastq.cpus` and `wgs.cpus`; these simulator settings are not
automatically derived from Nextflow's `task.cpus`. Cluster submission and memory
peaks for the full study have not yet been validated.

## Custom genomes and existing studies

The frozen panel is optional. Users may supply one local assembly per FASTA plus
`taxonomy.tsv`, or a prelinked `genomes.tsv` bundle. A list of accessions alone is
not accepted by the simulation CLI; first download the assemblies and provide
their taxonomy. See [the README's genome-first and linked-bundle sections](README.md#matched-short-read-wgs-experimental)
for schemas, copy-number checks, and allele restrictions.

`add-wgs` can preserve an existing study only when the supplied linked genomes
cover every original biological and contaminant ASV by exact sequence. Neither
genus names nor this new panel establish those matches automatically.

## Files that travel through Git

Commit the source, configuration, tests, public guides, and small
`references/airway_v1/` files (including lock, accession table, citations, and
validation report). Large downloaded genomes, generated results, Nextflow work,
and environment installations remain local. `docs/` remains ignored for the
unfinished application-note methods/results. Nothing in this guide requires
publishing those private drafts or committing generated FASTQs.
