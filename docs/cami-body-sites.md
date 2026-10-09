# Larger CAMI-derived demonstration for ASPIRE

This example starts with CAMI2 reference genomes and community compositions to
simulate an airway–oral comparison, with skin communities as biological controls
and separate extraction blanks. It produces V4 paired amplicon reads and synthetic
chemical measurements. There is no shotgun simulation, assembly or binning.
The study is CAMI2-derived, not an official CAMI dataset or measured human cohort
([Meyer et al., 2022](#references)).

## Study design

The primary showcase has 50 synthetic patients and 179 libraries:

| Cohort | Patients | Samples per patient |
| --- | --- | --- |
| Control | 25 | Healthy airway, oral and skin |
| Cancer | 25 | Tumor-side airway, contralateral airway, oral and skin |
| Extraction blanks | — | Four additional libraries |

This gives 75 airway, 50 oral and 50 skin libraries, plus four blanks. Skin is
used as the BIO control class for ASPIRE prevalence-based decontamination;
Airways and Oral form the biological analysis cohort. The clinical labels do not describe CAMI donors; they are constructed
so patient, disease-group and lung-side analysis paths can be tested.

For each site, ten original CAMI compositions are reused across the control and
cancer cohorts. A seeded per-patient ASV log effect (normal, SD 0.35) is shared
across that patient's samples. Independent per-library log effects (SD 0.2) add
within-patient variation; weights are renormalized and source zeros stay zero.
Tumor-side and contralateral airway samples share a source composition. Batch
is assigned per patient, with 13 and 12 patients per cohort in the two batches.
`reference/patient_log_effects.tsv` records the shared effects. Metadata records
the source composition, synthetic patient, disease group and lung side.

Body-site, cancer/control and lung-side abundance effects are subsequently
injected and recorded in `ground_truth_group_effects.tsv`. Lung-side effects apply
only to airway samples. These are simulation effects, not cancer biomarkers.
The earlier 64-library body-site-only recipe remains available as
`study/cami_airway_oral_skin.yaml`; the primary recipe is
`study/cami_patient_showcase.yaml`.

The four blanks have different expected contamination loads. The expected pairs
per contaminant ASV are 200, 600, 4,000 and 8,000, respectively. With three
contaminant ASVs this spans roughly 600 to 24,000 pairs, plus small amounts of
background. Poisson draws use the study seed, so rerunning the same configuration
reproduces the counts. These values deliberately exercise both low-depth controls
and deeper contaminated libraries; they are not fitted laboratory data.
Blank depth does not guarantee passing QC after read processing.

The recipe enables all applicable amplicon effects: differential abundance by
body site, latent microbial co-variation, microbial and chemical batch effects,
PCR bias, contaminants, extraction blanks, mitochondrial reads and chimeras.
Group/network effects are deliberately moderate so they perturb rather than
replace the starting site communities. All parameters are explicit in the recipe.
DADA2 validation is enabled. SparseDOSSA does not replace the supplied CAMI
compositions, and WGS-specific modules are disabled.

Chemistry includes acetone, ethanol, 2,3-butanedione, dimethyl sulfide, acetate
and lactate, with known ASV drivers, noise, zeros and batch effects. It is generated
from the biological communities before PCR and contamination. Values are provided
for the 175 biological libraries, not fabricated for extraction blanks. These are
constructed associations for testing, not predictions of microbial metabolism.

## Install and obtain the reference inputs

Follow the [DECOI installation guide](installation.md), including environment
creation and `python -m pip install .`. These preparation scripts are shipped on
`refactor/user-guide-installation`; use that branch for this example. For an
existing installation, update from its checkout:

```bash
cd ~/repos/DECOI
git switch refactor/user-guide-installation
git pull --ff-only
conda activate decoi
python -m pip install .
decoi check
```

The input is a local copy of the original CAMI II human-microbiome reference
collection ([dataset DOI](https://doi.org/10.4126/FRL01-006425518)). Obtain the
original genome FASTAs, genome-ID maps, taxonomy tables and short-read abundance
tables, then arrange them as below. The existing `mock_2022` collection uses this
layout. Simulated metagenomic FASTQs are not used for this amplicon recipe.
The large reference collection is not bundled with DECOI, and preparation does
not download missing references.

```text
~/data/mock_2022/
  MetaGs/CAMI_II_<site>/short_read/abundance*.tsv
  SAGs/CAMI_II_<site>/genome_to_id.tsv
  SAGs/CAMI_II_<site>/genome_taxa_info.tsv
  SAGs/CAMI_II_<site>/fasta/<original genome FASTAs>
```

`<site>` must include `Airways`, `Gastrointestinal`, `Oral`, `Skin` and
`Urogenital`. All five collections are indexed for provenance; only airway,
oral and skin compositions enter this study. The `SAGs` directory must contain
original reference genomes, not simulated SAG fragments. The genome basenames
must match `genome_to_id.tsv`. Preparation checks these inputs before writing
its output and records hashes and reference-retention audits.

## Prepare the 50-patient study

From the installed DECOI source checkout:

```bash
conda activate decoi
cd ~/repos/DECOI
CAMI_SOURCE="$HOME/data/mock_2022"
DEMO_ROOT="$HOME/data/aspire-demo"
mkdir -p "$DEMO_ROOT"

python scripts/prepare_cami_demo.py \
  --input "$CAMI_SOURCE" \
  --output "$DEMO_ROOT/aspire-cami-mock-inputs" \
  --recipe study/cami_patient_showcase.yaml \
  --threads 4
```

Use a fresh preparation directory. This command extracts V4 alleles, constructs
patient-linked mixtures and writes `config.yaml`, `study.yaml`, reference audits
and `aspire_metadata.yaml`; it does not generate reads. Three reference-derived
ASVs absent from the biological cohort become recorded synthetic contaminant
spike-ins. Their role in this simulation does not classify them as established
reagent contaminants.

## Generate reads, chemistry and truth tables

```bash
DEMO_ROOT="$HOME/data/aspire-demo"
decoi run \
  --config "$DEMO_ROOT/aspire-cami-mock-inputs/config.yaml" \
  --output "$DEMO_ROOT/aspire-cami-mock" --threads 8
```

The recipe starts with 50,000 read pairs per non-blank library, applies the
configured artifacts and MiSeq error model, and runs DADA2 validation. It creates
**179 libraries**: 125 Airways/Oral biological samples, 50 Skin biological controls
and four technical extraction blanks. No assembly, binning or shotgun simulation
is required. This is a larger demonstration, so read generation and downstream
analyses take substantially longer than the small installation test.

To resume an interrupted DECOI run, retain its `.decoi` directory and run:

```bash
decoi run \
  --config "$HOME/data/aspire-demo/aspire-cami-mock-inputs/config.yaml" \
  --output "$HOME/data/aspire-demo/aspire-cami-mock" --threads 8 --resume
```

Check the output and view the report:

```bash
DEMO_ROOT="$HOME/data/aspire-demo"
test -s "$DEMO_ROOT/aspire-cami-mock/dataset/mock_dataset/fastq_manifest.tsv"
test -s "$DEMO_ROOT/aspire-cami-mock/dataset/mock_dataset/sample_metadata.tsv"
test -s "$DEMO_ROOT/aspire-cami-mock/dataset/mock_dataset/chemistry.tsv"
decoi report -o "$DEMO_ROOT/aspire-cami-mock" --serve --no-browser
```

See [report viewing and remote access](installation.md) for remote access. Stop the
report server with Ctrl-C before continuing in the same terminal.

## Run and validate the larger ASPIRE demonstration

The handoff directory is:
`~/data/aspire-demo/aspire-cami-mock/dataset/mock_dataset/`.
It contains paired FASTQs, a manifest, metadata, chemistry, contaminant and
mitochondrial reference FASTAs, and ground-truth tables. Preserve the whole
bundle so the ASPIRE configurator and validator can inspect it.

Follow the [ASPIRE larger demonstration guide](https://hallamlab-aspire.readthedocs.io/en/latest/cami-mock.html)
for installation, configuration and the recorded benchmark. After installing the
current ASPIRE checkout, the complete handoff is:

```bash
cd ~/repos/ASPIRE
DEMO_ROOT="$HOME/data/aspire-demo"
DATASET="$DEMO_ROOT/aspire-cami-mock/dataset/mock_dataset"
RESULTS="$DEMO_ROOT/aspire-cami-mock-results"
CONFIG="$DEMO_ROOT/aspire-cami-mock.yml"

./examples/configure_mock_run.sh \
  --dataset "$DATASET" --output "$RESULTS" \
  --config-out "$CONFIG" --threads 32
./run_asv_pipeline.sh "$CONFIG"
./examples/validate_mock_run.sh --dataset "$DATASET" --results "$RESULTS"
```

The generated ASPIRE YAML maps `Type_Group` values explicitly: Airways/Oral are
biological, Skin is BIO and Control is TECH. Biological samples must meet the
5,000-read post-QC inclusion threshold; nonzero controls are retained at any
depth. Independent TECH and BIO decontam prevalence tests use a score threshold
of 0.1, and their contaminant-ASV union is removed before reference screening
and final microbial filtering. The final ASV gates require at least 0.1% relative
abundance in one retained biological sample and nonzero counts in at least 5%
of those samples. Cutadapt primer removal precedes fastp with fixed clipping zero.

The supplied reference FASTAs identify the sequences injected in this synthetic
study. They are filtering fixtures, not general contaminant databases.
`aspire_metadata.yaml` is a partial settings fragment for custom integrations;
the ASPIRE configurator above produces the complete runnable demonstration YAML.
Use the ASPIRE validator for truth recovery and workflow-output checks; DECOI's
DADA2 validation separately checks recovery of simulated amplicon sequences.

## Reference mapping and truth tables

Exact IUPAC-compatible V4 primer matches are sought on both genome strands.
All eligible loci contribute, including multiple alleles. Identical sequences
collapse into one ASV; expected mixture weights sum genome abundance multiplied
by allele copy count. Genomes without eligible loci are excluded and recorded.
Shared ASVs receive their common taxonomy prefix, with conflicting lower ranks
left blank. Duplicate accession copies are checked by hash.

Keep the prepared input bundle with the run:

- `preparation.yaml` and `study.yaml`: the recipe, source compositions and roles.
- `reference/genome_selection.tsv` and `genome_marker_loci.tsv`: genome hashes, exclusions and allele coordinates.
- `reference/sample_reference_retention.tsv`: source abundance retained after marker extraction.
- `reference/source_genome_abundances.tsv` and `source_amplicon_abundances.tsv`: unmodified source mixtures.
- `reference/amplicon_abundances.tsv`: the derived library mixtures.
- `reference/synthetic_contaminants.tsv` and `study_derivation.json`: contaminant identities and derivation parameters.

The run publishes reads, count tables before and after artifacts, metadata,
chemistry, reference fixtures, a report, and separate `ground_truth_*.tsv` tables
for every enabled effect. These effects exercise analysis tools; enabling them
does not establish that their distributions match a particular clinical study.

## References

Meyer, F., Fritz, A., Deng, Z.-L., Koslicki, D., Lesker, T. R., Gurevich, A.,
Robertson, G., Alser, M., Antipov, D., Beghini, F., et al. (2022).
**Critical Assessment of Metagenome Interpretation: the second round of challenges.**
*Nature Methods* **19**, 429–440.
[DOI: 10.1038/s41592-022-01431-4](https://doi.org/10.1038/s41592-022-01431-4).
[CAMI project website](https://cami-challenge.org/).
[CAMI II human microbiome dataset, DOI: 10.4126/FRL01-006425518](https://doi.org/10.4126/FRL01-006425518).
