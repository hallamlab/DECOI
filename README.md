# DECOI: **D**ata **E**mulator for **C**ommunity **O**mics tool Benchmark<strong><u>I</u></strong>ng v1.0.0

A Nextflow DSL2 workflow that generates a mock V4 16S amplicon study with correlated chemical measurements and preserved ground truth. This feature branch adds experimental matched short-read whole-genome shotgun (WGS) metagenomics from the same study.

## What DECOI produces

- Real, primer-defined V4 sequences extracted from SILVA 138.2 or supplied genomes.
- Reference taxonomy, source mappings, and sequence provenance.
- SparseDOSSA2 microbial abundance profiles.
- Optional group-level differential abundance.
- Optional ASV-specific PCR amplification bias.
- Optional low-abundance environmental/reagent contaminants.
- Optional two-parent chimeric amplicons.
- Sparse positive and negative ASV–compound relationships.
- Configurable, truth-recorded microbial and chemical batch effects.
- Untrimmed, primer-bearing paired-end FASTQs generated with the InSilicoSeq MiSeq model.
- Exact count checks for every R1/R2 pair.
- A machine-readable run manifest and a compact HTML report.
- A reusable tabular FASTQ manifest, paired-patient metadata, and FASTA reference fixtures for complete downstream ASPIRE testing.
- An optional DADA2 validation workflow.
- Optional paired WGS FASTQs from the same biological community, contaminants, and controls.
- Genome-copy/DNA-weighted WGS ground truth and a combined sample/assay manifest.
- A frozen airway-prioritized RefSeq panel: 547 biological genomes and three reserved reagent-contaminant representatives.

The simulator is intended for pipeline and statistical-method testing. The chemistry is correlated by construction but is not claimed to be metabolically realistic.

For setup on another machine and the full 64-sample paired test, follow
[SETUP_WGS.md](SETUP_WGS.md). See also the [reference panel](references/airway_v1/README.md),
[selection evidence/DOIs](references/airway_v1/CITATIONS.md), and [changelog](CHANGELOG.md).
The full real-genome sequencing run is pending; reference validation and the
small real-tool integration tests have passed.

## Primer-bearing raw reads

The default locus primers are:

```text
515F: GTGYCAGCMGCCGCGGTAA   19 nt
806R: GGACTACNVGGGTWTCTAAT  20 nt
```

The simulated molecule is:

```text
resolved 515F + reference-derived V4 insert + reverse-complement(resolved 806R)
```

IUPAC ambiguity codes are resolved reproducibly to A/C/G/T for each feature. R1 begins with a concrete 515F variant and R2 begins with a concrete 806R variant. These files contain the locus-specific primers, but not complete Illumina adapters, indices, EMP pads, or linkers.

## Installation

```bash
mamba env create -f environment.yml
conda activate mock16s-chem-v1
Rscript scripts/install_sparsedossa2.R
```

SparseDOSSA2 is installed separately from a pinned GitHub commit by that script;
it is not installed merely by creating the Conda environment. The environment
YAML is not a complete version lock. The [new-machine guide](SETUP_WGS.md)
includes dependency checks, tests, environment recording, and disk planning.

`main.nf` currently contains development-machine Conda paths. On another system,
use the direct Python route in that guide or its documented Nextflow process
override before running any of the Nextflow examples below.

## Amplicon-only full run

```bash
nextflow run main.nf \
  -profile standard \
  --config config/example.yaml \
  --study study/example_study.yaml \
  --outdir results
```

The first run downloads the DADA2-formatted SILVA 138.2 species training set and verifies its MD5 checksum. Later runs can reuse the Nextflow cache.

Use an existing download with:

```bash
nextflow run main.nf \
  --silva_fasta /path/to/silva_nr99_v138.2_toSpecies_trainset.fa.gz \
  --config config/example.yaml \
  --study study/example_study.yaml \
  --outdir results
```

Enable end-to-end DADA2 validation with:

```bash
nextflow run main.nf \
  --config config/example.yaml \
  --study study/example_study.yaml \
  --outdir results \
  --run_dada2 true
```

## Study request YAML

`study/example_study.yaml` defines a repeated-participant cohort layout:

```yaml
study_name: mock_airway_chemistry

cohorts:
  - name: control
    participant_prefix: CTRL
    n_participants: 10
    metadata:
      Case: Control
      disease_status: control
      site: airway
      lung_status: Healthy
    participant_metadata_cycle:
      - batch: plate_1
      - batch: plate_2
    sample_types:
      - name: Bronchial Brush
        code: BRUSH
      - name: BAL
        code: BAL

  - name: case
    participant_prefix: CASE
    n_participants: 10
    metadata:
      Case: Cancer
      disease_status: case
      site: airway
    participant_metadata_cycle:
      - batch: plate_1
      - batch: plate_2
    sample_types:
      - name: Bronchial Brush
        code: BRUSH
        metadata:
          lung_status: TumorSide
      - name: Bronchial Brush
        code: BRUSH_CONTRA
        metadata:
          lung_status: Contralateral
      - name: BAL
        code: BAL
        metadata:
          lung_status: TumorSide
      - name: BAL
        code: BAL_CONTRA
        metadata:
          lung_status: Contralateral
```

Every metadata field is copied into sample_metadata.tsv. Cohort designs also
emit Participant_ID, Case, and Type_Group. This example keeps each participant
within one batch, crosses both case groups over both batches, and gives cancer
participants paired tumor-side and contralateral bronchial-brush and BAL
samples. Excluding contralateral cancer samples leaves balanced 10-vs-10
patient-level comparisons for both sample types. Set
artifacts.group_differential_abundance.group_columns to impose known effects
for any generated metadata categories. The legacy independent groups layout
remains supported.

The current design expands to 60 biological samples. Four configured extraction
blanks give 64 sequenced libraries; chemistry is generated for biological samples
only. A study-design file overrides the fallback `simulation.n_samples` value.
The paired airway configuration uses the same design and artifact settings.

## Artifact configuration

Artifacts are explicit and independently switchable in `config/example.yaml`.
Selected current settings are shown below; the YAML also configures microbial
network modules, mitochondrial fixtures, and extraction controls:

```yaml
artifacts:
  group_differential_abundance:
    enabled: true
    group_columns: [Type_Group, Case]
    asvs_per_group: 6
    log_fold_change_sd: 0.25
    min_abs_log_fold_change: 2.50
    candidate_min_prevalence: 0.60
    disjoint_drivers: true
  microbiome_batch_effect:
    enabled: true
    column: batch
    asvs_per_batch: 10
    log_fold_change: 1.25
    direction: mixed
  pcr_bias:
    enabled: true
    log_mean: 0.0
    log_sd: 0.35
  contaminants:
    enabled: true
    n_asvs: 3
    prevalence: 0.45
    mean_reads: 40
  chimeras:
    enabled: true
    n_chimeras: 5
    fraction_of_reads: 0.02
  chemistry_batch_effect:
    enabled: true
    sd: 0.25
```

PCR bias, group effects, and microbial batch effects redistribute the original sample depth. Contaminants add reads. Chimeras transfer reads from their parent sequences so their creation does not inflate depth. Cohorts may define `participant_metadata_cycle` to cross technical variables such as batch with biological groups while keeping repeated samples from each participant together.

## Output layout

The layout below is for the SILVA/Nextflow route. For paired genome-first output,
including `wgs/` and `assay_manifest.tsv`, see [the paired output layout](SETUP_WGS.md#6-inspect-completion).
Direct Python execution writes the dataset directly into the requested output
directory, without Nextflow's `dataset/mock_dataset/` nesting.

```text
results/
├── reference/
│   ├── source/
│   │   └── silva_nr99_v138.2_toSpecies_trainset.fa.gz
│   └── v4/reference_v4/
│       ├── v4_asv_pool.fasta
│       ├── v4_asv_taxonomy.tsv
│       ├── v4_asv_source_map.tsv
│       └── reference_manifest.json
└── dataset/mock_dataset/
    ├── asv_counts_biological.tsv
    ├── asv_counts_post_pcr.tsv
    ├── asv_counts_final.tsv
    ├── asv_counts.tsv
    ├── asv_relative_abundance.tsv
    ├── asv_sequences.fasta
    ├── asv_taxonomy.tsv
    ├── chemistry.tsv
    ├── sample_metadata.tsv
    ├── ground_truth_feature_registry.tsv
    ├── ground_truth_asv_chem.tsv
    ├── ground_truth_group_effects.tsv
    ├── ground_truth_microbiome_batch_effects.tsv
    ├── ground_truth_network_modules.tsv
    ├── ground_truth_pcr_bias.tsv
    ├── ground_truth_chimeras.tsv
    ├── ground_truth_chemistry_batch.tsv
    ├── ground_truth_mitochondria.tsv
    ├── ground_truth_extraction_controls.tsv
    ├── fastq_validation.tsv
    ├── fastq_manifest.tsv
    ├── references/contaminants.fasta
    ├── references/mitochondria.fasta
    ├── ground_truth_reference_filters.tsv
    ├── manifest.json
    ├── report.html
    ├── iss_inputs/
    └── fastq/
```

## Matched short-read WGS (experimental)

Enable WGS to generate both assays for every sample in one mock study, including extraction blanks. Sample metadata, chemistry, group effects, microbial batch effects, biological ASV abundances, and contaminant identities are shared. Existing amplicon outputs stay at the dataset root; WGS outputs are under `wgs/`, with a combined `assay_manifest.tsv` linking both assays by `sample_id`.

### Genome-first reference preparation

A reusable airway-prioritized RefSeq panel is described in
[`references/airway_v1/README.md`](references/airway_v1/README.md). Its frozen
accession/checksum lock makes the genome pool reproducible across runs and
machines. Three reference genomes from reported reagent-contaminant genera are
reserved with `reference_role=contaminant`; biological community selection cannot
draw from that reserved pool. When a custom reference has no role annotations,
the original random contaminant-selection behavior remains available.
The [selection evidence and DOI table](references/airway_v1/CITATIONS.md)
separates airway ecology, reagent-contamination evidence, and benchmark-design
choices; it does not claim that every selected strain is a documented lung isolate.

For a new paired study, start with genome references and derive the V4 pool from them. The cohort remains `study/example_study.yaml`; `config/airway_paired.yaml` preserves its simulation, chemistry, artifact, and control settings and changes the reference source to genomes. The lung labels describe the study design: the inherited SparseDOSSA2 template is still `Stool`, and reference organisms are not automatically selected for lung ecology.

Place one assembly per `.fa`, `.fasta`, or `.fna` file (optionally gzip-compressed) under a source directory. DECOI discovers files recursively. A `taxonomy.tsv` at the source directory root supplies taxonomy without requiring an ASV list:

```text
fasta	Domain	Phylum	Class	Order	Family	Genus	Species
assembly_1.fna	Bacteria	...	...	...	...	...	...
```

Replace the ellipses with actual reference taxonomy. File paths are relative to the source directory. An optional `reference_role` column reserves rows as `biological` or `contaminant`; accession and panel priority metadata are retained in prepared references. Alternatively, `reference.taxonomy_reference_dir` can name an existing prepared SILVA reference directory: taxonomy is transferred only for exact V4 sequence IDs, without inferring genome identity from a taxonomic name. Source taxonomy takes precedence where supplied. Missing taxonomy is never fabricated; the configured domain and minimum-rank filters apply. Reference preparation reads local genomes; the separate panel builder downloads and freezes a public genome collection.

```bash
conda activate mock16s-chem-v1
python mock16s_chem.py --config config/airway_paired.yaml \
  --genome-dir /path/to/genome_sources \
  --reference-dir reference_airway_genomes prepare-reference
python mock16s_chem.py --config config/airway_paired.yaml \
  --reference-dir reference_airway_genomes simulate --output /path/to/paired_study_output
```

For Nextflow, supply `--genome_dir` explicitly. This stages the raw genomes and bypasses the SILVA download/preparation processes:

```bash
nextflow run main.nf --config config/airway_paired.yaml \
  --study study/example_study.yaml --genome_dir /path/to/genome_sources \
  --outdir /path/to/paired_study_output
```

Preparation searches both strands for full, exact IUPAC-compatible primer matches and retains inserts within the configured V4 length limits. It searches linear contigs; mismatched primer sites and loci spanning a circular origin are not supported. It writes `genome_selection.tsv` with an inclusion/exclusion reason for every assembly, `genome_marker_loci.tsv` with zero-based, half-open marker coordinates, the V4 FASTA/taxonomy tables, and a generated `genomes/genomes.tsv` plus genome FASTAs. The run subsequently samples its biological and contaminant features from this genome-derived pool. Assembly hashes provide stable genome IDs; marker sequence hashes provide ASV IDs. No genomic bases are synthesized or changed to force a marker match.

The initial one-to-one model excludes assemblies with multiple distinct V4 alleles and retains only the first genome in sorted source-path order when genomes share an allele. Multiple identical marker copies within a retained genome are supported and verified. These restrictions and the exact-primer requirement can substantially reduce the eligible pool. The airway configuration needs at least **503 eligible genomes** for 500 biological features plus three contaminants. Genome selection reasons, copy numbers, and reference hashes are retained for audit. This creates a new realization of the same study design; it cannot retroactively preserve unrelated SILVA ASVs from an old run.

An offline genome-first smoke test is available without any manually linked ASVs:

```bash
python scripts/create_paired_smoke.py /tmp/decoi-genome-first-demo --genome-first
python mock16s_chem.py --config /tmp/decoi-genome-first-demo/config.yaml prepare-reference
python mock16s_chem.py --config /tmp/decoi-genome-first-demo/config.yaml \
  simulate --output /tmp/decoi-genome-first-demo/output
```

### Using an existing linked reference bundle

Provide a directory containing genome FASTAs and a tab-separated `genomes.tsv`:

```text
ASV_ID	genome_id	fasta	marker_copies
ASV_<sequence_hash>	organism_1	organism_1.fa	1
```

`ASV_ID` must identify an exact V4 sequence in the prepared reference. FASTA paths are relative to this directory. Each ASV maps to one genome and each genome to one ASV in this initial implementation. `marker_copies` is the number of exact occurrences of that V4 sequence, on either strand, across the genome's contigs; DECOI verifies it. Contigs must contain at least 1,000 unambiguous A/C/G/T bases. Genome IDs permit letters, numbers, underscores, and hyphens. Keep accession/version information in an additional manifest column if useful; additional columns are retained in the output registry.

Biological and contaminant ASVs are selected only from mapped genomes that satisfy the existing taxonomy filters. Supply enough eligible genomes for `simulation.n_asvs` plus `artifacts.contaminants.n_asvs`. The reference bundle represents a deliberate restriction of the study's organism pool; enabling this restriction can change the selected community compared with an unrestricted SILVA run.

Add to your YAML configuration:

```yaml
wgs:
  enabled: true
  reference_dir: /absolute/path/to/genome_bundle
  read_pairs: 100000
  model: miseq
  cpus: 1
  gzip: true
```

For Nextflow, pass the bundle explicitly so its files are staged with the simulation task:

```bash
nextflow run main.nf \
  --config config/my_paired_study.yaml \
  --study study/example_study.yaml \
  --wgs_reference /absolute/path/to/genome_bundle \
  --outdir results/paired
```

The Python CLI also accepts `--wgs-reference /absolute/path/to/genome_bundle` before `simulate`. This argument enables WGS and overrides `wgs.reference_dir`.

WGS starts with pre-PCR biological counts plus the shared contaminant component. These counts are interpreted as marker equivalents: dividing by `marker_copies` gives cell equivalents, and multiplying by genome length gives DNA sampling weights. Read pairs are sampled at the configured WGS depth, then distributed over contigs in proportion to length. InSilicoSeq runs in `metagenomics` mode on those genome contigs without amplicon primers or read-through padding. WGS uses an independent random stream so its depth and sequencing settings do not change amplicon or chemistry generation. Blanks share their mapped contaminant/background composition; a completely empty mapped library produces empty FASTQs rather than an invented community.

WGS truth files include `genome_registry.tsv` (ASV IDs, shared feature UUIDs, copy numbers, genome lengths, FASTA checksums), `shared_marker_counts.tsv`, `cell_equivalents.tsv`, `expected_dna_fractions.tsv`, `genome_read_pairs.tsv`, and `contig_read_pairs.tsv`. These distinguish biological composition, expected DNA proportions, and realized sequencing counts. `fastq_manifest.tsv` and `fastq_validation.tsv` record files and verified paired-read totals. Chemistry remains defined against shared biological ASVs, which can be linked to genomes through the registry; it is not a genome-function-based chemical model.

This first version supports amplicon-only and paired amplicon/WGS runs. WGS-only execution, multiple genomes sharing one V4 allele, heterogeneous 16S alleles within a genome, functional annotation truth, dedicated WGS library biases, and full host genomes are not implemented. Amplicon chimeras and the short mitochondrial marker fixture are excluded from WGS and listed as such in the run manifest. WGS read depth is currently one fixed requested depth per nonempty library, including blanks. The synthetic smoke references below test software behavior and are not biological genome models.

### Adding WGS to an existing study

To add WGS to an already-completed airway study in place, use `add-wgs` after supplying a genome bundle that covers **every existing biological and contaminant ASV**:

```bash
conda activate mock16s-chem-v1
python mock16s_chem.py --config config/airway_paired.yaml \
  --wgs-reference /absolute/path/to/verified_airway_genomes \
  add-wgs --output results_control_bearing/dataset/mock_dataset
```

The cohort definition is `study/example_study.yaml`. For `add-wgs`, pass an already-linked bundle with `--wgs-reference`; existing biological counts, final contaminant/control counts, sample IDs, and the saved study seed are authoritative. Only WGS settings are taken from the new configuration. The command adds `wgs/` and `assay_manifest.tsv` and updates `manifest.json` and `report.html`. Existing amplicon FASTQs, count tables, metadata, and chemistry are preserved. Missing genome mappings fail before sequencing; a failed WGS simulation is cleaned up so it can be retried. Existing WGS output is never overwritten. This path does not rerun SparseDOSSA2 or choose a new community. Genome-first preparation cannot supply an arbitrary historical ASV unless that exact marker occurs in a supplied genome.

### Small paired-assay smoke test

With the DECOI environment active and SparseDOSSA2 installed, run from the repository root. Choose a new destination directory:

```bash
python scripts/create_paired_smoke.py /tmp/decoi-paired-demo
python mock16s_chem.py --config /tmp/decoi-paired-demo/config.yaml prepare-reference
python mock16s_chem.py --config /tmp/decoi-paired-demo/config.yaml \
  simulate --output /tmp/decoi-paired-demo/output
```

This creates eight small synthetic genomes with matching marker references, two biological samples, one extraction blank, correlated chemistry, and both sets of paired FASTQs. Inspect `/tmp/decoi-paired-demo/output/assay_manifest.tsv`; WGS requests 100 pairs per sample. No reference download is needed.

Unit tests run with `python -m pytest -q` in an environment containing pytest and the Python dependencies. To also exercise real Cutadapt, SparseDOSSA2, and InSilicoSeq, run `DECOI_RUN_INTEGRATION=1 python -m pytest -q`. The integration test compares paired-run amplicon counts, chemistry, metadata, feature registries, and decompressed amplicon FASTQs with an amplicon-only run using the same eligible reference pool.

## Ground-truth interpretation

- `asv_counts_biological.tsv`: shared biological counts after group effects, microbial network modules, and microbial batch effects, before PCR bias and injected features.
- `asv_counts_post_pcr.tsv`: counts after ASV-specific PCR efficiency.
- `asv_counts_final.tsv`: exact per-feature read-pair requests supplied to InSilicoSeq, including contaminants and chimeras.
- `ground_truth_feature_registry.tsv`: stable feature UUID, ASV ID, feature type, exact V4 sequence, sequence hash, representative SILVA record, taxonomy, and original SparseDOSSA2 feature mapping.
- `ground_truth_asv_chem.tsv`: exact nonzero ASV–compound coefficients.
- `ground_truth_group_effects.tsv`: imposed group-specific log fold changes.
- `ground_truth_network_modules.tsv`: feature membership and imposed microbial network module effects.
- `ground_truth_microbiome_batch_effects.tsv`: ASVs receiving known batch-specific effects, including batch, signed log fold-change, and direction.
- `ground_truth_pcr_bias.tsv`: per-ASV PCR efficiencies.
- `ground_truth_chimeras.tsv`: chimera parents and breakpoint.
- `ground_truth_chemistry_batch.tsv`: compound-specific batch effects.
- `ground_truth_mitochondria.tsv`: mitochondrial source record, exact sequence checksum, configured abundance, observed prevalence, and injected read count. By default DECOI uses a vendored, checksum-pinned segment of the human RefSeq mitochondrial genome, avoiding a runtime network dependency; `reference_fixtures.mitochondria.source_fasta` can select another local fixture and `source_url` remains available explicitly. DECOI synthetically adds the configured amplicon primers to this genuine mitochondrial template so downstream mitochondrial filtering can be tested; it does not assert natural amplification by those primers.
- `ground_truth_extraction_controls.tsv`: synthetic extraction-blank identifiers, read totals, and DNA concentrations. When `artifacts.extraction_controls.enabled` is true, these controls are enriched for the implanted contaminant ASVs and are included in the FASTQ manifest and metadata so prevalence- and frequency-based decontamination workflows can be tested without patient data.
- `fastq_validation.tsv`: expected and observed R1/R2 record counts.

## Reproducibility

`manifest.json` stores:

- schema and software version;
- study name;
- full resolved configuration;
- random seed;
- artifact settings;
- dimensions and total read pairs;
- Python, R, Cutadapt, and InSilicoSeq version strings;
- reference preparation metadata (SILVA or genome-derived);
- paired-run assay summaries when WGS is enabled.

Feature UUIDs use deterministic UUIDv5 values and remain stable for a given ASV identifier.

For genome-derived runs, prepared reference tables retain source checksums,
marker loci/copy numbers, taxonomy, and supplied accession/role metadata. The
frozen panel's exact inputs are pinned separately in `references/airway_v1/panel.lock.json`.
Archive its hash and the complete environment with each run. Stable references
do not guarantee byte-identical output across different simulator/library versions.

## DADA2 validation

The optional validation stage:

1. removes 515F and 806R with Cutadapt;
2. filters paired reads;
3. learns forward and reverse error models;
4. denoises each direction;
5. merges pairs;
6. removes bimeras;
7. writes `dada2_sequence_table.tsv` and `dada2_tracking.tsv`.

The supplied settings are reasonable defaults, not universal optimums. Adjust truncation lengths for the selected InSilicoSeq MiSeq model when needed.
This stage validates amplicon reads only; it is not a WGS classifier or assembly
validation. The Python `simulate` command does not launch this separate stage.

## Tests

```bash
pytest -q
```

The unit suite tests stable identifiers, primer resolution, study expansion and participant metadata cycling, depth-preserving group and microbial batch effects, contaminants, extraction controls, mitochondria, chimeras, and chemistry batch effects without requiring SILVA or a full R/InSilicoSeq installation. It also tests genome-derived marker extraction, unique mappings/copy counts, role separation, panel checksums, and WGS allocation. Set `DECOI_RUN_INTEGRATION=1` to include the real-tool paired tests. The development run passed 29 tests with this setting; that is not a full-size real-genome study benchmark.

## Current scope

The released v1.0 supports SILVA-derived V4 amplicons and InSilicoSeq MiSeq paired-end output. This feature branch additionally supports experimental paired WGS using explicitly mapped genome references, as described above. It does not yet model KEGG pathways, mechanistic metabolism, index hopping, complete adapter/index constructs, PacBio, or Nanopore reads.

## Citing DECOI and its dependencies

If you use DECOI in a publication, please cite DECOI itself once a project citation is available and cite the software and reference-data publications relevant to your run. Nextflow applies to orchestrated runs; Cutadapt applies to SILVA preparation or DADA2 primer removal. Genome-first Python runs do not use SILVA or Cutadapt for reference preparation. Cite DADA2 only when the optional validation stage is enabled with `--run_dada2 true`.

### Core workflow software

- **Nextflow:** Di Tommaso, P. *et al.* (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35, 316–319. [https://doi.org/10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)
- **SparseDOSSA2:** Ma, S. *et al.* (2021). A statistical model for describing and simulating microbial community profiles. *PLOS Computational Biology*, 17(9), e1008913. [https://doi.org/10.1371/journal.pcbi.1008913](https://doi.org/10.1371/journal.pcbi.1008913)
- **Cutadapt:** Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*, 17(1), 10–12. [https://doi.org/10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200)
- **InSilicoSeq:** Gourlé, H., Karlsson-Lindsjö, O., Hayer, J., and Bongcam-Rudloff, E. (2019). Simulating Illumina metagenomic data with InSilicoSeq. *Bioinformatics*, 35(3), 521–522. [https://doi.org/10.1093/bioinformatics/bty630](https://doi.org/10.1093/bioinformatics/bty630)
- **NumPy:** Harris, C. R. *et al.* (2020). Array programming with NumPy. *Nature*, 585, 357–362. [https://doi.org/10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)
- **pandas:** McKinney, W. (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*, 56–61. [https://doi.org/10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a)
- **Biopython:** Cock, P. J. A. *et al.* (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422–1423. [https://doi.org/10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)

### Optional validation software

- **DADA2** (only when `--run_dada2 true`): Callahan, B. J. *et al.* (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13, 581–583. [https://doi.org/10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)

### Reference data

For the genome-first airway panel, report the exact RefSeq assembly accessions
and lock SHA-256 and cite the applicable [airway/reagent evidence](references/airway_v1/CITATIONS.md).
The SILVA citations below apply only when that reference is used; they do not
describe the genome-first panel.

- **SILVA:** Quast, C. *et al.* (2013). The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. *Nucleic Acids Research*, 41(D1), D590–D596. [https://doi.org/10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)
- **Exact default training set:** SILVA 138.2 NR99 taxonomic training data formatted for DADA2 (`silva_nr99_v138.2_toSpecies_trainset.fa.gz`). [https://doi.org/10.5281/zenodo.14169026](https://doi.org/10.5281/zenodo.14169026)

### Dependency-audit note

This list was derived from the workflow definitions and source imports, not solely from the Conda environment. Python, R, PyYAML, `optparse`, `curl`, and GNU `md5sum` are directly used but do not have a single canonical peer-reviewed software DOI that DECOI can recommend. SciPy, SeqKit, pigz, wget, `jsonlite`, `remotes`, and the listed SparseDOSSA2 support packages are present in `environment.yml` for installation or transitive support but are not directly invoked by the current standard workflow, so they are not included as required citations above.
