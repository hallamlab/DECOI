# DECOI: **D**ata **E**mulator for **C**ommunity **O**mics tool Benchmark<strong><u>I</u></strong>ng v1.0.0

A Nextflow DSL2 workflow that generates a complete mock V4 16S amplicon study with correlated chemical measurements and preserved ground truth.

## What DECOI produces

- Real, primer-defined V4 sequences extracted from SILVA 138.2.
- SILVA taxonomy and source-record mappings.
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

The simulator is intended for pipeline and statistical-method testing. The chemistry is correlated by construction but is not claimed to be metabolically realistic.

## Primer-bearing raw reads

The default locus primers are:

```text
515F: GTGYCAGCMGCCGCGGTAA   19 nt
806R: GGACTACNVGGGTWTCTAAT  20 nt
```

The simulated molecule is:

```text
resolved 515F + SILVA-derived V4 insert + reverse-complement(resolved 806R)
```

IUPAC ambiguity codes are resolved reproducibly to A/C/G/T for each feature. R1 begins with a concrete 515F variant and R2 begins with a concrete 806R variant. These files contain the locus-specific primers, but not complete Illumina adapters, indices, EMP pads, or linkers.

## Installation

```bash
mamba env create -f environment.yml
mamba activate mock16s-chem-v1
```

## Full run

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
    n_participants: 6
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
    n_participants: 6
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
samples. Excluding contralateral cancer samples leaves balanced 6-vs-6
patient-level comparisons for both sample types. Set
artifacts.group_differential_abundance.group_columns to impose known effects
for any generated metadata categories. The legacy independent groups layout
remains supported.

## Artifact configuration

All artifacts are explicit and independently switchable in `config/example.yaml`:

```yaml
artifacts:
  group_differential_abundance:
    enabled: true
    group_columns: [Type_Group, Case]
    asvs_per_group: 8
    log_fold_change_sd: 0.25
    min_abs_log_fold_change: 1.75
    candidate_min_prevalence: 0.20
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
    ├── ground_truth_pcr_bias.tsv
    ├── ground_truth_chimeras.tsv
    ├── ground_truth_chemistry_batch.tsv
    ├── ground_truth_mitochondria.tsv
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

## Ground-truth interpretation

- `asv_counts_biological.tsv`: community counts after optional experimental-group effects, before technical artifacts.
- `asv_counts_post_pcr.tsv`: counts after ASV-specific PCR efficiency.
- `asv_counts_final.tsv`: exact per-feature read-pair requests supplied to InSilicoSeq, including contaminants and chimeras.
- `ground_truth_feature_registry.tsv`: stable feature UUID, ASV ID, feature type, exact V4 sequence, sequence hash, representative SILVA record, taxonomy, and original SparseDOSSA2 feature mapping.
- `ground_truth_asv_chem.tsv`: exact nonzero ASV–compound coefficients.
- `ground_truth_group_effects.tsv`: imposed group-specific log fold changes.
- `ground_truth_microbiome_batch_effects.tsv`: ASVs receiving known batch-specific effects, including batch, signed log fold-change, and direction.
- `ground_truth_pcr_bias.tsv`: per-ASV PCR efficiencies.
- `ground_truth_chimeras.tsv`: chimera parents and breakpoint.
- `ground_truth_chemistry_batch.tsv`: compound-specific batch effects.
- `ground_truth_mitochondria.tsv`: mitochondrial source record, exact sequence checksum, configured abundance, observed prevalence, and injected read count. By default DECOI uses a vendored, checksum-pinned segment of the human RefSeq mitochondrial genome, avoiding a runtime network dependency; `reference_fixtures.mitochondria.source_fasta` can select another local fixture and `source_url` remains available explicitly. DECOI synthetically adds the configured amplicon primers to this genuine mitochondrial template so downstream mitochondrial filtering can be tested; it does not assert natural amplification by those primers.
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
- SILVA source SHA-256 and reference preparation metadata.

Feature UUIDs use deterministic UUIDv5 values and remain stable for a given ASV identifier.

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

## Tests

```bash
pytest -q
```

The unit suite tests stable identifiers, primer resolution, study expansion and participant metadata cycling, depth-preserving group and microbial batch effects, contaminants, mitochondria, chimeras, and chemistry batch effects without requiring SILVA or a full R/InSilicoSeq installation.

## Current scope

v1.0 supports SILVA-derived V4 amplicons and InSilicoSeq MiSeq paired-end output. It does not yet model full genomes, KEGG pathways, mechanistic metabolism, index hopping, complete adapter/index constructs, PacBio, or Nanopore reads.

## Citing DECOI and its dependencies

If you use DECOI in a publication, please cite DECOI itself once a project citation is available and cite the software and reference-data publications below. The first group applies to every standard DECOI run. Cite DADA2 only when the optional validation stage is enabled with `--run_dada2 true`.

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

- **SILVA:** Quast, C. *et al.* (2013). The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. *Nucleic Acids Research*, 41(D1), D590–D596. [https://doi.org/10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)
- **Exact default training set:** SILVA 138.2 NR99 taxonomic training data formatted for DADA2 (`silva_nr99_v138.2_toSpecies_trainset.fa.gz`). [https://doi.org/10.5281/zenodo.14169026](https://doi.org/10.5281/zenodo.14169026)

### Dependency-audit note

This list was derived from the workflow definitions and source imports, not solely from the Conda environment. Python, R, PyYAML, `optparse`, `curl`, and GNU `md5sum` are directly used but do not have a single canonical peer-reviewed software DOI that DECOI can recommend. SciPy, SeqKit, pigz, wget, `jsonlite`, `remotes`, and the listed SparseDOSSA2 support packages are present in `environment.yml` for installation or transitive support but are not directly invoked by the current standard workflow, so they are not included as required citations above.
