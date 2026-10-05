# Inputs and study designs

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

