# Configuration and injected artifacts

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

