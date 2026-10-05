# Outputs and ground truth

## Output layout

The layout below is for the SILVA/Nextflow route. For paired genome-first output,
including `wgs/` and `assay_manifest.tsv`, see [the paired output layout](paired-study.md#6-inspect-completion).
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

