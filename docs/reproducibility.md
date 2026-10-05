# Reproducibility and scientific scope

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

## Current scope

DECOI supports SILVA-derived V4 amplicons and InSilicoSeq MiSeq paired-end output, including experimental paired WGS using explicitly mapped genome references, as described above. It does not yet model KEGG pathways, mechanistic metabolism, index hopping, complete adapter/index constructs, PacBio, or Nanopore reads.

