# Changelog

## Unreleased — paired short-read metagenomics

- Added experimental amplicon/WGS simulation within one study, sharing biological
  composition, sample IDs, chemistry, contaminants, and extraction controls.
- Added copy-number-adjusted cell equivalents, genome-length-weighted DNA
  sampling, an independent WGS random stream, and per-genome/contig truth tables.
- Added WGS FASTQ validation and a combined sample/assay manifest.
- Added genome-first preparation from local assemblies and supplied taxonomy,
  exact two-strand V4 extraction, copy verification, provenance, and exclusion
  audits. Initial support requires one V4 allele per genome and unique markers
  across genomes; contigs must contain at least 1 kb of unambiguous A/C/G/T.
- Added optional biological/contaminant reference roles and role-aware sampling.
- Added the frozen airway v1 panel with 547 biological references and three
  reagent-contaminant representatives, plus download/rehydration/checksum and
  prepared-reference audit tools, accession locks, and evidence/DOI documentation.
- Added `add-wgs` for preserving an existing study when exact linked genomes cover
  every original biological and contaminant ASV; this does not infer genome
  identity from taxonomy or replace missing mappings with synthetic genomes.
- Added genome-first Nextflow routing, a paired airway configuration, synthetic
  smoke fixtures, and regression/integration tests (29 passed with real tools
  in the development environment).
- Documented new-machine setup, pinned SparseDOSSA2 installation, full-test
  commands, storage requirements, output interpretation, and Nextflow environment
  overrides in `SETUP_WGS.md`. Private application-note drafts remain ignored.
- Status: all 550 panel genomes passed preparation and independent validation;
  the full 64-library real-genome sequencing run is pending. WGS-only runs,
  full host genomes, functional truth, and automatic WGS input cleanup are not
  implemented. The abundance template remains SparseDOSSA2 `Stool`.

## 1.0.0

- Added Nextflow DSL2 orchestration and automatic checksum-verified SILVA download.
- Added SILVA V4 extraction, taxonomy retention, and source mappings.
- Added SparseDOSSA2 count generation and correlated chemistry.
- Added study-design YAML expansion and schema.
- Added group abundance effects, PCR bias, contaminants, chimeras, and chemical batch effects.
- Added deterministic feature UUIDs and complete ground-truth registry.
- Added exact InSilicoSeq per-feature read-count inputs and MiSeq paired FASTQ validation.
- Added machine-readable manifest, HTML report, and optional DADA2 validation.
