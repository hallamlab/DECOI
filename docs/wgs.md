# Matched short-read metagenomics

## Matched short-read WGS (experimental)

Enable WGS to generate both assays for every sample in one mock study, including extraction blanks. Sample metadata, chemistry, group effects, microbial batch effects, biological ASV abundances, and contaminant identities are shared. Existing amplicon outputs stay at the dataset root; WGS outputs are under `wgs/`, with a combined `assay_manifest.tsv` linking both assays by `sample_id`.

### Genome-first reference preparation

A reusable airway-prioritized RefSeq panel is described in
[`references/airway_v1/README.md`](reference-panel.md). Its frozen
accession/checksum lock makes the genome pool reproducible across runs and
machines. Three reference genomes from reported reagent-contaminant genera are
reserved with `reference_role=contaminant`; biological community selection cannot
draw from that reserved pool. When a custom reference has no role annotations,
the original random contaminant-selection behavior remains available.
The [selection evidence and DOI table](reference-panel.md)
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

