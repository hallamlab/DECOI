# Frozen airway reference panel

## DECOI airway-prioritized reference panel v1

This panel supplies one fixed genome pool for paired V4 amplicon and short-read
metagenomics simulations. The lung cohort is defined separately in
`study/example_study.yaml`. Reference selection does not change the inherited
SparseDOSSA2 `Stool` abundance template or establish a calibrated lung ecosystem.

See [the citation and evidence table](citations.md) for verified publication DOIs,
their specific relevance, and the limits of genus-level evidence.
For installation, rehydration, and the full paired-study test on another system,
follow [SETUP_WGS.md](paired-study.md). Genome data are ignored by Git and must
be downloaded from the lock or transferred separately.

## Frozen v1 panel

The completed panel contains **550 bacterial genome/V4 pairs**: 547 biological
references from prioritized genera and three reserved contaminant references.
No background-tier genomes were needed. All 550 were retained during independent
DECOI reference preparation. The configured study can select 500 biological
features and the three contaminants without overlap.

The contaminants represent **Bradyrhizobium**, **Ralstonia**, and **Sphingomonas**;
their exact accessions and the evidence limitations are listed in `CITATIONS.md`.
The lock-file SHA-256 is
`38f773c52f19ec49c7f9201f743e59a4721a57689d7be937f9fb5e664261beaf`.

Selection/preparation used Python 3.11.15, Biopython 1.87, NumPy 2.4.6,
pandas 3.0.3, PyYAML 6.0.3, and curl 8.5.0 (Ubuntu package
8.5.0-2ubuntu10.15). Reference validation is recorded in `validation.json`.

## Selection policy

`policy.yaml` defines a target of 550 unique genome/marker pairs, including three
reserved reagent-contaminant representatives. Candidates come from a frozen
NCBI RefSeq assembly-summary snapshot and NCBI taxonomy snapshot. Only current,
complete bacterial/archaeal assemblies of 0.3–12 Mb are considered. Reference or
representative assemblies are preferred; at most two assemblies per species are
screened, with at most one retained. Candidates are ordered deterministically
using seed 42 and balanced across genera within each priority tier.

Airway/oral-associated genera and respiratory opportunists are prioritized before
background diversity. The priority list is an explicit benchmark design choice,
not evidence that every included species occurs in healthy lungs. It is informed
by primary studies of [supraglottic-associated lung communities](https://pmc.ncbi.nlm.nih.gov/articles/PMC3971609/)
and [bronchial communities in lung cancer cohorts](https://pmc.ncbi.nlm.nih.gov/articles/PMC8900294/).
Background members, if needed to reach the target, remain labelled `background`.

The contaminant pool draws three representatives from distinct genera among
Ralstonia, Bradyrhizobium, Sphingomonas, Methylobacterium, and Microbacterium.
Genus-level reagent/laboratory contamination evidence comes from
[Salter et al. (2014), DOI 10.1186/s12915-014-0087-z](https://doi.org/10.1186/s12915-014-0087-z).
This designation models contamination; it does not assert that the selected
reference strains were isolated from a reagent. `reference_role=contaminant`
excludes these genomes from biological community selection and reserves them for
contaminant injection and extraction controls in both assays.

All genomes must satisfy the current DECOI restrictions: unambiguous A/C/G/T
contigs at least 1 kb long, an exact primer-bounded V4 insert of the configured
length, one distinct V4 allele per assembly, and no V4 allele shared with another
retained assembly. Identical copies within an assembly are supported. Whole
assemblies are retained or excluded; bases and contigs are not edited to qualify.
These restrictions impose ascertainment bias and are appropriate to an initial
software benchmark, not a comprehensive census of airway organisms.

## Frozen files and local data

The completed panel is represented by `panel.lock.json` and `panel.tsv` beside
this file. The lock stores exact versioned accessions, download URLs, NCBI MD5
checksums, SHA-256 checksums of the compressed FASTAs, taxonomy, biological versus
contaminant roles, priority labels, and hashes of the selection inputs. Pin the
lock-file SHA-256 when publishing a run. Do not regenerate v1 from a later live
catalog: create a new panel version for any changed selection.

`build_config.yaml` archives the exact configuration supplied during selection;
its hash is recorded in the lock. Its original directory placeholders are not
the runtime paths. Use `config/airway_paired.yaml` for simulation. Reproducing the
selection itself also requires the archived catalog and taxonomy snapshots;
reusing the panel requires only the lock and the pinned genomes.

Large files live under ignored `data/reference_panels/airway_v1/`:

- `metadata/`: frozen catalog, taxonomy archive, and candidate plan;
- `cache/`: downloaded candidate genomes and NCBI MD5 files;
- `screening.jsonl`: screening outcomes, including excluded candidates;
- `raw/`: selected original genome archives and `taxonomy.tsv`;
- `prepared/`: linked genome/V4 reference produced by DECOI.

The genome sequences are from NCBI RefSeq. See the
[NCBI genome-download documentation](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/)
for the upstream assembly catalog and per-assembly files. Rehydration can only
succeed while the pinned upstream files remain available; local archival copies
and the lock provide the strongest long-term reproducibility.

## Commands

Activate the DECOI environment and run from the repository root. Build a *new*
panel directory using the recorded policy:

```bash
python scripts/build_genome_panel.py build --workers 4
```

`build` is for panel development, not installation of frozen v1. It refuses an
already frozen destination. Changing a panel requires a new policy/panel ID and
output directory. The historical configuration is `build_config.yaml`; the
current runtime config has updated directory paths. Use the rehydration command
below to reproduce the existing panel on another machine.

Verify an existing downloaded panel without network access:

```bash
python scripts/build_genome_panel.py verify \
  --lock references/airway_v1/panel.lock.json \
  --output data/reference_panels/airway_v1/raw
```

Rehydrate exactly the frozen accessions on another machine:

```bash
python scripts/build_genome_panel.py rehydrate \
  --lock references/airway_v1/panel.lock.json \
  --output data/reference_panels/airway_v1/raw
```

Prepare V4 and genome links in a fresh directory:

```bash
python mock16s_chem.py --config config/airway_paired.yaml \
  --genome-dir data/reference_panels/airway_v1/raw \
  --reference-dir data/reference_panels/airway_v1/prepared prepare-reference
```

Simulation uses that prepared reference with the existing study design. Creating
a paired dataset from this genome panel is a new realization of the study and
does not imply exact genomic matches to the historical SILVA-selected ASVs.

Audit the downloaded files and their prepared genome/marker links:

```bash
python scripts/audit_genome_panel.py \
  --lock references/airway_v1/panel.lock.json \
  --raw data/reference_panels/airway_v1/raw \
  --prepared data/reference_panels/airway_v1/prepared \
  --report data/reference_panels/airway_v1/validation.json
```

The current WGS renderer retains per-sample InSilicoSeq input FASTAs. A full
study can therefore use substantially more disk space than the reference panel:
budget for copies of represented genomes across samples as well as reads.


## Evidence supporting the airway and contaminant reference panel

## Scope of the evidence

This is an airway-prioritized **simulation benchmark**, not a published list of
550 validated lung isolates. Literature supports the genus-level ecological
rationale. The exact versioned RefSeq assemblies in `panel.lock.json` were chosen
by the reproducible selection and sequence-eligibility rules in `policy.yaml`.
Genus membership does not establish airway residence, pathogenicity, or reagent
origin for every species or strain in that genus.

The three core sources below informed selection. The supplementary sources add
context; they were documented after selection and did not change the frozen
accessions, selection policy, or simulated abundances.

## Core selection evidence

| Source | DOI | Supported choice and evidence level |
| --- | --- | --- |
| Segal LN et al. (2013). *Enrichment of lung microbiome with supraglottic taxa is associated with increased pulmonary inflammation*. Microbiome 1:19. | [10.1186/2049-2618-1-19](https://doi.org/10.1186/2049-2618-1-19) | BAL evidence for oral-associated lung communities, particularly **Prevotella** and **Veillonella**. Supports genus prioritization, not the selected RefSeq strains. |
| Marshall EA et al. (2022). *Distinct bronchial microbiome precedes clinical diagnosis of lung cancer*. Molecular Cancer 21:68. | [10.1186/s12943-022-01544-6](https://doi.org/10.1186/s12943-022-01544-6) | Bronchial-brushing evidence relevant to a lung-cancer study; **Veillonella**, **Streptococcus**, and **Prevotella** feature prominently. Does not establish every prioritized genus or prescribe the synthetic case/control effects. |
| Salter SJ et al. (2014). *Reagent and laboratory contamination can critically impact sequence-based microbiome analyses*. BMC Biology 12:87. | [10.1186/s12915-014-0087-z](https://doi.org/10.1186/s12915-014-0087-z) | Table 1 documents **Ralstonia**, **Bradyrhizobium**, **Sphingomonas**, **Methylobacterium**, and **Microbacterium** in blank controls. Direct genus-level support for the contaminant candidate list; supports modelling contamination in both amplicon and shotgun data. |

## Supplementary ecological context

| Source | DOI | Relevance and limitation |
| --- | --- | --- |
| Bassis CM et al. (2015). Study of upper-respiratory microbiotas as sources of lung and gastric microbiotas in healthy people. mBio 6:e00037-15. | [10.1128/mBio.00037-15](https://doi.org/10.1128/mBio.00037-15) | Supports the oral-to-lung immigration rationale; oral and lung communities overlap but are not identical. |
| Dickson RP et al. (2015). *Spatial Variation in the Healthy Human Lung Microbiome and the Adapted Island Model of Lung Biogeography*. Annals of the American Thoracic Society 12:821–830. | [10.1513/AnnalsATS.201501-029OC](https://doi.org/10.1513/AnnalsATS.201501-029OC) | Supports oral-associated taxa and lung biogeography, including **Prevotella**, **Veillonella**, and **Streptococcus**. Not a genome-panel prescription. |
| Dewhirst FE et al. (2010). *The Human Oral Microbiome*. Journal of Bacteriology 192:5002–5017. | [10.1128/JB.00542-10](https://doi.org/10.1128/JB.00542-10) | Primary taxonomic resource supporting the oral-community component. Oral occurrence alone is not evidence of lower-airway residence. |
| Fodor AA et al. (2012). Longitudinal study of adult cystic-fibrosis airway microbiota and antibiotic treatment. PLOS ONE 7:e45001. | [10.1371/journal.pone.0045001](https://doi.org/10.1371/journal.pone.0045001) | Observed **Pseudomonas**, **Burkholderia**, **Prevotella**, **Streptococcus**, **Rothia**, **Veillonella**, **Actinomyces**, and **Granulicatella**. Supports inclusion of disease-associated airway diversity, not a healthy-lung or lung-cancer abundance model. |

## Exact contaminant representatives

| Frozen assembly | Reference organism | Evidence used |
| --- | --- | --- |
| `GCF_053592855.1` | Bradyrhizobium sp. T-1 | Salter et al., genus-level blank-control evidence. |
| `GCF_024925465.1` | Ralstonia pseudosolanacearum | Salter et al., genus-level blank-control evidence. |
| `GCF_002151445.1` | Sphingomonas sp. KC8 | Salter et al., genus-level blank-control evidence. |

These are computational representatives of contaminant genera. The cited paper
does **not** identify these exact strains/assemblies as reagent contaminants.
No strain-level reagent-isolation claim is made. In this benchmark their
`reference_role=contaminant` designation reserves them for injection and controls
and excludes them from biological community sampling.

## Limits on interpretation

- The expanded priority-genus list is a benchmark-design choice; the sources
  above are not an individually curated evidence map for every genus or genome.
- Several genera can occur both biologically and as technical contamination.
  A simulated role is ground truth for this experiment, not a universal label
  for interpreting real samples.
- The papers do not specify the panel size, deterministic selection seed,
  per-genus caps, or one-to-one V4/genome restriction. Those are software-design
  decisions documented separately in the policy.
- Abundance, prevalence, contamination intensity, and disease associations are
  simulated settings, not estimates fitted to these publications. The inherited
  SparseDOSSA2 template remains `Stool`.
- Cite the ecological papers for the rationale, and cite/report the frozen
  accessions and lock SHA-256 for the exact sequence inputs. These serve different
  provenance purposes.
