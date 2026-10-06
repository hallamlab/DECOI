# Airway mock microbial community

DECOI includes a reference genome collection for generating mock communities
inspired by the human oral and respiratory microbiota. It lets you simulate V4
16S amplicons and shotgun metagenomic reads from the same organisms, with known
abundances and deliberately introduced contaminants.

The collection contains **550 bacterial genomes**. Of these, **547 are available
as community members** and **three are reserved for simulated contamination**.
The supplied study selects 500 community members and adds the three contaminant
references. Each genome has a corresponding V4 sequence, so the amplicon and
metagenomic outputs can be traced back to the same reference organisms.

This is a community designed for testing analysis methods. It is not intended to
reproduce the composition of a particular patient's airway or a typical healthy
lung. The study design is defined separately in `study/example_study.yaml`, and
abundances are generated using the SparseDOSSA2 `Stool` template with the effects
specified in the configuration. They have not been fitted to airway observations.

To use these genomes, follow the [full mock-community setup](paired-study.md).
The [tiny installation test](reviewer-test.md) uses synthetic sequences instead
and does not require this download.

## Which organisms are included, and why?

Selection gives priority to genera associated with the oral cavity and airways,
including respiratory opportunists. The rationale draws on studies of
[oral-associated lung communities](https://doi.org/10.1186/2049-2618-1-19) and
[bronchial communities in lung cancer](https://doi.org/10.1186/s12943-022-01544-6).
The studies and their relevance are listed below.

These studies support the choice of microbial groups. They do not establish
that every selected species or strain lives in the human airway. The expanded
genus list is a simulation-design choice, and all 547 community references came
from that list. Additional genera could have been included if needed to reach
the target number, but were not needed for this collection.

Three genomes represent potential reagent or laboratory contamination:
**Bradyrhizobium**, **Ralstonia**, and **Sphingomonas**. They were selected from a
candidate list that also included Methylobacterium and Microbacterium, based on
[Salter et al. (2014)](https://doi.org/10.1186/s12915-014-0087-z). In the simulation,
these three references are kept separate from the biological community and used
to add contamination and generate extraction controls in both assays. This role
does not imply that the particular reference strains were isolated from reagents.

## How were the reference genomes chosen?

The selection settings are recorded in `references/airway_v1/policy.yaml`.
Candidate genomes came from saved copies of the NCBI RefSeq assembly catalogue
and NCBI taxonomy. The selection considered complete bacterial and archaeal
assemblies between 0.3 and 12 Mb; all retained genomes were bacterial. NCBI
reference or representative assemblies were preferred.

At most two assemblies per species were screened, with no more than one retained.
Selection was balanced across genera within each group of candidates, using a
fixed random seed of 42 so the procedure could be reproduced.

Each retained assembly had to meet the sequence requirements used by DECOI:

- Contigs contain only unambiguous A, C, G and T bases and are at least 1 kb long.
- The configured primers identify a V4 region within the accepted length range.
- The assembly contains one distinct V4 sequence, although identical copies are allowed.
- That V4 sequence is not shared with another retained genome.

Assemblies were accepted or excluded as a whole; sequences were not edited to
make them qualify. These requirements make it possible to link amplicon features
to individual genomes, but they also exclude some organisms and reduce the
biological diversity represented by the mock community.

## Download and check the genomes

The genome sequences are downloaded separately from NCBI. GitHub contains the
selection settings, genome inventory and information needed to retrieve the same
sequence files. Activate the DECOI environment and run these commands from the
repository root.

Download the recorded genome versions:

```bash
python scripts/build_genome_panel.py rehydrate \
  --lock references/airway_v1/panel.lock.json \
  --output data/reference_panels/airway_v1/raw
```

Here, `rehydrate` is the script's command for downloading the existing collection.
It checks the downloaded files against the recorded checksums. To check an
existing download again without network access:

```bash
python scripts/build_genome_panel.py verify \
  --lock references/airway_v1/panel.lock.json \
  --output data/reference_panels/airway_v1/raw
```

Prepare the V4 sequences and their genome links in a new directory:

```bash
decoi prepare-reference --config config/airway_paired.yaml \
  -o data/reference_panels/airway_v1/prepared
```

Check the prepared sequences, genome identities, marker copy numbers and assigned
community/contaminant roles:

```bash
python scripts/audit_genome_panel.py \
  --lock references/airway_v1/panel.lock.json \
  --raw data/reference_panels/airway_v1/raw \
  --prepared data/reference_panels/airway_v1/prepared \
  --report data/reference_panels/airway_v1/validation.json
```

All 550 genomes passed this check. Continue with the
[full study instructions](paired-study.md) to generate reads and chemical
measurements. A new simulation from these genomes does not preserve the ASV
identities of an earlier, independently selected SILVA-based study.

WGS simulation retains a genome FASTA for each sample's read-generation step.
Allow disk space for these files as well as the final reads; a full study uses
more space than the downloaded reference collection alone.

## Keep track of the exact sequences used

The files under `references/airway_v1/` document the reference collection:

- `panel.tsv`: the genome inventory.
- `panel.lock.json`: exact assembly versions, download locations, taxonomy,
  checksums, and the assigned community or contaminant role.
- `policy.yaml`: the selection settings.
- `build_config.yaml`: the configuration used when selecting genomes.
- `validation.json`: the recorded reference checks.
- `CITATIONS.md`: the literature supporting the organism choices.

The filenames retain the term `panel` for compatibility with the download tools.
The SHA-256 checksum of `panel.lock.json` is
`38f773c52f19ec49c7f9201f743e59a4721a57689d7be937f9fb5e664261beaf`.
Record this checksum with your results to identify the exact genome inventory.

Selection originally used Python 3.11.15, Biopython 1.87, NumPy 2.4.6,
pandas 3.0.3, PyYAML 6.0.3 and curl 8.5.0. These describe the selection environment;
they are not an additional installation requirement.

Large files are stored under the Git-ignored `data/reference_panels/airway_v1/`
directory. `raw/` contains the downloaded genomes and taxonomy; `prepared/`
contains the linked genome/V4 references. A selection run can also create
`metadata/`, `cache/` and `screening.jsonl`, recording source catalogues, candidate
downloads and the reasons candidates were accepted or excluded.

Use the recorded accessions when repeating a study. Choosing genomes from a newer
NCBI catalogue creates a different reference collection, even with the same
selection settings. Keep local copies if long-term availability matters: the
original files may not remain downloadable indefinitely. See the
[NCBI genome-download documentation](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/)
for details of the source files.

## Choosing a different set of genomes

For most users, downloading the existing collection or providing their own
references is sufficient. To repeat the genome-selection procedure with different
settings, use the developer script's `build` command with a new collection ID and
output directory. It refuses to replace an existing completed collection.
Reproducing the original selection also requires the saved NCBI catalogue and
taxonomy versions; downloading the already selected genomes does not.

## What does the literature support?

The three studies below informed the choice of airway-associated genera and
potential contaminants. The additional readings provide ecological context;
they were documented after selection and did not change the chosen genomes or
their simulated abundances. Genus membership alone does not establish airway
residence, pathogenicity or reagent origin for an individual species or strain.

## Studies used to choose the reference organisms

| Source | DOI | Supported choice and evidence level |
| --- | --- | --- |
| Segal LN et al. (2013). *Enrichment of lung microbiome with supraglottic taxa is associated with increased pulmonary inflammation*. Microbiome 1:19. | [10.1186/2049-2618-1-19](https://doi.org/10.1186/2049-2618-1-19) | BAL evidence for oral-associated lung communities, particularly **Prevotella** and **Veillonella**. Supports genus prioritization, not the selected RefSeq strains. |
| Marshall EA et al. (2022). *Distinct bronchial microbiome precedes clinical diagnosis of lung cancer*. Molecular Cancer 21:68. | [10.1186/s12943-022-01544-6](https://doi.org/10.1186/s12943-022-01544-6) | Bronchial-brushing evidence relevant to a lung-cancer study; **Veillonella**, **Streptococcus**, and **Prevotella** feature prominently. Does not establish every prioritized genus or prescribe the synthetic case/control effects. |
| Salter SJ et al. (2014). *Reagent and laboratory contamination can critically impact sequence-based microbiome analyses*. BMC Biology 12:87. | [10.1186/s12915-014-0087-z](https://doi.org/10.1186/s12915-014-0087-z) | Table 1 documents **Ralstonia**, **Bradyrhizobium**, **Sphingomonas**, **Methylobacterium**, and **Microbacterium** in blank controls. Direct genus-level support for the contaminant candidate list; supports modelling contamination in both amplicon and shotgun data. |

## Further reading on airway microbial communities

| Source | DOI | Relevance and limitation |
| --- | --- | --- |
| Bassis CM et al. (2015). Study of upper-respiratory microbiotas as sources of lung and gastric microbiotas in healthy people. mBio 6:e00037-15. | [10.1128/mBio.00037-15](https://doi.org/10.1128/mBio.00037-15) | Supports the oral-to-lung immigration rationale; oral and lung communities overlap but are not identical. |
| Dickson RP et al. (2015). *Spatial Variation in the Healthy Human Lung Microbiome and the Adapted Island Model of Lung Biogeography*. Annals of the American Thoracic Society 12:821–830. | [10.1513/AnnalsATS.201501-029OC](https://doi.org/10.1513/AnnalsATS.201501-029OC) | Supports oral-associated taxa and lung biogeography, including **Prevotella**, **Veillonella**, and **Streptococcus**. Not a list of genomes to include. |
| Dewhirst FE et al. (2010). *The Human Oral Microbiome*. Journal of Bacteriology 192:5002–5017. | [10.1128/JB.00542-10](https://doi.org/10.1128/JB.00542-10) | Primary taxonomic resource supporting the oral-community component. Oral occurrence alone is not evidence of lower-airway residence. |
| Fodor AA et al. (2012). Longitudinal study of adult cystic-fibrosis airway microbiota and antibiotic treatment. PLOS ONE 7:e45001. | [10.1371/journal.pone.0045001](https://doi.org/10.1371/journal.pone.0045001) | Observed **Pseudomonas**, **Burkholderia**, **Prevotella**, **Streptococcus**, **Rothia**, **Veillonella**, **Actinomyces**, and **Granulicatella**. Supports inclusion of disease-associated airway diversity, not a healthy-lung or lung-cancer abundance model. |

## Genomes used to simulate contamination

| NCBI assembly accession | Reference organism | Evidence used |
| --- | --- | --- |
| `GCF_053592855.1` | Bradyrhizobium sp. T-1 | Salter et al., genus-level blank-control evidence. |
| `GCF_024925465.1` | Ralstonia pseudosolanacearum | Salter et al., genus-level blank-control evidence. |
| `GCF_002151445.1` | Sphingomonas sp. KC8 | Salter et al., genus-level blank-control evidence. |

These are computational representatives of contaminant genera. The cited paper
does **not** identify these exact strains/assemblies as reagent contaminants.
No strain-level reagent-isolation claim is made. In this benchmark their
`reference_role=contaminant` designation reserves them for simulated contamination and extraction controls
and excludes them from biological community sampling.

## Limits on interpretation

- The expanded priority-genus list is a benchmark-design choice; the sources
  above are not an individually curated evidence map for every genus or genome.
- Several genera can occur both biologically and as technical contamination.
  A simulated role is ground truth for this experiment, not a universal label
  for interpreting real samples.
- The papers do not specify the number of reference genomes, random seed,
  maximum genomes per genus, or the requirement for a unique V4 sequence per genome. Those are software-design
  decisions documented separately in the policy.
- Abundance, prevalence, contamination intensity, and disease associations are
  simulated settings, not estimates fitted to these publications. The inherited
  SparseDOSSA2 template remains `Stool`.
- Cite the ecological papers for the rationale, and cite/report the recorded
  accessions and the checksum of the genome inventory for the exact sequence inputs. These serve different
  reproducibility purposes.
