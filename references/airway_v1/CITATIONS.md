# Evidence supporting the airway and contaminant reference panel

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
