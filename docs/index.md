# DECOI: Data Emulator for Community Omics tool Benchmarking

DECOI generates a mock microbial study with known sequence identities, abundances,
artifacts and relationships to chemical measurements. SparseDOSSA2 supplies the
community profiles; DECOI applies recorded biological and technical effects,
then simulates reads with InSilicoSeq. Matched short-read metagenomics derives
from the same underlying community, rather than a separately sampled study.

The outputs support workflow testing and comparison of statistical methods
against implanted ground truth. Chemical associations are constructed correlations,
not a mechanistic metabolic model. Genome-derived WGS support is experimental;
small integration tests and a full-size study benchmark are different levels of evidence.

```{image} assets/workflow-brief.svg
:alt: DECOI inputs, reference preparation, simulated community and matched output assays
:class: decoi-primary-workflow
```

Start with [installation](installation.md) and the [reviewer test](reviewer-test.md).
The detailed [workflow](workflow.md) explains execution, assay linkage and data flow.

```{toctree}
:maxdepth: 2

installation
reviewer-test
inputs
configuration
workflow
outputs
wgs
paired-study
reference-panel
validation
reproducibility
troubleshooting
citations
```
