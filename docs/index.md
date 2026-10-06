# DECOI: Data Emulator for Community Omics tool Benchmarking

When you test a microbiome analysis workflow on real samples, you rarely know
exactly what it should recover. DECOI gives you a study where you do: the organisms,
their abundances, the effects you introduced, and their relationships to chemical
measurements are recorded alongside the simulated data. You can use that record
to check whether a workflow runs correctly and how well an analysis recovers the
patterns you put into the study.

You can generate paired-end V4 16S amplicon reads on their own, or add matched
shotgun metagenomic reads (WGS) from linked reference genomes. Both assays start
from the same simulated community and share sample IDs, so you can follow the
same study through different analysis workflows. DECOI also generates chemical
measurements with known microbial associations, plus the sample metadata and
reference tables needed to interpret the results.

The study can include participant groups, repeated samples, differences in
microbial abundance, microbial modules, and batch effects. You can also introduce
PCR bias, contaminants, extraction blanks, chimeras, and mitochondrial marker
sequences to test how your workflow handles imperfect data. These effects are
recorded in ground-truth tables. WGS uses genome-linked community abundances,
with marker-copy and genome-length adjustments; amplicon-only artifacts stay
with the amplicon assay.

Follow the [quickstart](installation.md) and run the [tiny test](reviewer-test.md)
to get started. You can then use the supplied [airway mock microbial community](reference-panel.md) or provide
your own references and study design. Current support covers V4 amplicons
and paired short-read WGS with the InSilicoSeq MiSeq model, run locally or through
Nextflow with a Slurm executor. Optional DADA2 validation is available for
amplicons. The full paired airway demo has passed generation and output-consistency
checks; WGS support remains experimental. Chemical associations are deliberately
constructed for testing and do not represent a model of microbial metabolism.

```{image} assets/workflow-brief.svg
:alt: DECOI inputs, reference preparation, simulated community and matched output assays
:class: decoi-primary-workflow
```

[Quickstart](installation.md) · [Run the test](reviewer-test.md) · [Command-line guide](cli.md)

The detailed [workflow](workflow.md) explains execution, assay linkage and data flow.

```{toctree}
:maxdepth: 2

installation
cli
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
