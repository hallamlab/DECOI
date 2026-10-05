# Software and reference citations

## Citing DECOI and its dependencies

If you use DECOI in a publication, please cite DECOI itself once a project citation is available and cite the software and reference-data publications relevant to your run. Nextflow applies to orchestrated runs; Cutadapt applies to SILVA preparation or DADA2 primer removal. Genome-first Python runs do not use SILVA or Cutadapt for reference preparation. Cite DADA2 only when the optional validation stage is enabled with `--run_dada2 true`.

### Core workflow software

- **Nextflow:** Di Tommaso, P. *et al.* (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35, 316–319. [https://doi.org/10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)
- **SparseDOSSA2:** Ma, S. *et al.* (2021). A statistical model for describing and simulating microbial community profiles. *PLOS Computational Biology*, 17(9), e1008913. [https://doi.org/10.1371/journal.pcbi.1008913](https://doi.org/10.1371/journal.pcbi.1008913)
- **Cutadapt:** Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*, 17(1), 10–12. [https://doi.org/10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200)
- **InSilicoSeq:** Gourlé, H., Karlsson-Lindsjö, O., Hayer, J., and Bongcam-Rudloff, E. (2019). Simulating Illumina metagenomic data with InSilicoSeq. *Bioinformatics*, 35(3), 521–522. [https://doi.org/10.1093/bioinformatics/bty630](https://doi.org/10.1093/bioinformatics/bty630)
- **NumPy:** Harris, C. R. *et al.* (2020). Array programming with NumPy. *Nature*, 585, 357–362. [https://doi.org/10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)
- **pandas:** McKinney, W. (2010). Data structures for statistical computing in Python. *Proceedings of the 9th Python in Science Conference*, 56–61. [https://doi.org/10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a)
- **Biopython:** Cock, P. J. A. *et al.* (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422–1423. [https://doi.org/10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)

### Optional validation software

- **DADA2** (only when `--run_dada2 true`): Callahan, B. J. *et al.* (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13, 581–583. [https://doi.org/10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)

### Reference data

For the genome-first airway panel, report the exact RefSeq assembly accessions
and lock SHA-256 and cite the applicable [airway/reagent evidence](reference-panel.md).
The SILVA citations below apply only when that reference is used; they do not
describe the genome-first panel.

- **SILVA:** Quast, C. *et al.* (2013). The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. *Nucleic Acids Research*, 41(D1), D590–D596. [https://doi.org/10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)
- **Exact default training set:** SILVA 138.2 NR99 taxonomic training data formatted for DADA2 (`silva_nr99_v138.2_toSpecies_trainset.fa.gz`). [https://doi.org/10.5281/zenodo.14169026](https://doi.org/10.5281/zenodo.14169026)

### Dependency-audit note

This list was derived from the workflow definitions and source imports, not solely from the Conda environment. Python, R, PyYAML, `optparse`, `curl`, and GNU `md5sum` are directly used but do not have a single canonical peer-reviewed software DOI that DECOI can recommend. SciPy, SeqKit, pigz, wget, `jsonlite`, `remotes`, and the listed SparseDOSSA2 support packages are present in `environment.yml` for installation or transitive support but are not directly invoked by the current standard workflow, so they are not included as required citations above.

