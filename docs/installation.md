# Installation

Use a Linux machine with mamba and Git available. This installation keeps Python,
R and the read simulators in one supporting environment.

```bash
git clone https://github.com/hallamlab/DECOI.git
cd DECOI
mamba env create -f environment.yml
conda activate decoi
Rscript scripts/install_sparsedossa2.R
```

The R setup command installs SparseDOSSA2 from the pinned upstream revision.
It is required: creating the mamba environment alone does not install that package.
Nextflow uses the activated environment; no development-machine path is assumed.
Run commands from the repository root unless a command explicitly gives another path.

Check the installation:

```bash
python --version
nextflow -version
iss --version
Rscript -e 'stopifnot(requireNamespace("SparseDOSSA2", quietly=TRUE))'
```

Continue with the [small reviewer dataset](reviewer-test.md), then the
[paired-study guide](paired-study.md) for the full frozen genome panel.
The environment file specifies versions but is not a complete dependency lock.
See [reproducibility](reproducibility.md) before comparing simulations between machines.
