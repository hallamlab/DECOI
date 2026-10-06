# Get and install DECOI

DECOI currently installs from GitHub into a mamba environment on Linux. Mamba
provides Python, R, Nextflow and the sequencing tools; pip installs the `decoi`
command and its bundled workflow. No administrator access is needed for this
installation. Start in a fresh terminal without a Python `.venv` active.

The controller is currently on `refactor/user-guide-installation`. The clone
command below selects that branch so you get the interface documented here.

## 1. Install mamba

Skip this step if `mamba --version` already works. Otherwise, install
[Miniforge](https://github.com/conda-forge/miniforge#install), which provides both
mamba and conda:

```bash
cd ~
curl -fLO "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash "Miniforge3-$(uname)-$(uname -m).sh"
```

Follow the installer prompts and accept shell initialization. Close and reopen
your terminal afterward. These commands assume Bash on Linux and an available
`curl` command; the Miniforge page also provides browser downloads.

## 2. Get DECOI

Install Git if it is not already available, then clone the source:

```bash
mamba install -n base -c conda-forge git
mkdir -p ~/repos
cd ~/repos
git clone --branch refactor/user-guide-installation https://github.com/hallamlab/DECOI.git
cd DECOI
```

## 3. Create the environment and install

```bash
mamba env create -n decoi -f environment.yml
conda activate decoi
Rscript scripts/install_sparsedossa2.R
python -m pip install .
```

Run these commands from the DECOI checkout. Use a new environment named `decoi`;
if that name already exists, activate it for an existing installation or choose
a different name for an independent fresh install.

The R command installs SparseDOSSA2 from the pinned upstream revision. It is
required: environment creation alone does not install it. The final pip command
installs the controller and runtime resources into the active environment.
After installation, `decoi` works from any directory; the checkout is not required
as your working directory.

The pip package includes the workflow, helper scripts, example configurations,
small fixture and frozen genome-panel metadata. It does not include the full
genome panel. Conda-channel and Quay distribution are not published by this setup
yet; use the source installation above.

## 4. Check and run the tiny test

```bash
decoi --version
decoi check
mkdir -p ~/data/decoi-testing
cd ~/data/decoi-testing
decoi test -o first-test
```

`check` should report the required modules and tools as available. `test` creates
eight tiny synthetic genomes locally and runs reference preparation and paired
amplicon/WGS simulation. It requests one CPU and 4 GB per task. Both stages should
finish successfully, producing three libraries per assay, 12 compressed FASTQs,
chemistry, truth tables and an HTML report. No external genome download is needed.

Use a new output directory for an independent test. To check checkpoint reuse:

```bash
decoi test -o first-test --resume
```

Both completed stages should show as cached. Preserve `first-test/.decoi/`,
which contains the work directory and Nextflow cache.

## 5. View the report

```bash
decoi report -o first-test --serve --no-browser
```

Open `http://localhost:8765/report.html`. If DECOI is running on a remote server,
run this command on your own computer in a separate terminal:

```bash
ssh -N -L 8765:localhost:8765 YOUR_USER@YOUR_SERVER
```

Then open the same localhost URL in your computer's browser. Ctrl-C stops the
report server or SSH tunnel. The report server listens on localhost only.

## Later sessions and updates

Activate the environment whenever you open a new terminal:

```bash
conda activate decoi
```

To update from your checkout after changes have been pushed:

```bash
cd ~/repos/DECOI
git pull --ff-only
python -m pip install .
decoi check
```

Run the small test in a fresh output directory to verify an update. If dependency
requirements change, follow the release instructions for updating the supporting
environment too. The environment YAML is not an exact package lock.

Continue with the [test output guide](test.md), [CLI guide](cli.md), or
[full paired airway study](paired-study.md). See [reproducibility](reproducibility.md)
before comparing simulations between machines.
