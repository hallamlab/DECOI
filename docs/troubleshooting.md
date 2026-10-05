# Troubleshooting

## SparseDOSSA2 is missing

Activate the `decoi` environment and run `Rscript scripts/install_sparsedossa2.R`
from the repository root. This package must be installed after environment creation.

## Nextflow cannot find a tool

Activate the same environment in the shell that launches Nextflow. The workflow
uses tools on PATH, not a developer's absolute Conda path. Slurm compute nodes
must see the environment, reference inputs and Nextflow work directory.

## Resource reservations exceed the machine

`--threads` controls simulation workers and the Nextflow simulation CPU request.
Memory reservations are configured in `nextflow.config`. Use an additional
Nextflow configuration file to request resources appropriate to your hardware.
Optional DADA2 validation has its own `--dada2_threads` setting.

## Reference or read-count validation fails

Preserve the task log and input reference. Genome/ASV links and marker copy numbers
must agree with exact sequences; a claimed match is not accepted without sequence
verification. Final FASTQ record counts must agree with the requested counts.
Do not replace a failed frozen-panel checksum with a newer accession.

## Report over SSH

On the remote machine serve the dataset directory with `python -m http.server
8765 --directory PATH_TO_DATASET`. On your computer run `ssh -N -L
8765:localhost:8765 USER@HOST`, then open `http://localhost:8765/report.html`.

Report problems through [GitHub issues](https://github.com/hallamlab/DECOI/issues).
Include the command, tool versions and relevant task logs; omit private study data.
