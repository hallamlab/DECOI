# Optional DADA2 validation

## DADA2 validation

The optional validation stage:

1. removes 515F and 806R with Cutadapt;
2. filters paired reads;
3. learns forward and reverse error models;
4. denoises each direction;
5. merges pairs;
6. removes bimeras;
7. writes `dada2_sequence_table.tsv` and `dada2_tracking.tsv`.

The supplied settings are reasonable defaults, not universal optimums. Adjust truncation lengths for the selected InSilicoSeq MiSeq model when needed.
This stage validates amplicon reads only; it is not a WGS classifier or assembly
validation. The Python `simulate` command does not launch this separate stage.


## Completed full paired study

The frozen airway-panel demo completed on a local Linux server using the direct
Python simulation command, seed 42 and 32 simulator workers for each assay.
The run used Python 3.11.16, InSilicoSeq 2.0.1, R 4.3.3 and Cutadapt 5.2.

| Output | Verified result |
| --- | ---: |
| Biological samples | 60 |
| Extraction controls | 4 |
| Biological sequence features | 500 |
| Final amplicon features | 509 |
| WGS reference genomes | 503 |
| Amplicon read pairs | 2,466,881 |
| WGS read pairs | 6,400,000 |
| Compressed FASTQ files | 256 |
| Chemical measurements | 60 samples × 6 compounds |

An independent output audit read every FASTQ in full and checked gzip integrity,
record structure, sequence/quality lengths and requested read counts. Sample IDs
agreed across both assays and metadata. Feature sequence hashes and ground-truth
identifiers were consistent. WGS shared counts matched the pre-PCR biological
state, with the configured contaminant/control allocations. Marker-copy and
genome-length transformations reconciled with the saved abundance tables, and
per-contig read allocations summed exactly to per-genome allocations. The report
and manifest agreed with the outputs.

This demonstrates successful full-dataset generation and internal consistency.
DADA2 validation was disabled; this run does not establish downstream ASV recovery,
biological realism, or full-study Slurm performance. The Nextflow route was tested
separately with the small paired test dataset, including checkpoint reuse.
