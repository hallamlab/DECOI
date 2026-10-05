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

