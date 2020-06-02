# ImmPrint

Immprint is a command-line tool built to identify the origine of T-cell samples. Given two samples of TCR sequences it's able to decide whethever they come from the same individual or not.

## Installation 

Clone or download the repository, then run:
```sh
cd Immprint
pip3 install .
```

## Usage

```sh
immprint sampleA.csv sampleB.csv
```

The two csv files should contain a column `cdr3_nt`, which gives the CDR3 sequence of the receptor in nucleotides. ImmPrint can also use the clonal frequency information if a `count` column is provided.

Additional arguments:

- `-S --no-I`: Don't use the recombination probability (faster but less precise).
- `--no-counts`: Don't use the counts even if provided.
- `-f / --full`: Use both chains of the receptor, in that case the elements of the `cdr3_nt` column should all have this specific format: `TRA:TGCA...ACCTTT;TRB:TGTGC...TCTTC`.

## Examples

The folder "examples" contain a few test datasets. From the base folder:

```sh
immprint examples/P1_a.csv examples/P1_b.csv 
```
Will return:



