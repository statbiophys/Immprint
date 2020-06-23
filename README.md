# ImmPrint

Immprint is a command-line tool built to identify the origine of T-cell samples. Given two samples of TCR sequences it's able to decide whethever they come from the same individual or not.

## Installation 

Immprint requires a recent version of python3 and pip3.
Clone or download the repository, then run:
```sh
cd Immprint
pip3 install ./olga3 .
```

## Usage

```sh
immprint sampleA.csv sampleB.csv
```

The two csv files should contain a column `cdr3_nt`, which gives the CDR3 sequence of the receptor in nucleotides. ImmPrint can also use the clonal frequency information if a `count` column is provided.

Additional arguments:

- `-h --help`: Print the help
- `-S --no-I`: Don't use the recombination probabilities (faster but less precise).
- `-n, --no-counts`: Don't use the counts even if provided.
- `-f / --full`: Use both chains of the receptor, in that case the elements of the `cdr3_nt` column should all have this specific format: `TRA:TGCA...ACCTTT;TRB:TGTGC...TCTTC`.
- `-g, --gamma`: modify the value of gamma, that influences the precision of Immprint (if using the recombination probabilities).

## Examples

The folder "examples" contain a few test datasets. From the base folder:
```sh
	immprint examples/P1_a.csv examples/P1_b.csv 
```
Or with full receptors:
```sh
	immprint --full examples/P1_full_a.csv examples/P1_full_b.csv
```

## Caveats



