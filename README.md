# unitig-caller
A wrapper around mantis to determine unitig presence/absence

For small numbers of queries use `--simple`. For larger numbers
(e.g. whole genome) use `--index` followed by `--call`.

## Install
Requires [squeakr](https://github.com/splatlab/squeakr) version 1.0
and [mantis](https://github.com/splatlab/mantis) version 0.2.0 to run
the `--index` and `--call` modes.

Use `unitig_caller` if installed through pip/conda, or
`python unitig_caller-runner.py` if using a clone of the code.

## Running simple mode
This uses code from [seer](https://github.com/johnlees/seer) to perform
simple string matches:
```
unitig-caller --mode simple --strains strain_list.txt --unitigs queries.txt --output unitig_calls.txt
```

`--strains` is a file listing input assemblies, name followed by location
of fasta file (tab separated), each on a new line. No header row.

`--unitigs` is a list of the unitig sequences to call. The unitigs need
to be in the first column (tab separated). A header row is assumed, so
output from [pyseer](https://github.com/mgalardini/pyseer) etc can be directly used.

`unitig_calls.txt` will contain unitig calls in seer/pyseer k-mer format.

## Running with mantis
First create an index of all the sequences:
```
unitig_caller --index --strains strain_list.txt --output mantis_idx --cpus 8
```
Note this will require about 20Mb temporary disk space per strain, and the final step
may be fairly RAM intensive. If this fails you can re-run without having to
generate all the intermediate files again. Use of multiple CPUs is recommended, if available.

You can then use this index to rapidly call many unitig queries:
```
unitig_caller --call --unitigs queries.txt --output unitig_calls.txt --mantis-index mantis_idx
```
Which will write to `unitig_calls.txt` as in simple mode.

## Usage
```
usage: unitig-caller [-h] (--index | --call | --simple) [--strains STRAINS]
                     [--unitigs UNITIGS] [--output OUTPUT]
                     [--mantis-index MANTIS_INDEX] [--cpus CPUS] [--overwrite]
                     [--squeakr SQUEAKR] [--mantis MANTIS] [--version]

Call unitigs in a population

optional arguments:
  -h, --help            show this help message and exit

Mode of operation:
  --index               Index sequences, before calling.
  --call                Make calls from an indexed dataset.
  --simple              Use string matching to make calls.

Input/output:
  --strains STRAINS     List of strains to index
  --unitigs UNITIGS     List of unitigs to call
  --output OUTPUT       Prefix for output
  --mantis-index MANTIS_INDEX
                        Directory containing mantis index (produced by index
                        mode)

Other:
  --cpus CPUS           Number of CPUs to use. [default = 1]
  --overwrite           Overwrite existing output
  --squeakr SQUEAKR     Location of squeakr executable [default = squeakr]
  --mantis MANTIS       Location of mantis executable [default = mantis]
  --version             show program's version number and exit
```