# unitig-caller
Determines presence/absence of sequence elements in bacterial sequence
data. Currently uses assemblies as inputs (reads also possible by changing
some options).

For small numbers of queries use `--simple`. For larger numbers
(e.g. whole genome) use `--index` followed by `--call`.

Simple mode uses string matching after loading seqs into main memory.
Index/call mode is a wrapper around squeakr and mantis.

## Install

Use `unitig_caller` if installed through pip/conda, or
`python unitig_caller-runner.py` if using a clone of the code.

### With conda (recommended)
Get it from [bioconda](http://bioconda.github.io/):
```
conda install unitig-caller
```

If you haven't set this up, first install
[miniconda](https://docs.conda.io/en/latest/miniconda.html). Then
add the correct channels:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

### With pip
Get it from PyPI:
```
pip install unitig_caller
```

Requires [squeakr](https://github.com/splatlab/squeakr) version 1.0
and [mantis](https://github.com/splatlab/mantis) version 0.2.0 to run
the `--index` and `--call` modes.

### From source
Requires `cmake`, `pthreads` and a C++11 compiler (e.g. gcc >=4.8), in addition
to the pip requirements.
```
git clone https://github.com/johnlees/unitig-caller
python setup.py install
```

## Usage

### Running simple mode
This uses code from [seer](https://github.com/johnlees/seer) to perform
simple string matches:
```
unitig-caller --mode simple --strains strain_list.txt --unitigs queries.txt --output calls.txt
```

`--strains` is a file listing input assemblies, name followed by location
of fasta file (tab separated), each on a new line. No header row.

`--unitigs` is a list of the unitig sequences to call. The unitigs need
to be in the first column (tab separated). A header row is assumed, so
output from [pyseer](https://github.com/mgalardini/pyseer) etc can be directly used.

`calls_unitigs.txt` will contain unitig calls in seer/pyseer k-mer format.

### Running with mantis
First create an index of all the sequences:
```
unitig_caller --index --strains strain_list.txt --output mantis_idx --cpus 8
```
Note this will require about 20Mb temporary disk space per strain, and the final step
may be fairly RAM intensive. If this fails you can re-run without having to
generate all the intermediate files again. Use of multiple CPUs is recommended, if available.

You can then use this index to rapidly call many unitig queries:
```
unitig_caller --call --unitigs queries.txt --output calls --mantis-index mantis_idx
```
Which will write to `calls_unitigs.txt` as in simple mode.

### Option reference
```
usage: unitig-caller [-h] (--index | --call | --simple) [--strains STRAINS]
                     [--unitigs UNITIGS] [--output OUTPUT]
                     [--mantis-index MANTIS_INDEX] [--approximate]
                     [--kmer-size KMER_SIZE] [--count-cutoff COUNT_CUTOFF]
                     [--log-slots LOG_SLOTS] [--cpus CPUS] [--overwrite]
                     [--squeakr SQUEAKR] [--mantis MANTIS] [--version]

Call unitigs in a population dataset

optional arguments:
  -h, --help            show this help message and exit

Mode of operation:
  --index               Index sequences, before calling
  --call                Make calls from an indexed dataset
  --simple              Use string matching to make calls

Input/output:
  --strains STRAINS     List of strains to index
  --unitigs UNITIGS     List of unitigs to call
  --output OUTPUT       Prefix for output [default = 'mantis_index']
  --mantis-index MANTIS_INDEX
                        Directory containing mantis index (produced by index
                        mode) [default = 'mantis_index']

mantis/squeakr options:
  --approximate         Use approximate count mode [default = exact]
  --kmer-size KMER_SIZE
                        K-mer size for counts [default = 28]
  --count-cutoff COUNT_CUTOFF
                        Minimum k-mer count to be included [default = 1]
  --log-slots LOG_SLOTS
                        Starting log(number of count slots). Automatically
                        increased if necessary [default = 22]

Other:
  --cpus CPUS           Number of CPUs to use [default = 1]
  --overwrite           Overwrite existing output
  --squeakr SQUEAKR     Location of squeakr executable [default = squeakr]
  --mantis MANTIS       Location of mantis executable [default = mantis]
  --version             show program's version number and exit
```

## Timing
Testing on 603 *S. pneumoniae* assemblies:

| Test          | Mode          | Cores | Time    | Memory | Disk  | +tmp disk |
| ------------- | ------------- |-------|---------|--------|-------|---------- |
| Index         | Index         | 1     | 25.0min | 12.7Gb | 287Mb | 10.7Gb    |
|               | Index         | 4     | 13.3min | 13.0Gb | 287Mb | 10.7Gb    |
| 100 unitigs   | Simple        | 1     | 27.4min | 1.1Gb  | None  | None      |
|               | Simple        | 4     | 9.3min  | 1.1Gb  | None  | None      |
|               | Call          | 1     | 0.6s    | 34.1Mb | Index | Index     |
| 100k unitigs  | Call          | 1     | 2.9min  | 362Mb  | Index | Index     |