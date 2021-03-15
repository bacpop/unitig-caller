# unitig-caller
[![Dev build Status](https://dev.azure.com/jlees/unitig-caller/_apis/build/status/johnlees.unitig-caller?branchName=master)](https://dev.azure.com/jlees/unitig-caller/_build/latest?definitionId=1&branchName=master)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/unitig-caller/badges/version.svg)](https://anaconda.org/bioconda/unitig-caller)

Determines presence/absence of sequence elements in bacterial sequence
data using Bifrost Build and Query functions. Uses assemblies and/or reads as inputs.

The implementation of unitig-caller is a wrapper around [Bifrost](https://github.com/pmelsted/bifrost)
which formats files for use with pyseer, as well as an implementation which calls sequences
using an FM-index.

Build mode creates a compact de Bruijn graph using Bifrost. Query mode converts the .gfa
file produced by Build mode to a .fasta, using an associated colours file to query
the presence of unitigs in the source genomes used to build the original de Bruijn graph.

Simple mode finds presence of unitigs in a new population using an FM-index.

## Install

Use `unitig-caller` if installed through pip/conda, or
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
pip install unitig-caller
```

Requires [bifrost](https://github.com/pmelsted/bifrost) version 1.0.3 installed, and accessible
via PATH (see steps for installation at Bifrost github page).

### From source
Requires `cmake`, `pthreads`, `pybind11` and a C++17 compiler (e.g. gcc >=7.3), in addition
to the pip requirements.
```
git clone https://github.com/johnlees/unitig-caller --recursive
python setup.py install
```

## Usage

There are three ways to use this package:
1. Build a population graph to extract unitigs for GWAS with pyseer like [unitig-counter](https://github.com/johnlees/unitig-counter) (`--call`).
2. Find existing unitigs in a new population using a graph (`--query`).
3. Find existing unitigs in a new population using an index (`--simple`).

For 1), run `--call` mode.

Both 2) and 3) give the same results with different index tools, both finding unitigs so pyseer models can be applied to a new population.

For 2) Run `--query` mode, specifying *new population* input fastas file names in a text file (one file per line), with `--unitigs` from *the original population*.

For 3), run `--simple` mode giving the new genomes as `--refs` and the `--unitigs` from *the original population*.

These modes are detailed below

### Running Call mode
This uses Bifrost Build to generate a compact coloured de Bruijn graph, and return colours of unitigs within.

#### If no pre-built Bifrost graph exists
```
unitig-caller --call --refs refs.txt --reads reads.txt --out out_prefix
```

`--refs` is a required .txt file listing paths of input ASSEMBLIES
(.fasta or .fastq), each on a new line. No header row.

`--reads` is an optional .txt file listing paths to READ files (.fasta or .fastq), each on new line. No header row

`--kmer` can be specified for the kmer size used to built the graph. By default this is 31 bp.

#### If pre-built Bifrost graph exists

```
unitig-caller --call --graph graph.gfa --colours graph.bfg_colors --out out_prefix
```

`--graph` is a pre-built bifrost graph .gfa, and `--colours` is it's associated colours file.

#### For both call modes

`--out` is the prefix for output files.

Call mode automatically generates a .pyseer file containing unitigs found within the graph and their graph. Rtab or pyseer
formats can be specified with `--rtab` and `--pyseer` respectively.

### Running Query mode
Queries existing unitigs in a Bifrost graph.

#### If no pre-built Bifrost graph exists
```
unitig-caller --query --refs refs.txt --reads reads.txt --unitigs query_unitigs.fasta --out out_prefix
```

`refs` and `reads` are the same arguments as in `call`.

`--kmer` can be specified for the kmer size used to built the graph. By default this is 31 bp.

#### If pre-built Bifrost graph exists

```
unitig-caller --query --graph graph.gfa --colours graph.bfg_colors --unitigs query_unitigs.fasta --out out_prefix
```

#### For both query modes

`--unitigs` is .fasta file or text file with unitig sequences (one sequence per line, with header line).

`--out` is the prefix for output files.

Query mode automatically generates a .pyseer file containing unitigs found within the graph and their graph. Rtab or pyseer
formats can be specified with `--rtab` and `--pyseer` respectively.

### Running simple mode
This uses suffix arrays (FM-index) provided by [SeqAn3](https://www.seqan.de/) to perform
string matches:
```
unitig-caller --simple --refs strain_list.txt --unitigs queries.txt --output calls
```

`--refs` is a required file listing input assemblies, the same as `refs` in `call`.

`--unitigs` is a required list of the unitig sequences to call. The unitigs need
to be in the first column (tab separated). A header row is assumed, so
output from [pyseer](https://github.com/mgalardini/pyseer) etc can be directly used.

`calls_pyseer.txt` will contain unitig calls in seer/pyseer k-mer format.

By default FM-indexes are saved in the same location as the assembly files so that they can
be quickly loaded by subsequent runs. To turn this off use `--no-save-idx`.

### Option reference
```
usage: unitig-caller [-h] (--call | --query | --simple) [--refs REFS]
                     [--reads READS] [--graph GRAPH] [--colours COLOURS]
                     [--unitigs UNITIGS] [--pyseer] [--rtab] [--out OUT]
                     [--kmer KMER] [--write-graph]
                     [--no-save-idx] [--threads THREADS] [--version]

Call unitigs in a population dataset

optional arguments:
  -h, --help         show this help message and exit

Mode of operation:
  --call             Build a DBG and call colours of unitigs within
  --query            Query unitig colours in reference genomes/DBG
  --simple           Use FM-index to make calls

Unitig-caller input/output:
  --refs REFS        Ref file to used to build DBG or use with --simple
  --reads READS      Read file to used to build DBG
  --graph GRAPH      Existing graph in GFA format
  --colours COLOURS  Existing bifrost colours file in .bfg_colors format
  --unitigs UNITIGS  Text or fasta file of unitigs to query (--query or --simple)
  --pyseer           Output pyseer format
  --rtab             Output rtab format
  --out OUT          Prefix for output [default = 'unitig_caller']

Bifrost options:
  --kmer KMER        K-mer size for graph building/querying [default = 31]
  --write-graph      Output DBG built with unitig-caller

Simple mode options:
  --no-save-idx      Do not save FM-indexes for reuse

Other:
  --threads THREADS  Number of threads to use [default = 1]
  --version          show program's version number and exit
```

## Citation

If you use this, please cite the Bifrost paper:

Holley G., Melsted, P. Bifrost – Highly parallel construction and indexing of colored and compacted de Bruijn graphs.
bioRxiv 695338 (2019). doi: https://doi.org/10.1101/695338
