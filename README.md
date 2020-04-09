# unitig-caller

Determines presence/absence of sequence elements in bacterial sequence
data using Bifrost Build and Query functions. Uses assemblies and/or reads as inputs.

The implementation of unitig-caller is a wrapper around Bifrost (https://github.com/pmelsted/bifrost)

Build mode creates a compact de Bruijn graph using Bifrost. Query mode converts the .gfa
file produced by Build mode to a .fasta, using an associated colours file to query
the presence of unitigs in the source genomes used to build the original de Bruijn graph.

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
Requires `cmake`, `pthreads` and a C++11 compiler (e.g. gcc >=4.8), in addition
to the pip requirements.
```
git clone https://github.com/samhorsfield96/unitig-caller
python setup.py install
```

## Usage

### Running Build mode
This uses Bifrost Build to generate a compact de Bruijn graph. By default this a
coloured compact de Bruijn graph.
```
unitig-caller --build --seq refs.txt --addit_seq reads.txt --output out_prefix
```

`--seq` is a .txt file listing paths of input assemblies or read files (.fasta or .fastq),
each on a new line. Must be specified as either 'refs.txt' for assemblies
or 'reads.txt' for read files. No header row.

`--addit_seq` is an optional .txt file listing paths to additional sequence files of different type
to those specified in --seq (e.g. if 'refs.txt' is given in --seq, then 'reads.txt' will
be given in --addit_seq and vice versa), each on new line. No header row.

`--output` is the prefix for output files.

By default de Bruijn graphs are coloured, with an accompanying .bfg_colors being
generated alongside the .gfa file. To turn this off, use `--no_colour`. Note, Query mode
cannot be run without a .bfg_colors file.

To generate a clean de Bruijn graph (clip tips and delete isolated contigs shorter
than k k-mers in length), specify `--clean`.

### Running Query mode
Before running Query mode, generate a coloured compact de Bruijn graph using build mode.
Then run the Query command as below.
```
unitig-caller --query --input in_prefix --output out_prefix
```

`--input` is the prefix for the .gfa and .bfg_colors files generated from Build mode.

`--output` is the prefix for output files.

To query unitigs that are not present in the graph, specify `--fasta` and pass a path to a .fasta file containing unitigs to query in the graph.

By default, Query mode generates files in rtab format. To generate an additional file that is compatible with pyseer, specify `--pyseer`.

The sensitivity of querying can be altered by passing a float argument to `--ratiok`
(between 0 and 1, default 1.0), which determines the threshold proportion of k-mers of a
specific colour present in a unitig for colour classification. Specifying `--inexact` will
search the graph for both exact and inexact k-mers (1 substitution or indel) from queries.
Lowering ''--ratiok' and/or specifying '--inexact' will result in more colour hits per unitig,
but will increase probability of false positives and run-time.


### Option reference
```
usage: unitig-caller [-h] (--build | --query) [--seq SEQ]
                     [--addit_seq ADDIT_SEQ] [--no_colour] [--clean]
                     [--input INPUT] [--ratiok RATIOK] [--inexact] [--pyseer]
                     [--output OUTPUT] [--kmer_size KMER_SIZE]
                     [--minimizer_size MINIMIZER_SIZE] [--threads THREADS]
                     [--bifrost BIFROST] [--version]

Call unitigs in a population dataset

optional arguments:
  -h, --help            show this help message and exit

Mode of operation:
  --build               Build coloured/uncoloured de Bruijn graph using
                        Bifrost
  --query               Query unitig presence/absence across input genomes

Build Input/output:
  --seq SEQ             List of input files in .txt format (reads or
                        references)
  --addit_seq ADDIT_SEQ
                        List of additional input files in .txt format of
                        different type to those in --seq. [default = None]
  --no_colour           Specify for uncoloured de Bruijn Graph [default =
                        False]
  --clean               Clean DBG (clip tips and delete isolated contigs
                        shorter than k k-mers in length) [default = False]

Query Input/output:
  --input INPUT         Prefix for graph and colour files from build exection
                        (file name without extension)
  --ratiok RATIOK       ratio of k-mers from queries that must occur in the
                        graph to be considered as belonging to colour [default
                        = 1.0]
  --inexact             Graph is searched with exact and inexact k-mers (1
                        substitution or indel) from queries[default = False]
  --pyseer              Generate file compatible with pyseer analysis [default
                        = False]

Shared options:
  --output OUTPUT       Prefix for output [default = 'bifrost']
  --kmer_size KMER_SIZE
                        K-mer size for graph building/querying[default = 31]
  --minimizer_size MINIMIZER_SIZE
                        Minimizer size to be used for k-mer hashing [default =
                        23]

Other:
  --threads THREADS     Number of threads to use [default = 1]
  --bifrost BIFROST     Location of bifrost executable [default = Bifrost]
  --version             show program's version number and exit
```

## Citation

If you use this, please cite the Bifrost paper:

Holley G., Melsted, P. Bifrost – Highly parallel construction and indexing of colored and compacted de Bruijn graphs.
bioRxiv 695338 (2019). doi: https://doi.org/10.1101/695338
