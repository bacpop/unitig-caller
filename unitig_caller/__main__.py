# Copyright 2019 John Lees

'''Wrapper around Bifrost to detect presence of sequence elements'''

import os, sys

import tempfile
from multiprocessing import Pool
from functools import partial

#import map_strings

#from .__init__ import __version__

#from .bifrost import check_bifrost_version
#from .bifrost import run_bifrost_build
#from .bifrost import gfa_to_fasta
#from .bifrost import run_bifrost_query

from __init__ import __version__

from bifrost import check_bifrost_version
from bifrost import run_bifrost_build
from bifrost import gfa_to_fasta
from bifrost import run_bifrost_query
from bifrost import rtab_format
from bifrost import pyseer_format

def get_options():
    import argparse

    description = 'Call unitigs in a population dataset'
    parser = argparse.ArgumentParser(description=description,
                                     prog='unitig-caller')

    modeGroup = parser.add_argument_group('Mode of operation')
    mode = modeGroup.add_mutually_exclusive_group(required=True)
    mode.add_argument('--build',
                        action='store_true',
                        default=False,
                        help='Build coloured/uncoloured de Bruijn graph using Bifrost ')
    mode.add_argument('--query',
                        action='store_true',
                        default=False,
                        help='Query unitig presence/absence across input genomes ')

    buildio = parser.add_argument_group('Build Input/output')
    buildio.add_argument('--seq',
                    help='List of input files in .txt format (reads or references) ')
    buildio.add_argument('--addit_seq',
                    default=None,
                    help='List of additional input files in .txt format of different type to those in --seq. '
                         '[default = None]')
    buildio.add_argument('--no_colour',
                        action='store_true',
                        default=False,
                        help='Specify for uncoloured de Bruijn Graph '
                             '[default = False]')
    buildio.add_argument('--clean',
                        action='store_true',
                        default=False,
                        help='Clean DBG (clip tips and delete isolated contigs shorter than k k-mers in length) '
                             '[default = False]')

    queryio = parser.add_argument_group('Query Input/output')
    queryio.add_argument('--input',
                    help='Prefix for graph and colour files from build exection (file name without extension)  ')
    queryio.add_argument('--ratiok',
                        type=float,
                        default=0.8,
                        help='ratio of k-mers from queries that must occur in the graph to be considered as belonging to colour'
                             '[default = 0.8]')
    queryio.add_argument('--inexact',
                        action='store_true',
                        default=False,
                        help='Graph is searched with exact and inexact k-mers (1 substitution or indel) from queries'
                             '[default = False]')
    queryio.add_argument('--pyseer',
                        action='store_true',
                        default=False,
                        help='Generate file compatible with pyseer analysis.'
                             '[default = False]')

    shared = parser.add_argument_group('Shared options')
    shared.add_argument('--output',
                    default='bifrost',
                    help='Prefix for output '
                         '[default = \'bifrost\']')
    shared.add_argument('--kmer_size',
                        type=int,
                        default=31,
                        help='K-mer size for graph building/querying'
                             '[default = 31]')
    shared.add_argument('--minimizer_size',
                        type=int,
                        default=23,
                        help='Minimizer size to be used for k-mer hashing '
                             '[default = 23]')

    other = parser.add_argument_group('Other')
    other.add_argument('--threads',
                        type=int,
                        default=1,
                        help='Number of threads to use '
                             '[default = 1]')
    other.add_argument('--bifrost',
                        default='Bifrost',
                        help='Location of bifrost executable '
                             '[default = Bifrost]')
    #other.add_argument('--version', action='version',
    #                   version='%(prog)s '+__version__)

    return parser.parse_args()


def main():
    options = get_options()

    if options.build:
        sys.stderr.write("Building de Bruijn graph with Bifrost\n")

        run_bifrost_build(options.seq, options.output, options.addit_seq, options.no_colour, options.clean, options.kmer_size, options.minimizer_size, options.threads, options.bifrost)

    elif options.query:

        sys.stderr.write("Creating .fasta query file from unitigs in Bifrost graph\n")

        graph_file = options.input + ".gfa"
        query_file = options.input + "_unitigs.fasta"
        colour_file = options.input+ ".bfg_colors"
        tsv_file = options.output + ".tsv"

        gfa_to_fasta(graph_file)

        sys.stderr.write("Querying unitigs in Bifrost graph\n")

        run_bifrost_query(graph_file, query_file, colour_file, options.output, options.ratiok, options.kmer_size, options.minimizer_size, options.threads, options.inexact, options.bifrost)

        sys.stderr.write("Generating rtab file\n")

        rtab_format(tsv_file)

        if options.pyseer == True:
            sys.stderr.write("Generating pyseer file\n")
            pyseer_format(tsv_file, query_file)

    sys.exit(0)

if __name__ == "__main__":
    main()
