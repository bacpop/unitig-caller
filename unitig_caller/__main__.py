# Copyright 2019-2020 John Lees, Sam Horsfield

'''Detect presence of sequence elements using graphs/indexes'''

import os, sys

import tempfile
from multiprocessing import Pool
from functools import partial

import map_strings

from .__init__ import __version__

from .bifrost import check_bifrost_version
from .bifrost import run_bifrost_build
from .bifrost import gfa_to_fasta
from .bifrost import run_bifrost_query
from .bifrost import rtab_format
from .bifrost import pyseer_format

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
    mode.add_argument('--simple',
                        action='store_true',
                        default=False,
                        help='Use FM-index to make calls ')

    io = parser.add_argument_group('Unitig-caller input/output')
    io.add_argument('--refs',
                    help='Ref file to use to --build bifrost graph (or with --simple)',
                    default=None)
    io.add_argument('--reads',
                    help='Read file to use to --build bifrost graph',
                    default=None)
    io.add_argument('--graph-prefix',
                    help='Prefix of bifrost graph to --query',
                    default=None)
    io.add_argument('--unitigs',
                    help='fasta file of unitigs to query (--query or --simple)',
                    default=None)
    io.add_argument('--output',
                    default='unitig_caller',
                    help='Prefix for output '
                         '[default = \'unitig_caller\']')

    buildio = parser.add_argument_group('Build Input/output')
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
    queryio.add_argument('--ratiok',
                        type=float,
                        default=1.0,
                        help='ratio of k-mers from queries that must occur in the graph to be considered as belonging to colour '
                             '[default = 1.0]')
    queryio.add_argument('--inexact',
                        action='store_true',
                        default=False,
                        help='Graph is searched with exact and inexact k-mers (1 substitution or indel) from queries '
                             '[default = False]')

    bifrost = parser.add_argument_group('Bifrost options')
    bifrost.add_argument('--kmer_size',
                        type=int,
                        default=31,
                        help='K-mer size for graph building/querying '
                             '[default = 31]')
    bifrost.add_argument('--minimizer_size',
                        type=int,
                        default=23,
                        help='Minimizer size to be used for k-mer hashing '
                             '[default = 23]')

    simple = parser.add_argument_group('Simple mode options')
    simple.add_argument('--no-save-idx',
                default=False,
                action='store_true',
                help='Do not save FM-indexes for reuse')

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
    other.add_argument('--version', action='version',
                       version='%(prog)s '+__version__)

    return parser.parse_args()


def main():
    options = get_options()

    if options.build:
        # Read input1 (and input2 if specified) as reads.txt or refs.txt. Call `Bifrost build`

        sys.stderr.write("Building de Bruijn graph with Bifrost\n")
        run_bifrost_build(options.refs,
                          options.output,
                          options.reads,
                          options.no_colour,
                          options.clean,
                          options.kmer_size,
                          options.minimizer_size,
                          options.threads,
                          options.bifrost)

        sys.stderr.write("Creating .fasta file from unitigs in Bifrost graph\n")
        graph_file = options.output + ".gfa"
        gfa_to_fasta(graph_file)

    elif options.query:
        # Read files with prefix input1 as gfa, colours file and fasta file if input2 not specified,
        # or if input2 specified, use this file for unitig querying. Call `Bifrost query`

        graph_file = options.graph_prefix + ".gfa"
        colour_file = options.graph_prefix + ".bfg_colors"
        tsv_file = options.output + ".tsv"
        rtab_file = options.output + ".rtab"

        if options.unitigs == None:
            query_file = options.graph_prefix + "_unitigs.fasta"

        else:
            query_file = options.unitigs

        sys.stderr.write("Querying unitigs in Bifrost graph\n")
        run_bifrost_query(graph_file,
                          query_file,
                          colour_file,
                          options.output,
                          options.ratiok,
                          options.kmer_size,
                          options.minimizer_size,
                          options.threads,
                          options.inexact,
                          options.bifrost)

        sys.stderr.write("Generating rtab file\n")
        rtab_format(tsv_file, rtab_file)

        if os.path.exists(tsv_file):
            os.remove(tsv_file)

        sys.stderr.write("Generating pyseer file\n")
        pyseer_format(rtab_file, query_file)

    elif options.simple:
        # Read input into lists, as in 'index' and 'call'

        if options.refs == None or options.unitigs == None:
            sys.stderr.write("Please specify a strains-list file as input 1 and unitigs file as input 2\n")
        else:
            names_in = []
            fasta_in = []
            with open(options.refs, 'r') as strain_file:
                for strain_line in strain_file:
                    (strain_name, strain_fasta) = strain_line.rstrip().split("\t")
                    names_in.append(strain_name)
                    fasta_in.append(strain_fasta)

            unitigs = []
            with open(options.unitigs, 'r') as unitig_file:
                unitig_file.readline() # header
                for unitig_line in unitig_file:
                    unitig_fields = unitig_line.rstrip().split("\t")
                    unitigs.append(unitig_fields[0])

            # call c++ code to map (also writes output file)
            map_strings.call(fasta_in,
                             names_in,
                             unitigs,
                             options.output + "_pyseer.txt",
                             not options.no_save_idx,
                             options.threads)

    sys.exit(0)

if __name__ == "__main__":
    main()
