# Copyright 2019 John Lees

'''Wrapper around Bifrost to detect presence of sequence elements'''

import os, sys

import tempfile
from multiprocessing import Pool
from functools import partial

#import map_strings

from .__init__ import __version__

from .bifrost import check_bifrost_version
from .bifrost import run_bifrost_build

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

    io = parser.add_argument_group('Input/output')
    io.add_argument('--input',
                    help='List of input files (reads or references) ')
    io.add_argument('--addit_input',
                    default=None,
                    help='List of additional input files of different type to --input '
                         '[default = None]')
    io.add_argument('--output',
                    default='bifrost_graph',
                    help='Prefix for output '
                         '[default = \'bifrost_graph\']')

    bifrost = parser.add_argument_group('bifrost options')
    bifrost.add_argument('--coloured',
                        default=True,
                        help='Specify for coloured/uncoloured de Bruijn Graph '
                             '[default = True]')
    bifrost.add_argument('--kmer_size',
                        type=int,
                        default=31,
                        help='K-mer size for graph building'
                             '[default = 31]')
    bifrost.add_argument('--minimizer_size',
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
        sys.stderr.write("Building de Bruijn Graph with Bifrost\n")

        run_bifrost_build(options.input, options.output, options.addit_input, options.coloured, options.kmer_size, options.minimizer_size, options.threads, options.bifrost)

    sys.exit(0)

if __name__ == "__main__":
    main()
