# Copyright 2019-2020 John Lees, Sam Horsfield

'''Detect presence of sequence elements using graphs/indexes'''

import os, sys

import tempfile
from multiprocessing import Pool
from functools import partial

import unitig_query

from .__init__ import __version__

from .bifrost import rtab_format
from .bifrost import pyseer_format

def get_options():
    import argparse

    description = 'Call unitigs in a population dataset'
    parser = argparse.ArgumentParser(description=description,
                                     prog='unitig-caller')

    modeGroup = parser.add_argument_group('Mode of operation')
    mode = modeGroup.add_mutually_exclusive_group(required=True)
    mode.add_argument('--call',
                        action='store_true',
                        default=False,
                        help='Build a DBG and call colours of unitigs within ')
    mode.add_argument('--query',
                        action='store_true',
                        default=False,
                        help='Query unitig colours in reference genomes/DBG ')
    mode.add_argument('--simple',
                        action='store_true',
                        default=False,
                        help='Use FM-index to make calls ')

    io = parser.add_argument_group('Unitig-caller input/output')
    io.add_argument('--refs',
                    help='Ref file to used to build DBG or use with --simple',
                    default=None)
    io.add_argument('--reads',
                    help='Read file to used to build DBG',
                    default=None)
    io.add_argument('--graph',
                    help='Existing graph in GFA format',
                    default=None)
    io.add_argument('--colours',
                    help='Existing bifrost colours file in .bfg_colors format',
                    default=None)
    io.add_argument('--unitigs',
                    help='Text or fasta file of unitigs to query (--query or --simple)',
                    default=None)
    io.add_argument('--pyseer',
                    action='store_true',
                    help='Output pyseer format',
                    default=False)
    io.add_argument('--rtab',
                    action='store_true',
                    help='Output rtab format',
                    default=False)
    io.add_argument('--out',
                    default='unitig_caller',
                    help='Prefix for output '
                         '[default = \'unitig_caller\']')

    bifrost = parser.add_argument_group('Bifrost options')
    bifrost.add_argument('--kmer',
                        type=int,
                        default=31,
                        help='K-mer size for graph building/querying '
                             '[default = 31]')
    bifrost.add_argument('--write-graph',
                        action='store_true',
                        default=False,
                        help='Output DBG built with unitig-caller')

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
    other.add_argument('--version', action='version',
                       version='%(prog)s '+__version__)

    return parser.parse_args()


def main():
    options = get_options()

    if options.call or options.query:
        if options.call and options.graph != None and options.colours != None and options.refs == None and options.reads == None:
            # Read input1 (and input2 if specified) as reads.txt or refs.txt. Call `Bifrost build`

            sys.stderr.write("Calling unitigs within input genomes...\n")
            unitig_map, input_colour_pref = unitig_query.call_unitigs_existing(options.graph, options.colours, True, "NA", options.threads)

        elif options.call and (options.refs != None or options.reads != None) and options.graph == None and options.colours == None:

            sys.stderr.write("Building a DBG and calling unitigs within...\n")
            if options.refs != None and options.reads == None:
                unitig_map, input_colour_pref = unitig_query.call_unitigs_build(options.refs, options.kmer, True, "NA", options.threads, True, options.write_graph)
            elif options.refs == None and options.reads != None:
                unitig_map, input_colour_pref = unitig_query.call_unitigs_build(options.reads, options.kmer, True, "NA", options.threads, False, options.write_graph)
            elif options.refs != None and options.reads != None:
                unitig_map, input_colour_pref = unitig_query.call_unitigs_build(options.refs, options.kmer, True, "NA", options.threads, False, options.write_graph, options.reads)

        elif options.query and options.graph != None and options.colours != None and options.unitigs != None and options.refs == None and options.reads == None:

            sys.stderr.write("Querying unitigs within existing DBG...\n")
            unitig_map, input_colour_pref = unitig_query.call_unitigs_existing(options.graph, options.colours, False, options.unitigs, options.threads)

        elif options.query and (options.refs != None or options.reads != None) and options.unitigs != None and options.graph == None and options.colours == None:
            sys.stderr.write("Building a DBG and querying unitigs within...\n")
            if options.refs != None and options.reads == None:
                unitig_map, input_colour_pref = unitig_query.call_unitigs_build(options.refs, options.kmer, False, options.unitigs, options.threads, True, options.write_graph)
            elif options.refs == None and options.reads != None:
                unitig_map, input_colour_pref = unitig_query.call_unitigs_build(options.reads, options.kmer, False, options.unitigs, options.threads, False, options.write_graph)
            elif options.refs != None and options.reads != None:
                unitig_map, input_colour_pref = unitig_query.call_unitigs_build(options.refs, options.kmer, False, options.unitigs, options.threads, False, options.write_graph, options.reads)

        else:
            print("Error: incorrect number of input files specified. Please only specify the below combinations:\n"
                "For --call:\n"
                "   - Bifrost GFA and Bifrost colours file (--graph & --colours)\n"
                "   - List of reference files (--refs)\n"
                "   - List of read files (--reads)\n"
                "   - A list of reference files and a list of read files (--refs & --reads)\n"
                "For --query:\n"
                "   - One of the above combinations with a text file with header, one file per line/fasta file of query unitigs (--unitigs)\n")
            sys.exit(1)

        if options.pyseer or (not options.pyseer and not options.rtab):
            sys.stderr.write("Generating pyseer file...\n")
            pyseer_format(unitig_map, input_colour_pref, options.out + ".pyseer")

        if options.rtab:
            sys.stderr.write("Generating rtab file...\n")
            rtab_format(unitig_map, input_colour_pref, options.out + ".rtab")

    elif options.simple:
        # Read input into lists, as in 'index' and 'call'

        if options.refs == None or options.unitigs == None:
            sys.stderr.write("Error: Please specify a list of assemblies as --refs and unitigs file as --unitigs\n")
            sys.exit(1)
        else:
            names_in = []
            fasta_in = []
            with open(options.refs, 'r') as strain_file:
                for strain_line in strain_file:
                    strain_line = strain_line.rstrip()
                    names_in.append(os.path.splitext(os.path.basename(strain_line))[0])
                    fasta_in.append(strain_line)

            unitigs = []
            with open(options.unitigs, 'r') as unitig_file:
                unitig_file.readline() # header
                for unitig_line in unitig_file:
                    if unitig_line[0] != '>':
                        unitig_fields = unitig_line.rstrip().split("\t")
                        unitigs.append(unitig_fields[0])

            # call c++ code to map (also writes output file)
            unitig_query.call(fasta_in,
                             names_in,
                             unitigs,
                             options.out + ".pyseer",
                             not options.no_save_idx,
                             options.threads)

    else:
        sys.stderr.write("Error: Please specify one of: build, query or simple modes. Use -h or --help for more information.\n")
        sys.exit(1)

    sys.exit(0)

if __name__ == "__main__":
    main()
