# Copyright 2019 John Lees

'''Wrapper around mantis to detect presence of sequence elements'''

import sys

import tempfile
from multiprocessing import Pool
from functools import partial

from .__init__ import __version__

from .squeakr import check_squeakr_version
from .squeakr import squeakr_multi_wrapper

from .mantis import check_mantis_version
from .mantis import run_mantis_build
from .mantis import run_mantis_index
from .mantis import run_mantis_query

def get_options():
    import argparse

    description = 'Call unitigs in a population'
    parser = argparse.ArgumentParser(description=description,
                                     prog='unitig-caller')

    parser.add_argument('--mode',
                        choices=['index', 'call'],
                        required=True,
                        help='\'index\' sequences, or \'call\' on '
                             'an indexed dataset.')

    io = parser.add_argument_group('Input/output')
    io.add_argument('--strains',
                    help='List of strains to index')
    io.add_argument('--unitigs',
                    help='List of unitigs to call')
    io.add_argument('--output',
                    default='mantis_index',
                    help='Prefix for output')
    io.add_argument('--mantis-index',
                    default='mantis_index',
                    help='Directory containing mantis index '
                         '(produced by index mode)')


    other = parser.add_argument_group('Other')
    other.add_argument('--cpus',
                        type=int,
                        default=1,
                        help='Number of CPUs to use.')
    other.add_argument('--overwrite',
                        action='store_true',
                        default=False,
                        help='Overwrite existing output')
    other.add_argument('--squeakr',
                        default='squeakr',
                        help='Location of squeakr executable'
                             '[default = squeakr]')
    other.add_argument('--mantis',
                        default='mantis',
                        help='Location of mantis executable'
                             '[default = mantis]')
    other.add_argument('--version', action='version',
                       version='%(prog)s '+__version__)

    return parser.parse_args()


def main():
    options = get_options()

    mantis_major, mantis_minor, mantis_patch = check_mantis_version(options.mantis)
    if (mantis_major == 0 and mantis_minor < 2):
        sys.stderr.write("Requires mantis version 0.2 or higher")

    if options.mode == "index":
        sys.stderr.write("Creating counts with squeakr\n")

        if (check_squeakr_version(options.squeakr) < 1):
            sys.stderr.write("Requires squeakr version 1.0 or higher")

        squeakr_options = {'overwrite': options.overwrite,
                           'exact': True,
                           'kmer_size': 28,
                           'count_cutoff': 1,
                           'log_slots': 22,
                           'num_threads': 1,
                           'squeakr_exe': options.squeakr}

        strains_in = []
        output_list = tempfile.NamedTemporaryFile('w')
        with open(options.strains, 'r') as strain_file:
            for strain_line in strain_file:
                (strain_name, strain_fasta) = strain_line.rstrip().split("\t")
                strains_in.append((strain_name, strain_fasta))
                output_list.write(strain_name + "/" + strain_name + ".squeakr\n")

        with Pool(processes=options.cpus) as pool:
            pool.map(partial(squeakr_multi_wrapper,
                             overwrite=options.overwrite,
                             squeakr_options=squeakr_options),
                     strains_in)

        output_list.close()

        sys.stderr.write("Building mantis index\n")
        run_mantis_build(output_list, options.prefix, log_slots=22, mantis_exe=options.mantis)
        run_mantis_index(options.prefix, options.cpus, delete_RRR=True, mantis_exe=options.mantis)

    elif options.mode == "call":
        sys.stderr.write("Running mantis queries")

        unitig_list_file = tempfile.NamedTemporaryFile('w')
        unitigs = []
        with open(options.unitigs, 'r') as unitig_file:
            unitig_file.readline() # header
            for unitig_line in unitig_file:
                unitig_fields = unitig_line.rstrip().split("\t")
                unitigs.append(unitig_fields[0])
                unitig_list_file.write(unitig_fields[0] + "\n")

        # Run query and format output as pyseer k-mers/unitigs
        with open(options.prefix + "_unitigs.txt", 'w') as call_output:
            for unitig, samples in zip(unitigs, run_mantis_query(unitig_list_file,
                                                                 options.mantis_index,
                                                                 options.mantis)):
                if len(samples) > 0:
                    out_array = [unitig] + ["|"] + [x + ":1" for x in samples]
                    call_output.write(" ".join(out_array) + "\n")

        unitig_list_file.close()

    sys.exit(0)

if __name__ == "__main__":
    main()