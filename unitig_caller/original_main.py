# Copyright 2019 John Lees

'''Wrapper around mantis to detect presence of sequence elements'''

import os, sys

import tempfile
from multiprocessing import Pool
from functools import partial

import map_strings

from .__init__ import __version__

from .squeakr import check_squeakr_version
from .squeakr import squeakr_multi_wrapper

from .mantis import check_mantis_version
from .mantis import run_mantis_build
from .mantis import run_mantis_index
from .mantis import run_mantis_query

def get_options():
    import argparse

    description = 'Call unitigs in a population dataset'
    parser = argparse.ArgumentParser(description=description,
                                     prog='unitig-caller')

    modeGroup = parser.add_argument_group('Mode of operation')
    mode = modeGroup.add_mutually_exclusive_group(required=True)
    mode.add_argument('--index',
                        action='store_true',
                        default=False,
                        help='Index sequences, before calling')
    mode.add_argument('--call',
                        action='store_true',
                        default=False,
                        help='Make calls from an indexed dataset')
    mode.add_argument('--simple',
                        action='store_true',
                        default=False,
                        help='Use FM-index to make calls')

    io = parser.add_argument_group('Input/output')
    io.add_argument('--strains',
                    help='List of strains to index')
    io.add_argument('--unitigs',
                    help='List of unitigs to call')
    io.add_argument('--output',
                    default='mantis_index',
                    help='Prefix for output '
                         '[default = \'mantis_index\']')
    io.add_argument('--no-save-idx',
                    default=False,
                    action='store_true',
                    help='Do not save FM-indexes for reuse')
    io.add_argument('--mantis-index',
                    default='mantis_index',
                    help='Directory containing mantis index '
                         '(produced by index mode) '
                         '[default = \'mantis_index\']')

    mantis = parser.add_argument_group('mantis/squeakr options')
    mantis.add_argument('--approximate',
                        action='store_true',
                        default=False,
                        help='Use approximate count mode '
                             '[default = exact]')
    mantis.add_argument('--kmer-size',
                        type=int,
                        default=31,
                        help='K-mer size for counts '
                             '[default = 31]')
    mantis.add_argument('--count-cutoff',
                        type=int,
                        default=1,
                        help='Minimum k-mer count to be included '
                             '[default = 1]')
    mantis.add_argument('--log-slots',
                        type=int,
                        default=22,
                        help='Starting log(number of count slots). '
                             'Automatically increased if necessary '
                             '[default = 22]')

    other = parser.add_argument_group('Other')
    other.add_argument('--cpus',
                        type=int,
                        default=1,
                        help='Number of CPUs to use '
                             '[default = 1]')
    other.add_argument('--overwrite',
                        action='store_true',
                        default=False,
                        help='Overwrite existing output')
    other.add_argument('--squeakr',
                        default='squeakr',
                        help='Location of squeakr executable '
                             '[default = squeakr]')
    other.add_argument('--mantis',
                        default='mantis',
                        help='Location of mantis executable '
                             '[default = mantis]')
    other.add_argument('--version', action='version',
                       version='%(prog)s '+__version__)

    return parser.parse_args()


def main():
    options = get_options()

    if not options.simple:
        mantis_major, mantis_minor, mantis_patch = check_mantis_version(options.mantis)
        if (mantis_major == 0 and mantis_minor < 2):
            sys.stderr.write("Requires mantis version 0.2 or higher\n")

    if options.index:
        sys.stderr.write("Creating counts with squeakr\n")

        if (check_squeakr_version(options.squeakr) < 1):
            sys.stderr.write("Requires squeakr version 1.0 or higher\n")

        squeakr_options = {'overwrite': options.overwrite,
                           'exact': not options.approximate,
                           'kmer_size': options.kmer_size,
                           'count_cutoff': options.count_cutoff,
                           'log_slots': options.log_slots,
                           'num_threads': 1, # a pool is used instead
                           'squeakr_exe': options.squeakr}

        strains_in = []
        output_list = tempfile.NamedTemporaryFile('w')
        with open(options.strains, 'r') as strain_file:
            for strain_line in strain_file:
                (strain_name, strain_fasta) = strain_line.rstrip().split("\t")
                strains_in.append((strain_name, strain_fasta))
                output_list.write(strain_name + "/" + strain_name + ".squeakr\n")
        output_list.file.flush()

        with Pool(processes=options.cpus) as pool:
            pool.map(partial(squeakr_multi_wrapper,
                             overwrite=options.overwrite,
                             squeakr_options=squeakr_options),
                     strains_in)

        sys.stderr.write("Building mantis index\n")
        run_mantis_build(output_list.name, options.output, log_slots=22, mantis_exe=options.mantis)
        run_mantis_index(options.output, options.cpus, delete_RRR=False, mantis_exe=options.mantis)

        output_list.close()

    elif options.call:
        sys.stderr.write("Running mantis queries\n")

        unitig_list_file = tempfile.NamedTemporaryFile('w')
        unitigs = []
        with open(options.unitigs, 'r') as unitig_file:
            unitig_file.readline() # header
            for unitig_line in unitig_file:
                unitig_fields = unitig_line.rstrip().split("\t")
                unitigs.append(unitig_fields[0])
                unitig_list_file.write(unitig_fields[0] + "\n")
        unitig_list_file.file.flush()

        # Run query and format output as pyseer k-mers/unitigs
        mantis_index = os.path.join(options.mantis_index, '')
        with open(options.output + "_unitigs.txt", 'w') as call_output:
            for unitig, samples in zip(unitigs, run_mantis_query(unitig_list_file.name,
                                                                 mantis_index,
                                                                 options.mantis)):
                if len(samples) > 0:
                    out_array = [unitig] + ["|"] + [x + ":1" for x in samples]
                    call_output.write(" ".join(out_array) + "\n")

        unitig_list_file.close()

    elif options.simple:
        # Read input into lists, as in 'index' and 'call'
        names_in = []
        fasta_in = []
        with open(options.strains, 'r') as strain_file:
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
                         options.output + "_unitigs.txt",
                         not options.no_save_idx,
                         options.cpus)

    sys.exit(0)

if __name__ == "__main__":
    main()
