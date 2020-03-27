# Copyright 2019 John Lees

'''Functions to run Bifrost'''

import os, sys
import subprocess
import re
import tempfile

def check_bifrost_version(exe='Bifrost'):
    """Checks that Bifrost can be run, and it's version.

    Args:
        exe (str)
            Location of executable
    Returns:
        version (tuple)
            Tuple of major, minor, patch bifrost version
    """
    p = subprocess.Popen([exe], shell=True, stdout=subprocess.PIPE)
    major_version = 0
    minor_version = 0
    patch_version = 0
    install_error = "No Bifrost install available"
    version_line = p.stdout.readline()
    version_match = re.search(r"^Bifrost (\d+).(\d+).(\d+)$", version_line.rstrip().decode())
    major_version = int(version_match.group(1))
    minor_version = int(version_match.group(2))
    patch_version = int(version_match.group(3))
    if major_version != 0 or minor_version != 0 or patch_version !=0:
        return (major_version, minor_version, patch_version)
    else:
        return install_error

def run_bifrost_build(in_file, add_in_file = None, out_file, coloured = True, kmer_size=31, minimizer_size = 23, threads=1, bifrost_exe='Bifrost'):
    """Runs mantis build on a set of squeakr files.

    Args:
        in_file (str)
            path for .txt file containing list of paths to fasta/fastq files. Must be specified as either 'reads.txt' or 'refs.txt'
        add_in_file (str)
            additional .txt files containing list of paths to fasta/fastq files. Must be specified as either 'reads.txt' or 'refs.txt'
        out_file (str)
            prefix for output files
        coloured (bool)
            Colour compacted DBG
            [default = True]
        kmer_size (int)
            k-mer size used for DBG construction (maximum = 31)
            [default = 31]
        minimizer_size (int)
            minimizer size used for k-mer hashing (maximum = kmer_size)
            [default = 21]
        threads (int)
            number of threads to use
            [default = 1]
        bifrost_exe (str)
            Location of bifrost executable
            [default = 'bifrost']
    """
    bifrost_cmd = bifrost_exe + " build"

    if add_in_file == None:
        if "reads.txt" in str(in_file):
            bifrost_cmd += " -s " + str(in_file)
        elif "refs.txt" in str(in_file):
            bifrost_cmd += " -r " + str(in_file)
        else:
            pass
    else:
        if "reads.txt" in str(in_file) and "refs.txt" in str(add_in_file):
            bifrost_cmd += " -s " + str(in_file)
            bifrost_cmd += " -r " + str(add_in_file)
        elif "refs.txt" in str(in_file) and "reads.txt" in str(add_in_file):
            bifrost_cmd += " -r " + str(in_file)
            bifrost_cmd += " -s " + str(add_in_file)
        else:
            pass
    bifrost_cmd += " -o " + out_file
    bifrost_cmd += " -k " + str(kmer_size)
    bifrost_cmd += " -m " + str(minimizer_size)
    bifrost_cmd += " -t " + str(threads)
    if coloured == True:
        bifrost_cmd += " -c"

    if any(item in str(bifrost_cmd) for item in [' -r ', ' -s ']:
        subprocess.run(bifrost_cmd, shell=True, check=True)
    else:
        return "Please submit input files as 'reads.txt' or 'refs.txt' only"
