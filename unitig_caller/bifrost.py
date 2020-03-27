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

def run_bifrost_build(in_file, add_in_file = None, out_file, coloured = True, kmer_size=31, threads=1, bifrost_exe='Bifrost'):
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
            k-mer size used for DBG construction (maximum 31)
            [default = 31]
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
    bifrost_cmd += " -t " + str(threads)
    if coloured == True:
        bifrost_cmd += " -c"

    if any(item in str(bifrost_cmd) for item in [' -r ', ' -s ']:
        subprocess.run(bifrost_cmd, shell=True, check=True)
    else:
        return "Please submit input files as 'reads.txt' or 'refs.txt' only"

def run_mantis_index(mantis_dir, num_threads=1, delete_RRR=False, mantis_exe='mantis'):
    """Runs mantis mst on a built index.

    Args:
        mantis_dir (str)
            Location of input directory (-p)
        num_threads (int)
            Number of threads to use when counting (-t)

            [default = 1]
        delete_RRR (bool)
            Remove the previous color class RRR representation (-d)

            [default = False]
        mantis_exe (str)
            Location of mantis executable

            [deafault = 'mantis']
    """
    mantis_cmd = mantis_exe + " mst"
    mantis_cmd += " -p " + mantis_dir + "/"
    mantis_cmd += " -t " + str(num_threads)

    # Note, delete will fail without trailing slash
    if delete_RRR:
        mantis_cmd += " -d"
    else:
        mantis_cmd += " -k"

    subprocess.run(mantis_cmd, shell=True, check=True)

def run_mantis_query(query_file, mantis_index, mantis_exe='mantis'):
    """Runs mantis call. Returns iterable with strains
    for each input sequence.

    Args:
        query_file (str)
            Location of input sequences
        mantis_index (str)
            Location of index to search
        mantis_exe (str)
            Location of mantis executable

            [deafault = 'mantis']
    """
    query_out_file = tempfile.mkstemp()[1]

    mantis_cmd = mantis_exe + " query"
    mantis_cmd += " -p " + mantis_index
    mantis_cmd += " -o " + query_out_file
    mantis_cmd += " " + query_file

    subprocess.run(mantis_cmd, shell=True, check=True)

    max_match = 0
    samples = []
    with open(query_out_file, 'r') as query_file:
        for query_line in query_file:
            new_seq = re.search(r"^seq\d+\t(\d+)", query_line.rstrip())
            if new_seq:
                if max_match > 0:
                    yield(samples)
                    samples = []
                max_match = int(new_seq.group(1))
            else:
                (squeakr_file, matches) = query_line.rstrip().split("\t")
                if int(matches) == max_match:
                    sample_name = re.search(r"\/(.+?)\.squeakr$", squeakr_file)
                    if sample_name:
                        samples.append(sample_name.group(1))
                    else:
                        sys.stderr.write("Error matching sample name in " + squeakr_file + "\n")

    yield(samples)
    os.remove(query_out_file)
