# Copyright 2019 John Lees

'''Functions to run mantis'''

import os, sys
import subprocess
import re
import tempfile

def check_mantis_version(exe='mantis'):
    """Checks that mantis can be run, and its version.

    Args:
        exe (str)
            Location of executable
    Returns:
        version (tuple)
            Tuple of major, minor, patch mantis version
    """
    p = subprocess.Popen([exe + ' -v'], shell=True, stdout=subprocess.PIPE)
    major_version = 0
    minor_version = 0
    patch_version = 0
    for line in iter(p.stdout.readline, ''):
        if line != '':
            version_match = re.search(r"^mantis (\d+)\.(\d+)\.(\d+)$", line.rstrip().decode())
            major_version = int(version_match.group(1))
            minor_version = int(version_match.group(2))
            patch_version = int(version_match.group(3))
            break

    return (major_version, minor_version, patch_version)

def run_mantis_build(squeakr_files, out_dir, log_slots=22, mantis_exe='mantis'):
    """Runs mantis build on a set of squeakr files.

    Args:
        squeakr_files (str)
            File listing input squeakr file
        out_dir (str)
            Output directory
        log_slots (int)
            Log of number of count slots (-s)

            [default = 22]
        mantis_exe (str)
            Location of mantis executable

            [deafault = 'mantis']
    """
    mantis_cmd = mantis_exe + " build"
    mantis_cmd += " -s " + str(log_slots)
    mantis_cmd += " -i " + squeakr_files
    mantis_cmd += " -o " + out_dir

    subprocess.run(mantis_cmd, shell=True, check=True)

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
            new_seq = re.search(r"^seq\d\t(\d+)", query_line.rstrip())
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
