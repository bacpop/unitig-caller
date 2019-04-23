# Copyright 2019 John Lees

'''Functions to run mantis'''

import os, sys
import subprocess
import re

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
            major_version = version_match.group(1)
            minor_version = version_match.group(2)
            patch_version = version_match.group(3)
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

def run_mantis_index(mantis_dir, num_threads=1, delete_RRR=True, mantis_exe='mantis'):
    """Runs mantis mst a built index.

    Args:
        mantis_dir (str)
            Location of input directory (-p)
        num_threads (int)
            Number of threads to use when counting (-t)

            [default = 1]
        delete_RRR (bool)
            Remove the previous color class RRR representation (-d)

            [default = True]
        mantis_exe (str)
            Location of mantis executable

            [deafault = 'mantis']
    """
    mantis_cmd = mantis_exe + " mst"
    mantis_cmd += " -p " + mantis_dir
    mantis_cmd += " -t " + str(num_threads)
    if delete_RRR:
        mantis_cmd += " -d"
    else:
        mantis_cmd += " -k"

    subprocess.run(mantis_cmd, shell=True, check=True)
