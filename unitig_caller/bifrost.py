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

def run_bifrost_build(in_file, out_file, addit_in_file = None, no_colour = False, clean = False, kmer_size = 31, minimizer_size = 23, threads = 1, bifrost_exe='Bifrost'):
    """Runs Bifrost build on set of refernce genomes/read files.

    Args:
        in_file (str)
            path for .txt file containing list of paths to fasta/fastq files. Must be specified as either 'reads.txt' or 'refs.txt'
        addit_in_file (str)
            additional .txt files containing list of paths to fasta/fastq files. Must be specified as either 'reads.txt' or 'refs.txt'
        out_file (str)
            prefix for output files
        no_colour (bool)
            Do not generate coloured compacted DBG
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

    if addit_in_file == None:
        if "reads.txt" in str(in_file):
            bifrost_cmd += " -s " + str(in_file)
        elif "refs.txt" in str(in_file):
            bifrost_cmd += " -r " + str(in_file)
        else:
            pass
    else:
        if "reads.txt" in str(in_file) and "refs.txt" in str(addit_in_file):
            bifrost_cmd += " -s " + str(in_file)
            bifrost_cmd += " -r " + str(addit_in_file)
        elif "refs.txt" in str(in_file) and "reads.txt" in str(addit_in_file):
            bifrost_cmd += " -r " + str(in_file)
            bifrost_cmd += " -s " + str(addit_in_file)
        else:
            pass
    bifrost_cmd += " -o " + out_file
    bifrost_cmd += " -k " + str(kmer_size)
    bifrost_cmd += " -m " + str(minimizer_size)
    bifrost_cmd += " -t " + str(threads)
    if no_colour == False:
        bifrost_cmd += " -c"
    if clean == True:
        bifrost_cmd += " -i -d "
    print(bifrost_cmd)
    if any(item in str(bifrost_cmd) for item in [' -r ', ' -s ']):
        subprocess.run(bifrost_cmd, shell=True, check=True)
    else:
        return "Please submit input files as 'reads.txt' or 'refs.txt' only"

def gfa_to_fasta(in_file):
    """Converts .gfa file to .fasta, using same file name.

    Args:
        in_file (str)
            path for .gfa file to convert.
    """
    base = os.path.splitext(in_file)[0]
    out_file = base + ".fasta"
    with open(in_file, "r") as f, open(out_file, "w") as o:
        for line in f:
            parsed_line = re.split(r'\t+', line.rstrip('\t'))
            if parsed_line[0] == "S":
                    o.write(">" + str(parsed_line[1]) + "\n" + str(parsed_line[2]) + "\n")

def run_bifrost_query(graph_file, query_file, colour_file, out_file, ratio_k = 0.8, kmer_size = 31, minimizer_size = 23, threads = 1, inexact = False, bifrost_exe='Bifrost'):
    """Runs query of unitigs on coloured Bifrost graph.

    Args:
        graph_file (str)
            input .gfa file from "Bifrost build" execution
        query_file (str)
            .fasta file/.txt file containing sequences to be queried
        colour_file (str)
            input colour file from "Bifrost build" execution
        out_file (str)
            prefix for output file
        ratio_k (float)
            ratio of k-mers from queries that must occur in the graph
            [default = 0.8]
        kmer_size (int)
            k-mer size used for query search
            [default = 31]
        minimizer_size (int)
            minimizer size used for k-mer hashing (maximum = kmer_size)
            [default = 21]
        threads (int)
            number of threads to use
            [default = 1]
        inexact (bool)
            Graph is searched with exact and inexact k-mers (1 substitution or indel) from queries
            [default = False]
        bifrost_exe (str)
            Location of bifrost executable
            [default = 'bifrost']
    """
    bifrost_cmd = bifrost_exe + " query"

    bifrost_cmd += " -g " + graph_file
    bifrost_cmd += " -q " + query_file
    bifrost_cmd += " -o " + out_file
    bifrost_cmd += " -e " + str(ratio_k)
    bifrost_cmd += " -f " + colour_file
    bifrost_cmd += " -k " + str(kmer_size)
    bifrost_cmd += " -m " + str(minimizer_size)
    bifrost_cmd += " -t " + str(threads)

    if inexact == True:
        bifrost_cmd += " -n"
    print(bifrost_cmd)
    subprocess.run(bifrost_cmd, shell=True, check=True)
