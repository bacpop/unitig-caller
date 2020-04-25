# Copyright 2019 John Lees

'''Functions to run squeakr'''

import os, sys
import subprocess
import re

def check_squeakr_version(exe='squeakr'):
    """Checks that squeaker can be run, and its version.

    Args:
        exe (str)
            Location of executable
    Returns:
        version (int)
            Major version of squeakr. Zero if failure
    """
    p = subprocess.Popen([exe + ' help -v'], shell=True, stdout=subprocess.PIPE)
    version = 0
    for line in iter(p.stdout.readline, ''):
        if line != '':
            version_match = re.search(r"^version (\d+)\.(\d+)$", line.rstrip().decode())
            version = int(version_match.group(1))
            break

    return version

def fasta_to_fastq(fasta_in, fastq_out, fake_qual="I"):
    """Converts a fasta to a fastq for use with squeakr

    Args:
        fasta_in (str)
            Location of fasta input file
        fastq_out (str)
            Fastq output file
        fake_qual (str)
            Quality score to pad with (ignored by squeakr)
    """

    # Could use biopython, but probably not worth it
    with open(fasta_in, 'r') as seq_in, open(fastq_out, 'w') as seq_out:
        sequence = ""
        for fasta_line in seq_in:
            header = re.search(r"^>(.+)$", fasta_line.rstrip())
            if header is not None and sequence != "":
                seq_out.write("@" + header.group(1) + "\n")
                seq_out.write(sequence + "\n")
                seq_out.write("+\n")
                seq_out.write(fake_qual * len(sequence) + "\n")
                sequence = ""
            else:
                sequence += fasta_line.rstrip()

def run_squeakr(fastq_in, out_file, exact=True, kmer_size=28, count_cutoff=1,
                log_slots=22, num_threads=1, squeakr_exe='squeakr'):
    """Runs squeakr on a single input file.

    Args:
        fastq_in (str)
            Location of input fastq file
        out_file (str)
            Output squeakr file
        exact (bool)
            Run in exact mode (-e) (as opposed to approximate)

            [default = True]
        kmer_size (int)
            K-mer size to count (-k)

            [default = 28]
        count_cutoff (int)
            Minimum number of k-mers to add to count (-c)

            [default = 1]
        log_slots (int)
            Log of number of count slots (-s)

            [default = 22]
        num_threads (int)
            Number of threads to use when counting (-t)

            [default = 1]
        squeakr_exe (str)
            Location of squeakr executable

            [deafault = 'squeakr']
    """
    squeakr_cmd = squeakr_exe + " count"
    squeakr_cmd += " -k " + str(kmer_size)
    squeakr_cmd += " -s " + str(log_slots)
    squeakr_cmd += " -t " + str(num_threads)
    squeakr_cmd += " -c " + str(count_cutoff)
    if exact:
        squeakr_cmd += " -e"
    squeakr_cmd += " -o " + out_file + " " + fastq_in

    subprocess.run(squeakr_cmd,
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL,
                   shell=True,
                   check=True)

def squeakr_multi_wrapper(names, overwrite, squeakr_options):
    """Turns fasta into fastq, runs squeakr, cleans up tmp files.

    Args:
        names (tuple)
            Tuple of strings (Name for output, Location of input file)
        overwrite (bool)
            Overwrite output
        squeakr_options (dict)
            Dictionary of options for :func:`~run_squeakr`
    Returns:
        version (int)
            Major version of squeakr. Zero if failure
    """
    output_name, fasta_in = names
    tmp_fastq_file = output_name + "/" + output_name + ".fastq"
    output_file = output_name + "/" + output_name + ".squeakr"
    if overwrite or not os.path.isfile(output_file):
        if not os.path.isdir(output_name):
            try:
                os.makedirs(output_name)
            except OSError:
                sys.stderr.write("Cannot create output directory " + output_name + "\n")
                sys.exit(1)

        fasta_to_fastq(fasta_in, tmp_fastq_file)
        run_squeakr(tmp_fastq_file,
                    output_file,
                    exact=squeakr_options['exact'],
                    kmer_size=squeakr_options['kmer_size'],
                    count_cutoff=squeakr_options['count_cutoff'],
                    log_slots=squeakr_options['log_slots'],
                    num_threads=squeakr_options['num_threads'],
                    squeakr_exe=squeakr_options['squeakr_exe'])
        os.remove(tmp_fastq_file)
