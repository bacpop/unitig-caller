#!/usr/bin/env python
# Copyright 2018 John Lees and Nick Croucher

"""Clean test files"""

import os
import sys
import shutil

def deleteDir(dirname):
    if os.path.isdir(dirname):
        shutil.rmtree(dirname)    

sys.stderr.write("Cleaning up tests\n")

# clean up
outputFiles = [
    "12673_8#24.contigs_velvet.fa.fm",
    "12673_8#26.contigs_velvet.fa.fm",
    "12673_8#27.contigs_velvet.fa.fm",
    "12673_8#28.contigs_velvet.fa.fm",
    "12673_8#29.contigs_velvet.fa.fm",
    "test_build.gfa",
    "test_build.bfg_colors",
    "test_build_unitigs.fasta",
    "test_query.rtab",
    "test_query_pyseer.txt",
    "test_unitigs_simple.fasta",
    "simple_calls_pyseer.txt"
]

for delFile in outputFiles:
    if os.path.isfile(delFile):
        os.remove(delFile)


