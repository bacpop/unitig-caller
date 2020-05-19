#!/usr/bin/env python

"""Tests for unitig-caller"""

import subprocess
import os
import sys


sys.stderr.write("Testing bifrost build\n")
subprocess.run("unitig-caller --build --refs refs.txt --output test_build", shell=True, check=True)

sys.stderr.write("Testing bifrost query\n")
subprocess.run("unitig-caller --query --graph-prefix test_build --unitigs test_unitigs.fasta --output test_query", shell=True, check=True) # checks if can be run
subprocess.run("python test-calls.py --method rtab --test test_query.rtab --expected bifrost_result.rtab", shell=True, check=True) # checks results match
subprocess.run("python test-calls.py --method pyseer --test test_query.pyseer --expected bifrost_result.pyseer", shell=True, check=True) # checks results match

sys.stderr.write("Testing simple mode\n")
subprocess.run("unitig-caller ---simple --refs simple_refs.txt --unitigs test_unitigs_simple.fasta --output simple_calls", shell=True, check=True) # checks if can be run
subprocess.run("python test-calls.py --method pyseer --test simple_calls_pyseer.txt --expected simple_results.pyseer", shell=True, check=True) # checks results match

sys.stderr.write("Tests completed\n")

