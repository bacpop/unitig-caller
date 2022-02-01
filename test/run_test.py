#!/usr/bin/env python

"""Tests for unitig-caller"""

import subprocess
import os
import sys


sys.stderr.write("Testing Bifrost call\n")
subprocess.run("unitig-caller --call --refs refs.txt --out test_call --pyseer --rtab", shell=True, check=True)
sys.stderr.write("Testing rtab\n")
subprocess.run("python test-calls.py --method rtab --test test_call.rtab --expected bifrost_call.rtab", shell=True, check=True) # checks results match
sys.stderr.write("Testing pyseer\n")
subprocess.run("python test-calls.py --method pyseer --test test_call.pyseer --expected bifrost_call.pyseer", shell=True, check=True) # checks results match

sys.stderr.write("Testing Bifrost query\n")
subprocess.run("unitig-caller --query --refs refs.txt --unitigs test_unitigs.fasta --out test_query --pyseer --rtab", shell=True, check=True) # checks if can be run
subprocess.run("python test-calls.py --method rtab --test test_query.rtab --expected bifrost_query.rtab --strict", shell=True, check=True) # checks results match
subprocess.run("python test-calls.py --method pyseer --test test_query.pyseer --expected bifrost_query.pyseer --strict", shell=True, check=True) # checks results match

sys.stderr.write("Testing simple mode\n")
subprocess.run("unitig-caller --simple --refs refs.txt --unitigs test_unitigs_simple.fasta --out simple_calls", shell=True, check=True) # checks if can be run
subprocess.run("python test-calls.py --method pyseer --test simple_calls.pyseer --expected simple_results.pyseer --strict", shell=True, check=True) # checks results match

sys.stderr.write("Tests completed\n")
