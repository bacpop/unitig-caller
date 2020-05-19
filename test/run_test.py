#!/usr/bin/env python

"""Tests for unitig-caller"""

import subprocess
import os
import sys
import shutil
import argparse


# python tests
# create sketches
sys.stderr.write("Testing bifrost build\n")
subprocess.run("unitig-caller --build --refs refs.txt --output test_build", shell=True, check=True)
# calculate distances
sys.stderr.write("Testing bifrost query\n")
subprocess.run("unitig-caller --query --graph-prefix test_bifrost --unitigs query_unitigs.fasta --output test_query", shell=True, check=True) # checks if can be run
subprocess.run("python test-intersection.py", shell=True, check=True) # checks results match
# Joining
sys.stderr.write("Testing simple mode\n")
subprocess.run("unitig-caller ---simple --refs simple_refs.txt --unitigs queries.txt --output calls", shell=True, check=True) # checks if can be run
subprocess.run("python test-intersection.py", shell=True, check=True) # checks results match

sys.stderr.write("Tests completed\n")

