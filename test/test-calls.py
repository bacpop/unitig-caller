#!/usr/bin/env python

import os, sys
import argparse

description = 'Test unitig calls match'
parser = argparse.ArgumentParser(description=description,
                                    prog='test-calls')

parser.add_argument('--method', help="File type", required=True)
parser.add_argument('--test', help="Calls produced", required=True)
parser.add_argument('--expected', help="Calls expected", required=True)

options = parser.parse_args()

if (options.method == "rtab"):
    with open(options.expected, 'r') as expected_rtab, open(options.test, 'r') as test_rtab:
        expected_samples = expected_rtab.readline().rstrip().split("\t")[1:]
        test_samples = test_rtab.readline().rstrip().split("\t")[1:]
        if (len(expected_samples) != len(test_samples)):
            sys.stderr.write("Length of rtab header mismatches\n")
            raise RuntimeError("Length of rtab header mismatches")
        elif len(set(expected_samples).intersection(test_samples)) != len(expected_samples):
            sys.stderr.write("Different samples contained\n")
            raise RuntimeError("Different samples contained")

        for expected_line, test_line in zip(expected_rtab, test_rtab):
            expected_calls = expected_line.rstrip().split("\t")[1:]  
            test_calls = test_line.rstrip().split("\t")[1:]
            for sample_idx, call in enumerate(expected_calls):
                if call != test_calls[test_samples.index(expected_samples[sample_idx])]:
                    sys.stderr.write("Calls mismatch for " + expected_samples[sample_idx] + "\n")
                    raise RuntimeError("Calls mismatch for " + expected_samples[sample_idx])

elif (options.method == "pyseer"):
    d = {}
    with open(options.expected, 'r') as expected_pyseer:
        for call in expected_pyseer:
            call_fields = call.rstrip().split(" ")
            d[call_fields[0]] = set()
            for sample in call_fields[2:]:
                d[call_fields[0]].add(sample)

    with open(options.test, 'r') as test_pyseer:
        for call in test_pyseer:
            call_fields = call.rstrip().split(" ")
            if not call_fields[0] in d.keys():
                sys.stderr.write("No call for " + call_fields[0] + "\n")
                raise RuntimeError("No call for " + call_fields[0])
            test_calls = set()
            for sample in call_fields[2:]:
                test_calls.add(sample)

            if test_calls != d[call_fields[0]]:
                sys.stderr.write("Calls mismatch for " + call_fields[0] + "\n")
                raise RuntimeError("Calls mismatch for " + call_fields[0])

else:
    raise RuntimeError("Unknown method type")

