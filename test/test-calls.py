#!/usr/bin/env python

import os, sys
import argparse

description = 'Test unitig calls match'
parser = argparse.ArgumentParser(description=description,
                                    prog='test-calls')

parser.add_argument('--method', help="File type", required=True)
parser.add_argument('--test', help="Calls produced", required=True)
parser.add_argument('--expected', help="Calls expected", required=True)
parser.add_argument('--strict', action="store_true", default=False, help="Ensure unitigs exactly match")

options = parser.parse_args()

if (options.method == "rtab"):
    with open(options.expected, 'r') as expected_rtab, open(options.test, 'r') as test_rtab:
        expected_samples = expected_rtab.readline().rstrip().split("\t")[1:]
        test_samples = test_rtab.readline().rstrip().split("\t")[1:]
        if (len(expected_samples) != len(test_samples)):
            sys.stderr.write("Length of rtab header mismatches\n")
            raise RuntimeError("Length of rtab header mismatches")
        elif set(expected_samples) != set(test_samples):
            sys.stderr.write("Different samples contained\n")
            raise RuntimeError("Different samples contained")

        expected_dict = {}
        test_dict = {}
        for expected_line, test_line in zip(expected_rtab, test_rtab):
            expected_calls = expected_line.rstrip().split("\t")
            test_calls = test_line.rstrip().split("\t")
            expected_dict[expected_calls[0]] = expected_calls[1:]
            test_dict[test_calls[0]] = test_calls[1:]
        if (len(expected_dict) != len(test_dict)):
            raise RuntimeError(f"Number of unitigs differs expected {len(expected_dict)} found {len(test_dict)}")
        for unitig, calls in expected_dict.items():
            if unitig in test_dict:
                test_calls = expected_dict[unitig]
                for sample_idx, call in enumerate(calls):
                    if call != test_calls[test_samples.index(expected_samples[sample_idx])]:
                        sys.stderr.write("Calls mismatch for " + expected_samples[sample_idx] + "\n")
                        raise RuntimeError("Calls mismatch for " + expected_samples[sample_idx])
            else:
                sys.stderr.write(f"No call for {unitig}\n")
                if options.strict:
                    raise RuntimeError("Unitigs mismatch")

elif (options.method == "pyseer"):
    d = {}
    with open(options.expected, 'r') as expected_pyseer:
        for call in expected_pyseer:
            call_fields = call.rstrip().split(" ")
            d[call_fields[0]] = set()
            for sample in call_fields[2:]:
                d[call_fields[0]].add(sample)

    num_calls = 0
    with open(options.test, 'r') as test_pyseer:
        for call in test_pyseer:
            num_calls += 1
            call_fields = call.rstrip().split(" ")
            if not call_fields[0] in d.keys():
                sys.stderr.write("No call for " + call_fields[0] + "\n")
                if options.strict:
                    raise RuntimeError("Unitigs mismatch")
            else:
                test_calls = set()
                for sample in call_fields[2:]:
                    test_calls.add(sample)

                if test_calls != d[call_fields[0]]:
                    sys.stderr.write("Calls mismatch for " + call_fields[0] + "\n")
                    raise RuntimeError("Calls mismatch for " + call_fields[0])

    if (len(d) != num_calls):
        raise RuntimeError(f"Number of unitigs differs expected {len(d)} found {num_calls}")

else:
    raise RuntimeError("Unknown method type")

