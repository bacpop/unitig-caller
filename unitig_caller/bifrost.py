# Copyright 2019-2020 John Lees, Sam Horsfield


'''Functions to run Bifrost'''

import os, sys
import subprocess
import re
import tempfile

def rtab_format(csv_tmp, input_colour_pref, outfile):
    """Creates RTAB file.

    Args:
        csv_tmp (filename)
            csv mapping each unitig to a binary colour vector.
        input_colour_pref (list)
            list of input colour file prefixs.
        outfile (str)
            output prefix for pyseer file
    """
    with open(outfile, "w") as o, open(csv_tmp, "r") as i:
        o.write("Unitig_sequence")
        for colour_name in input_colour_pref:
            o.write("\t" + colour_name)
        o.write("\n")
        for line in i:
            split_line = line.rstrip().split(",")
            unitig = split_line[0]
            colours_vector = split_line[1:]
            o.write(unitig)
            for colour in colours_vector:
                colour = int(colour == True)
                o.write("\t" + str(colour))
            o.write("\n")

def pyseer_format(csv_tmp, input_colour_pref, outfile):
    """Creates pyseer file.

    Args:
        csv_tmp (filename)
            csv mapping each unitig to a binary colour vector.
        input_colour_pref (list)
            list of input colour file prefixs.
        outfile (str)
            output prefix for pyseer file
    """
    with open(outfile, "w") as o, open(csv_tmp, "r") as i:
        for line in i:
            split_line = line.rstrip().split(",")
            unitig = split_line[0]
            colours_vector = split_line[1:]
            o.write(unitig + " " + "|")
            for index, colour in enumerate(colours_vector):
                if colour:
                    o.write(" " + input_colour_pref[index] + ":1")
            o.write("\n")
