# Copyright 2019-2020 John Lees, Sam Horsfield


'''Functions to run Bifrost'''

import os, sys
import subprocess
import re
import tempfile

def rtab_format(unitig_map, input_colour_pref, outfile):
    """Creates RTAB file.

    Args:
        unitig_map (dictionary)
            dictionary mapping each unitig to a binary colour vector.
        input_colour_pref (list)
            list of input colour file prefixs.
        outfile (str)
            output prefix for pyseer file
    """
    with open(outfile, "w") as o:
        o.write("Unitig_sequence")
        for index, colour_name in enumerate(input_colour_pref):
            o.write("\t" + colour_name)
            if index == len(input_colour_pref) - 1:
                o.write("\n")
        for unitig, colours_vector in unitig_map.items():
            o.write(unitig)
            for index, colour in enumerate(colours_vector):
                colour = int(colour == True)
                o.write("\t" + str(colour))
                if index == len(colours_vector) - 1:
                    o.write("\n")

def pyseer_format(unitig_map, input_colour_pref, outfile):
    """Creates pyseer file.

    Args:
        unitig_map (dictionary)
            dictionary mapping each unitig to a binary colour vector.
        input_colour_pref (list)
            list of input colour file prefixs.
        outfile (str)
            output prefix for pyseer file
    """
    with open(outfile, "w") as o:
        for unitig, colours_vector in unitig_map.items():
            o.write(unitig + " " + "|")
            for index, colour in enumerate(colours_vector):
                if colour:
                    o.write(" " + input_colour_pref[index] + ":1")
                if index == len(colours_vector) - 1:
                    o.write("\n")
