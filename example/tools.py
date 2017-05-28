# Utility functions for manipulating input/output

#XXX Should move this into a separate analysis tools area

import h5py
import numpy
from matplotlib import pyplot

def extract_params(path=""):
    """Read parameters from input.txt into a dictionary of strings"""
    param_file = open(path + "input.txt", "r")
    param_dict = {}
    for line in param_file:
        if not line.isspace() and line[0] != '#':
            line_list = line.split("=")
            param_dict[line_list[0].strip()] = line_list[1].strip()
    param_file.close()
    return param_dict

def read_hdf5(filename, dims):
    f = h5py.File(filename, 'r')
    array = numpy.array(f[f.keys()[0]])
    f.close()
    return array.reshape(dims[::-1]) # note reversed order here

