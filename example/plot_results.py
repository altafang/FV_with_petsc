# Example script for reading results
import numpy
import h5py
import pylab
import matplotlib
from matplotlib import pyplot
import subprocess

###################

def read_hdf5(filename, dims):
    f = h5py.File(filename, 'r')
    array = numpy.array(f[f.keys()[0]])
    f.close()
    return array.reshape(dims[::-1]) # note reversed order here

def extract_params(path=""):
    param_file = open(path + "input.txt", "r")
    param_dict = {}
    for line in param_file:
        if not line.isspace() and line[0] != '#':
            line_list = line.split("=")
            param_dict[line_list[0].strip()] = line_list[1].strip()
    param_file.close()
    return param_dict

###################

path = ""

# Read parameters from the input file
param_dict = extract_params(path)

NX = int(param_dict['NX'])
NY = int(param_dict['NY'])
NZ = int(param_dict['NZ'])

dims = [NX, NY, NZ]

###################

phi = read_hdf5(path + "phi_%03d.h5" % 0, dims)

pyplot.imshow(phi[NZ/2,:,:])
pyplot.savefig("halfway_down_slice.png")
pyplot.show()




