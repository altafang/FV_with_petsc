# Example script for reading results
import numpy
import h5py
import pylab
import matplotlib
from matplotlib import pyplot
import subprocess
from tools import *

def plot_results(path = ""):

    # Read parameters from the input file
    param_dict = extract_params(path)

    NX = int(param_dict['NX'])
    NY = int(param_dict['NY'])
    NZ = int(param_dict['NZ'])

    dims = [NX, NY, NZ]

    phi = read_hdf5(path + "phi.h5", dims)

    pyplot.imshow(phi[NZ/2,:,:])
    pyplot.savefig("halfway_down_slice.png")
    pyplot.show()


if __name__ == '__main__':
    plot_results()

