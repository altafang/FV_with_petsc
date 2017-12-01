# Example script for setting up sigma field
import h5py
import numpy
from matplotlib import pyplot

# XXX Assumes that all testing scripts are run from within their respective directories!
import sys
sys.path.insert(0, '../')
from tools import *

def make_sigma(path="", plot=True):

    # Read parameters from the input file
    param_dict = extract_params(path)

    NX = int(param_dict['NX'])
    NY = int(param_dict['NY'])

    file = path + "sigma.h5"
    f = h5py.File(file, 'w')

    # Constant
    sigma = numpy.ones((NY, NX))

    # Write to file
    dset = f.create_dataset("sigma", data=sigma)
    f.close()

    # Debugging: plot
    if plot:
        pyplot.imshow(sigma)
        pyplot.savefig("sigma.png")
        pyplot.show()

def make_source(path=""):
    
    # Read parameters from the input file
    param_dict = extract_params(path)
    
    NX = int(param_dict['NX'])
    NY = int(param_dict['NY'])
    
    file = path + "source.h5"
    f = h5py.File(file, 'w')
    
    source = numpy.zeros((NY, NX))
    
    # Write to file
    dset = f.create_dataset("source", data=source)
    f.close()

if __name__ == '__main__':
    make_sigma(plot=False)
    make_source()
