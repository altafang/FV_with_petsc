# Example script for setting up sigma field
import h5py
import numpy
from matplotlib import pyplot
from tools import *

def make_sigma(path="", plot=True):

    # Read parameters from the input file
    param_dict = extract_params(path)

    NX = int(param_dict['NX'])
    NY = int(param_dict['NY'])
    NZ = int(param_dict['NZ'])

    file = path + "sigma.h5"
    f = h5py.File(file, 'w')

    # Constant
    sigma = numpy.ones((NZ, NY, NX))

    # Write to file
    dset = f.create_dataset("sigma", data=sigma)
    f.close()

    # Debugging: plot
    if plot:
        pyplot.imshow(sigma[:,:, NX/2])
        pyplot.savefig("sigma.png")
        pyplot.show()
        
        pyplot.imshow(sigma[NZ/2,:,:])
        pyplot.savefig("sigma_slicez.png")
        pyplot.show()
        
        pyplot.plot(range(NZ), sigma[:,NY/2,NX/2], 'ko')
        pyplot.savefig("sigma_slice.png")
        pyplot.show()

def make_source(path=""):
    
    # Read parameters from the input file
    param_dict = extract_params(path)
    
    NX = int(param_dict['NX'])
    NY = int(param_dict['NY'])
    NZ = int(param_dict['NZ'])
    
    file = path + "source.h5"
    f = h5py.File(file, 'w')
    
    source = numpy.zeros((NZ, NY, NX))
    
    # Write to file
    dset = f.create_dataset("source", data=source)
    f.close()

if __name__ == '__main__':
    make_sigma(plot=False)
    make_source()
