# Example script for reading results
import numpy
import h5py
import pylab
import matplotlib
from matplotlib import pyplot
import subprocess
import numpy.testing
from tools import *

def validate_results(path = "", plot=True):

    # Read parameters from the input file
    param_dict = extract_params(path)

    NX = int(param_dict['NX'])
    NY = int(param_dict['NY'])
    NZ = int(param_dict['NZ'])

    dims = [NX, NY, NZ]

    phi = read_hdf5(path + "phi.h5", dims)

    # Check that it is linear
    theoretical_1D = numpy.linspace(1, 0, NZ+2)[1:-1]
    theoretical_3D = numpy.broadcast_to(theoretical_1D.reshape(NZ,1,1), (NZ, NY, NX))
    # Assert if error is too large
    numpy.testing.assert_allclose(phi, theoretical_3D, rtol=1.e-3)

    if plot:
        pyplot.imshow(phi[:,:,NX/2])
        pyplot.savefig("slice2D_NXhalf.png")
        pyplot.show()
        
        # Check that it is linear
        pyplot.plot(range(NZ), theoretical_1D, 'r')
        pyplot.plot(range(NZ), phi[:,NY/2,NX/2], 'k-o')
        pyplot.savefig("slice1D_middle.png")
        pyplot.show()

if __name__ == '__main__':
    validate_results()

