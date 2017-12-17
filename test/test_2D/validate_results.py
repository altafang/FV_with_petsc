# Example script for reading results
import numpy
import h5py
from matplotlib import pyplot
import subprocess
import numpy.testing

# XXX Assumes that all testing scripts are run from within their respective directories!
import sys
sys.path.insert(0, '../')
from tools import *

def validate_results(path="", plot=True, variation='X'):

    if not ((variation == 'X') or (variation == 'Y')):
        print "Invalid dimension of variation, must be X or Y"
        return

    # Read parameters from the input file
    param_dict = extract_params(path)

    NX = int(param_dict['NX'])
    NY = int(param_dict['NY'])

    dims = [NX, NY]

    phi = read_hdf5(path + "phi.h5", dims)

    # Check that it is linear
    if variation == 'X':
        theoretical_1D = numpy.linspace(1, 0, NX+2)[1:-1]
        theoretical_2D = numpy.broadcast_to(theoretical_1D.reshape(1, NX), (NY, NX))
    else:
        theoretical_1D = numpy.linspace(1, 0, NY+2)[1:-1]
        theoretical_2D = numpy.broadcast_to(theoretical_1D.reshape(NY, 1), (NY, NX))
            
    # Assert if error is too large
    try:
        rtol = 1.e-3
        numpy.testing.assert_allclose(phi, theoretical_2D, rtol=rtol)
        # If no exception raised so far, then it passed
        print "Ok: solution is within a factor of %f of theoretical." % rtol
    except Exception as e:
        print e

    if plot:
        pyplot.imshow(phi)
        pyplot.savefig("phi.png")
        pyplot.show()
        
        # Check that it is linear
        if variation == 'X':
            pyplot.plot(range(NX), theoretical_1D, 'r')
            pyplot.plot(range(NX), phi[NY/2,:], 'k-o')
        else:
            pyplot.plot(range(NY), theoretical_1D, 'r')
            pyplot.plot(range(NY), phi[:, NX/2], 'k-o')        
        pyplot.savefig("slice1D_middle.png")
        pyplot.show()

if __name__ == '__main__':
    validate_results(plot=True, variation='Y')

