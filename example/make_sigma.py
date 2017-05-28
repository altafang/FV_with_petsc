# Example script for setting up sigma field
import h5py
import numpy
import numpy.random
from matplotlib import pyplot

###############

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

##############

path = ""
plot = True

# Read parameters from the input file
param_dict = extract_params(path)

NX = int(param_dict['NX'])
NY = int(param_dict['NY'])
NZ = int(param_dict['NZ'])

# Parameters for the features & values of sigma
LEN = 30
SPACING = 20
CUTOFF_VAL = 1.e-3

#############

file = path + "sigma_000.h5"
f = h5py.File(file, 'w')

# Start with everything as ones
sigma = numpy.ones((NZ, NY, NX))

if LEN > 0:
    # Put in a rectangular feature
    centery = NY/2
    centerx = NX/2
    
    sigma[SPACING:-SPACING,centery-LEN/2:centery+LEN/2,centerx-LEN/2:centerx+LEN/2] = CUTOFF_VAL

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
