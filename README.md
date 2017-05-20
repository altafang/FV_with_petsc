# FV_with_petsc

Numerically solve Poisson's equation using the finite volume method on a regular grid in 3D with the option of using multiple processors in parallel.

This package is very much under development. The goal is to eventually include a significantly broader set of features.

## Dependencies
- PETSc
    - Installation of PETSc can take a while and be tricky. Detailed instructions coming soon
- hdf5
- MPI

## Usage
Inputs:
- input.txt
- an hdf5 file specifying the conductivity field

Outputs:
- an hdf5 file with the solved field
