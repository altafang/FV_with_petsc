# FV_with_petsc

This C++ code numerically solves Poisson's equation using the finite volume method on a regular grid in 3D.
Multiple processors in parallel may be used, as specified at runtime.

To compile, run `make` in the src/ directory. To run on 2 processors, for example, execute the following command:

```
mpiexec -np 2 ./main
```

Note: this package is very much under development. The goal is to eventually include a significantly broader set of features.

## Dependencies
- PETSc
    - Installation of PETSc can take a while and be tricky. Detailed instructions coming soon
- hdf5
- MPI

## Usage
Inputs:
- input.txt -- specifies NX, NY, NZ, DELTA_X, PHI_APPLIED
- sigma.h5 -- specifies the conductivity field
- source.h5 -- specifies the source field

Outputs:
- phi.h5 -- an hdf5 file with the solved field
