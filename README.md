# FV_with_petsc

This C++ code numerically solves Poisson's equation using the finite volume method on a regular grid in 3D.
Multiple processors in parallel may be used, as specified at runtime.

Note: this package is very much under development. The goal is to eventually include a significantly broader set of features.

## Dependencies

- PETSc
    - Installation of PETSc can take a while and be tricky. Detailed instructions coming soon
- hdf5
- MPI

## Usage

To compile, run `make` in the src/ directory. To run on 2 processors, for example, execute the following command:

```
mpiexec -np 2 ./main
```

Inputs:
- input.txt -- text file that specifies NX, NY, NZ, DELTA_X, PHI_APPLIED
- sigma.h5 -- hdf5 file that specifies the conductivity field
- source.h5 -- hdf5 that specifies the source field

Outputs:
- phi.h5 -- hdf5 file with the solved field

## Motivation/Applications

This package is a wrapper that uses PETSc to quickly solve simple linear finite difference/finite volume PDEs. 
Thus, the user does not need to worry about PETSc syntax and implementation details.
Also, large systems can be solved efficiently and rapidly because computation can be done in parallel over many processors, and the code is in C++.
For example, this package can be used to numerically calculate the effective conductivity of a material with an arbitrary microstructure.

