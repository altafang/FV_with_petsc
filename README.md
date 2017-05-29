# FV_with_petsc

This C++ code numerically solves Poisson's equation using the finite volume method on a regular grid in 3D.
The solve can be done in parallel on multiple processors as specified at runtime. 

The domain is only decomposed along the 'z' direction, which is the direction of applied constant-value boundary conditions. 
Lateral boundary conditions may be either zero-flux or periodic.

Note: this package is very much under development. The goal is to eventually include a significantly broader set of features.

## Dependencies

- PETSc
    - Installation of PETSc can take a while and be tricky. Instructions are coming soon.
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

## Motivation and example applications

This package uses PETSc to quickly solve simple linear finite difference/finite volume PDEs. The idea is to wrap the specific PETSc calls so that the user does not need to worry about learning PETSc syntax and usage details.

For example, this package can be used to numerically calculate the effective conductivity of a material with an arbitrary microstructure.

