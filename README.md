# FV_with_petsc

This C++ code uses PETSc to numerically solves Poisson's equation using the finite volume method on a regular grid in 3D.
The solve can be done in parallel on multiple processors as specified at runtime. 
The idea is to wrap the specific PETSc calls so that the user does not need to worry about learning PETSc syntax and usage details.
For example, this package can be used to numerically calculate the effective conductivity of a material with an arbitrary microstructure.
PETSc is an open-source toolkit for scientific computation.

The domain is only decomposed along the 'z' direction, which is the direction of applied constant-value boundary conditions. 
Lateral boundary conditions may be either zero-flux or periodic.

The generalized Poisson equation solved here is:

```
div (sigma grad phi) = source
```

Note: this package is very much under development. The goal is to eventually include a significantly broader set of features.

## Dependencies

- a C++ compiler such as gcc
- PETSc
    - Installation of PETSc may take a while and may be tricky. It is essential though. See the Installing PETSc section below.
- hdf5
- MPI

## Usage

To compile, run `make`. The executable `solve_poisson` will be in the `bin/` directory. To run on 2 processors, for example, execute the following command:

```
mpiexec -np 2 ./solve_poisson
```

Inputs:
- input.txt -- text file that specifies NX, NY, NZ, DELTA_X, PHI_UPPER, PHI_LOWER, X_BC_TYPE, Y_BC_TYPE
where X_BC_TYPE and Y_BC_TYPE are either periodic or zeroflux 
- sigma.h5 -- hdf5 file that specifies the conductivity field
- source.h5 -- hdf5 that specifies the source field

Outputs:
- phi.h5 -- hdf5 file with the solved field

For example(s), see the `test/` directory.

## Installing PETSc

There are many installation configuration options, as described in the [PETSc installation documentation](http://www.mcs.anl.gov/petsc/documentation/installation.html).
Here I give rough instructions that worked for me, but it may.

1. First download the latest `.tar.gz` release from [here](http://www.mcs.anl.gov/petsc/download/index.html). 
2. Untar it by executing the command `tar -xvf <tarball filename here>`. 
3. Make or identify a directory where you want the installation to go, for example `/home/<username>/software/`.
4. Enter the untarred directory and execute `./configure ` with the options described below.
5. Add to your `~/.bash_profile` a line that says something like `export PETSC_DIR=/home/<username>/software/petsc-3.6.3/` where 3.6.3 is replaced with your version number.

If you are on a shared cluster that already has some software for scientific computing installed, try the following `./configure ` options:

```
./configure --prefix=/home/<username>/software/petsc-3.6.3 --with-debugging=0 COPTFLAGS=
-O3 CXXOPTFLAGS=-O3 --download-hypre --with-hdf5=1 --download-sprng
```

If you don't think you already have certain packages installed, or if the above ends up not working out, try adding the following flags (or a subset them) to those above:

```
--with-cc=gcc --with-cxx=g++ --with-fc=gfortran-mp-4.9 --download-fblaslapack --download-mpich --download-hdf5
```

Some key points to note are:
- Be sure to specify `--download-hypre` since this code relies on the HYPRE preconditioner.
- `--with-debugging=0 COPTFLAGS=-O3 CXXOPTFLAGS=-O3` turns off debugging and turns on optimization, which makes the code run faster.
- `--download-sprng` is optional, but it configures PETSc to enable use of the SPRNG parallel random number generator, which may be useful if you need random numbers (e.g. for some phase-field simulations).