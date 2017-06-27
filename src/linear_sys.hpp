#ifndef LINEAR_SYS_H
#define LINEAR_SYS_H

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>

// LinearSys is a struct since its members are public and it only has constructor and destructor, 
// no other member functions.
// LinearSys holds and sets up a linear system that can be solved with PETSc.
struct LinearSys
{
    Mat A;
    Vec b;
    KSP ksp;
    PC pc;
    
    // Constructor and destructor
    LinearSys(const int &total_N, const int &n_stencil_nonzero);
    ~LinearSys();
};

#endif
