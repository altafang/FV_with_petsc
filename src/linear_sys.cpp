// The code here is adapted from
// http://www.mcs.anl.gov/petsc/petsc-3.6/src/ksp/ksp/examples/tutorials/ex2.c.html

#include "linear_sys.hpp"

// constructor
LinearSys::LinearSys(const int &total_N, const int &n_stencil_nonzero)
{
    MatCreate(PETSC_COMM_WORLD,&A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,total_N,total_N);
    MatSetFromOptions(A);
    
    MatMPIAIJSetPreallocation(A,n_stencil_nonzero,NULL,n_stencil_nonzero,NULL);
    MatSeqAIJSetPreallocation(A,n_stencil_nonzero,NULL);
    
    // Create parallel vectors.
    // - We form 1 vector from scratch and then duplicate as needed.
    // - When u sing VecCreate(), VecSetSizes and VecSetFromOptions()
    // in this example, we specify only the
    // vector's global dimension; the parallel partitioning is determined
    // at runtime.
    // - When solving a linear system, the vectors and matrices MUST
    // be partitioned accordingly.  PETSc automatically generates
    // appropriately partitioned matrices and vectors when MatCreate()
    // and VecCreate() are used with the same communicator.
    // - The user can alternatively specify the local vector and matrix
    // dimensions when more sophisticated partitioning is needed
    // (replacing the PETSC_DECIDE argument in the VecSetSizes() statement
    // below).
    
    VecCreate(PETSC_COMM_WORLD,&b);
    VecSetSizes(b,PETSC_DECIDE,total_N);
    VecSetFromOptions(b);
    
    // Create linear solver context
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    
    KSPGetPC(ksp,&pc);
    // Specify preconditioner to use
    PCSetType(pc,PCHYPRE);
    KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    // Here, choose KSPGCR, a solver that works well under many conditions
    KSPSetType(ksp, KSPGCR);
    
    // optional: tell solvers not to zero out x every time.
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    PCSetInitialGuessNonzero(pc, PETSC_TRUE);
    
    // Allow options to be set from the command line.
    KSPSetFromOptions(ksp);
    PCSetFromOptions(pc);
}

// destructor
LinearSys::~LinearSys()
{
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);
}
