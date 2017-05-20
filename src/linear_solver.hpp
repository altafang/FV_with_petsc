// The code here is adapted from
// http://www.mcs.anl.gov/petsc/petsc-3.6/src/ksp/ksp/examples/tutorials/ex2.c.html

#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include "field.hpp"
#include "nonlocal_field.hpp"

class LinearSolver
{
    public:
        LinearSolver(DM *da, NonLocalField *phi, NonLocalField *sigma);
        ~LinearSolver();
        void run_solver();
    
    private:
        DM *da; // Pointer to da. Also holds nx, ny, nz
        NonLocalField *phi; // The field we are solving for
        NonLocalField *sigma;
    
        Mat A;
        KSP ksp;
        PC pc;
        Vec b;
        int nx, ny, nz;
};

// constructor
LinearSolver::LinearSolver(DM *da, NonLocalField *phi, NonLocalField *sigma):
    da(da), phi(phi), sigma(sigma)
{
    DMDAGetInfo(*da, NULL, &nx, &ny, &nz, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    
    MatCreate(PETSC_COMM_WORLD,&A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,nx*ny*nz,nx*ny*nz);
    MatSetFromOptions(A);
    
    int n_stencil_nonzero = 7; // For 3D
    MatMPIAIJSetPreallocation(A,n_stencil_nonzero,NULL,n_stencil_nonzero,NULL);
    MatSeqAIJSetPreallocation(A,n_stencil_nonzero,NULL);
    
    // Create parallel vectors.
    // - We form 1 vector from scratch and then duplicate as needed.
    // - When using VecCreate(), VecSetSizes and VecSetFromOptions()
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
    VecSetSizes(b,PETSC_DECIDE,nx*ny*nz);
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
LinearSolver::~LinearSolver()
{
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);
}

void LinearSolver::run_solver()
{
    int i,j,k,Ii,J,Istart,Iend;
    double v;
    
    // XXX FIXME: types of BCs are currently hard-coded in
    
    // Currently, all PETSc parallel matrix formats are partitioned by
    // contiguous chunks of rows across the processors.  Determine which
    // rows of the matrix are locally owned.
    
    MatGetOwnershipRange(A,&Istart,&Iend);
    
    // Set matrix elements for the 3-D, seven-point stencil in parallel.
    // - Each processor needs to insert only elements that it owns
    // locally (but any non-local elements will be sent to the
    // appropriate processor during matrix assembly).
    // - Always specify global rows and columns of matrix entries.
    // 
    // Note: this uses the less common natural ordering that orders first
    // all the unknowns for x = h then for x = 2h etc; Hence you see J = Ii +- n
    // instead of J = I +- m as you might expect. The more standard ordering
    // would first do all variables for y = h, then y = 2h etc.
    
    double thiscoeff, coeffpxhalf, coeffpyhalf, coeffpzhalf, coeffmxhalf, coeffmyhalf, coeffmzhalf;
    for (Ii=Istart; Ii<Iend; Ii++)
    {
        i = Ii/(nx*ny); j = (Ii - i*(nx*ny))/nx; k = Ii - j*nx - i*nx*ny;
        
        // Careful with signs: here, positive on the diagonals.
        // Retrieve coeff values, including neighboring values.
        thiscoeff = sigma->local_array[i][j][k];
        coeffpxhalf = 0.5*(thiscoeff + sigma->local_array[i][j][k+1]);
        coeffmxhalf = 0.5*(thiscoeff + sigma->local_array[i][j][k-1]);
        coeffpyhalf = 0.5*(thiscoeff + sigma->local_array[i][j+1][k]);
        coeffmyhalf = 0.5*(thiscoeff + sigma->local_array[i][j-1][k]);
        coeffpzhalf = 0.5*(thiscoeff + sigma->local_array[i+1][j][k]);
        coeffmzhalf = 0.5*(thiscoeff + sigma->local_array[i-1][j][k]);
        
        // Note: assumes certain boundary conditions.
        if (i>0)
        {
            J = Ii - nx*ny;
            v = -coeffmzhalf;
            MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        if (i<nz-1)
        {
            J = Ii + nx*ny;
            v = -coeffpzhalf;
            MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        if (j>0)
        {
            J = Ii - nx;
            // zero-flux lateral boundary condition
            if (j == ny-1)
            {
                v = -coeffmyhalf-coeffpyhalf;
            }
            else
            {
                v = -coeffmyhalf;
            }
            v = -coeffmyhalf;
            MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        if (j<ny-1)
        {
            J = Ii + nx;
            // zero-flux lateral boundary condition
            if (j == 0)
            {
                v = -coeffmyhalf-coeffpyhalf;
            }
            else
            {
                v = -coeffpyhalf;
            }
            v = -coeffpyhalf;
            MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        if (k>0)
        {
            J = Ii - 1;
            // zero-flux lateral boundary condition
            if (k == nx-1)
            {
                v = -coeffmxhalf-coeffpxhalf;
            }
            else
            {
                v = -coeffmxhalf;
            }
            v = -coeffmxhalf;
            MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        if (k<nx-1)
        {
            J = Ii + 1;
            // zero-flux lateral boundary condition
            if (k == 0)
            {
                v = -coeffmxhalf-coeffpxhalf;
            }
            else
            {
                v = -coeffpxhalf;
            }
            v = -coeffpxhalf;
            MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        
//        // Periodic boundary conditions in x and y
//        if (j == 0)
//        {
//            J = Ii + nx*ny - nx;
//            v = -coeffmyhalf;
//            MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
//        }
//        if (j == ny - 1)
//        {
//            J = Ii - nx*ny + nx;
//            v = -coeffpyhalf;
//            MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
//        }
//        if (k == 0)
//        {
//            J = Ii + nx - 1;
//            v = -coeffmxhalf;
//            MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
//        }
//        if (k == nx - 1)
//        {
//            J = Ii - nx + 1;
//            v = -coeffpxhalf;
//            MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
//        }
        
        // Diagonal term
        v = coeffmxhalf + coeffmyhalf + coeffmzhalf + coeffpxhalf + coeffpyhalf + coeffpzhalf;
        
        MatSetValues(A,1,&Ii,1,&Ii,&v,INSERT_VALUES);
    }
    
    // Assemble matrix, using the 2-step process:
    // MatAssemblyBegin(), MatAssemblyEnd()
    // Computations can be done while messages are in transition
    // by placing code between these two statements.
    
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    
    // Set right-hand-side vector.
    VecGetOwnershipRange(b,&Istart,&Iend);
    for (Ii=Istart; Ii<Iend; Ii++)
    {
        i = Ii/(nx*ny); j = (Ii - i*(nx*ny))/nx; k = Ii - j*nx - i*nx*ny;
        
        // Top side
        if (i == 0)
        {
            thiscoeff = sigma->local_array[i][j][k];
            coeffmzhalf = 0.5*(thiscoeff + sigma->local_array[i-1][j][k]);
            v = 0. + coeffmzhalf*phi->bc->upper_BC_val;
        }
        // Bottom side
        else if (i == nz - 1)
        {
            thiscoeff = sigma->local_array[i][j][k];
            coeffpzhalf = 0.5*(thiscoeff + sigma->local_array[i+1][j][k]);
            v = 0. + coeffpzhalf*phi->bc->lower_BC_val;
        }
        // Center points
        else
        {
            v = 0.;
        }
        VecSetValue(b, Ii, v, INSERT_VALUES);
    }
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    
    // Actually perform the solve
    KSPSetOperators(ksp,A,A);
    KSPSolve(ksp,b,phi->global_vec);
    
    return;
}

#endif
