// The code here is adapted from
// http://www.mcs.anl.gov/petsc/petsc-3.6/src/ksp/ksp/examples/tutorials/ex2.c.html

#include "linear_solver.hpp"

// constructor
LinearSolver::LinearSolver(DM *da, NonLocalField *phi, NonLocalField *sigma, Field *source):
    da(da), phi(phi), sigma(sigma), source(source)
{
    int nx, ny, nz;
    DMDAGetInfo(*da, NULL, &nx, &ny, &nz, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    
    int total_N = nx*ny*nz;
    MatCreate(PETSC_COMM_WORLD,&A);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,total_N,total_N);
    MatSetFromOptions(A);
    
    int n_stencil_nonzero = 7; // For 3D
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
LinearSolver::~LinearSolver()
{
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);
}

void LinearSolver::run_solver(double DELTA_X)
{
    int nx, ny, nz;
    DMBoundaryType x_BC_type, y_BC_type; // Lateral boundary condition types
    DMDAGetInfo(*da, NULL, &nx, &ny, &nz, NULL, NULL, NULL, NULL, NULL, &x_BC_type, &y_BC_type, NULL, NULL);

    int i,j,k,Ii,J,Istart,Iend;
    double v;
    
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
    
    double coeffpxhalf, coeffpyhalf, coeffpzhalf, coeffmxhalf, coeffmyhalf, coeffmzhalf;
    for (Ii=Istart; Ii<Iend; Ii++)
    {
        i = Ii/(nx*ny); j = (Ii - i*(nx*ny))/nx; k = Ii - j*nx - i*nx*ny;
        
        // Careful with signs: here, positive on the diagonals.
        // Retrieve coeff values, including neighboring values.
        coeffpxhalf = 0.5*(sigma->local_array[i][j][k] + sigma->local_array[i][j][k+1])/(DELTA_X*DELTA_X);
        coeffmxhalf = 0.5*(sigma->local_array[i][j][k] + sigma->local_array[i][j][k-1])/(DELTA_X*DELTA_X);
        coeffpyhalf = 0.5*(sigma->local_array[i][j][k] + sigma->local_array[i][j+1][k])/(DELTA_X*DELTA_X);
        coeffmyhalf = 0.5*(sigma->local_array[i][j][k] + sigma->local_array[i][j-1][k])/(DELTA_X*DELTA_X);
        coeffpzhalf = 0.5*(sigma->local_array[i][j][k] + sigma->local_array[i+1][j][k])/(DELTA_X*DELTA_X);
        coeffmzhalf = 0.5*(sigma->local_array[i][j][k] + sigma->local_array[i-1][j][k])/(DELTA_X*DELTA_X);
        
        // Fill in the finite-volume stencil
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
            v = -coeffmyhalf;
            MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        if (j<ny-1)
        {
            J = Ii + nx;
            v = -coeffpyhalf;
            MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        if (k>0)
        {
            J = Ii - 1;
            v = -coeffmxhalf;
            MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        if (k<nx-1)
        {
            J = Ii + 1;
            v = -coeffpxhalf;
            MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        
        // Periodic boundary conditions in y
        if (y_BC_type == DM_BOUNDARY_PERIODIC)
        {
            if (j == 0)
            {
                J = Ii + nx*ny - nx;
                v = -coeffmyhalf;
                MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
            if (j == ny - 1)
            {
                J = Ii - nx*ny + nx;
                v = -coeffpyhalf;
                MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
        }
        // Zero-flux boundary conditions in y. GHOSTED instead of MIRROR in 3D
        else if (y_BC_type == DM_BOUNDARY_GHOSTED)
        {
            if (j == 0)
            {
                J = Ii + nx;
                v = -coeffmyhalf-coeffpyhalf;
                MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
            if (j == ny - 1)
            {
                J = Ii - nx;
                v = -coeffmyhalf-coeffpyhalf;
                MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
        }
        
        // Periodic boundary conditions in x
        if (x_BC_type == DM_BOUNDARY_PERIODIC)
        {
            if (k == 0)
            {
                J = Ii + nx - 1;
                v = -coeffmxhalf;
                MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
            if (k == nx - 1)
            {
                J = Ii - nx + 1;
                v = -coeffpxhalf;
                MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
        }
        // Zero-flux boundary conditions in x. GHOSTED instead of MIRROR in 3D
        else if (x_BC_type == DM_BOUNDARY_GHOSTED)
        {
            if (k == 0)
            {
                J = Ii + 1;
                v = -coeffmxhalf-coeffpxhalf;
                MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
            if (k == nx - 1)
            {
                J = Ii - 1;
                v = -coeffmxhalf-coeffpxhalf;
                MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
        }
        
        // Derivative boundary conditions in z
        if (phi->bc->upper_BC_type == derivativeBC)
        {
            if (i == 0)
            {
                J = Ii + nx*ny;
                v = -coeffmzhalf-coeffpzhalf;
                MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
        }
        // Derivative boundary conditions in z
        if (phi->bc->lower_BC_type == derivativeBC)
        {
            if (i == nz-1)
            {
                J = Ii - nx*ny;
                v = -coeffmzhalf-coeffpzhalf;
                MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
        }
        
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
            coeffmzhalf = 0.5*(sigma->local_array[i][j][k] + sigma->local_array[i-1][j][k])/(DELTA_X*DELTA_X);
            v = -source->global_array[i][j][k];
            
            if (phi->bc->upper_BC_type == constantBC)
            {
                v += coeffmzhalf*phi->bc->upper_BC_val;
            }
            else
            {
                v += 2.*DELTA_X*coeffmzhalf*phi->bc->upper_BC_val;
            }
        }
        // Bottom side
        else if (i == nz - 1)
        {
            coeffpzhalf = 0.5*(sigma->local_array[i][j][k] + sigma->local_array[i+1][j][k]);
            v = -source->global_array[i][j][k];
            
            if (phi->bc->lower_BC_type == constantBC)
            {
                v += coeffpzhalf*phi->bc->lower_BC_val;
            }
            else
            {
                v += 2.*DELTA_X*coeffpzhalf*phi->bc->lower_BC_val;
            }
        }
        // Center points
        else
        {
            v = -source->global_array[i][j][k];
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