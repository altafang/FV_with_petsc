#include "poisson_solver.hpp"
#include <unordered_map>

// Helper function: convert user-input lateral BC type strings to appropriate actual BC options 
DMBoundaryType& get_BC_type(std::string& type_name, bool warn=false)
{
    std::unordered_map<std::string, DMBoundaryType> lateral_BCs = {{"zeroflux", DM_BOUNDARY_GHOSTED}, {"periodic", DM_BOUNDARY_PERIODIC}};

    if (lateral_BCs.count(type_name) > 0)
    {
        return lateral_BCs[type_name];
    }
    else
    {
        if (warn)
        {
            PetscPrintf(PETSC_COMM_WORLD, "Warning: %s is not a valid lateral boundary condition type, using periodic instead\n", type_name.c_str());
        }
        return lateral_BCs["periodic"];
    }
}

void PoissonSolver::read_input(std::string input_file)
{
    std::map<std::string, std::string> params;
    read_parameters(params, input_file);
    
    unpack(params, "NX", NX);
    unpack(params, "NY", NY);
    unpack(params, "NZ", NZ);
    
    unpack(params, "DELTA_X", DELTA_X);
    
    unpack(params, "PHI_UPPER", PHI_UPPER);
    unpack(params, "PHI_LOWER", PHI_LOWER);
    
    unpack(params, "X_BC_TYPE", X_BC_TYPE);
    unpack(params, "Y_BC_TYPE", Y_BC_TYPE); 
}

// Constructor
PoissonSolver::PoissonSolver(std::string input_file)
{
    // Read input parameters from text file called "input.txt"
    // and unpack the input parameters.
    read_input(input_file);
    
    // Set up the distributed array used for all fields.
    // Lateral boundary conditions may be either periodic or zero-flux
    // For periodic, use DM_BOUNDARY_PERIODIC
    // For zero-flux, typically one should use DM_BOUNDARY_MIRROR
    // But, in 3D DM_BOUNDARY_MIRROR is not yet implemented so use DM_BOUNDARY_GHOSTED and manually fill the ghost cells    
    DMDACreate3d(PETSC_COMM_WORLD, get_BC_type(X_BC_TYPE, true), get_BC_type(Y_BC_TYPE, true), DM_BOUNDARY_GHOSTED, DMDA_STENCIL_STAR, \
                 NX, NY, NZ, 1, 1, PETSC_DECIDE, 1, 1, PETSC_NULL, PETSC_NULL, PETSC_NULL, &da);
    
    // Initialize with z boundary conditions
    phi = new NonLocalField(&da, constantBC, PHI_LOWER, constantBC, PHI_UPPER);
    PetscObjectSetName((PetscObject)phi->global_vec, "phi");
    
    // Zero derivative boundary conditions in z
    sigma = new NonLocalField(&da, derivativeBC, 0., derivativeBC, 0.);
    PetscObjectSetName((PetscObject)sigma->global_vec, "sigma");
    
    // Read sigma in from hdf5 file
    sigma->read_from_file("sigma.h5");
    sigma->send_global_to_local();
    
    source = new Field(&da);
    PetscObjectSetName((PetscObject)source->global_vec, "source");

    // Read source term in from hdf5 file
    source->read_from_file("source.h5");

    // Get coords of where this thread is in the global domain.
    int xs, ys, zs, xm, ym, zm, i, j, k;
    DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);

    // Fill in initial values.
    for (i = zs; i < zs+zm; i++)
    {
        for (j = ys; j < ys+ym; j++)
        {
            for (k = xs; k < xs+xm; k++)
            {
                phi->global_array[i][j][k] = 0.;
            }
        }
    }
    
    // Make a LinearSys
    linear_sys = new LinearSys(&da);
}

// Destructor
PoissonSolver::~PoissonSolver()
{
    delete linear_sys;
    delete phi;
    delete sigma;
    delete source;
    DMDestroy(&da);
}

// The code here is adapted from
// http://www.mcs.anl.gov/petsc/petsc-3.6/src/ksp/ksp/examples/tutorials/ex2.c.html
void PoissonSolver::run_solver(std::string output_file)
{
    // Make ghost rows available
    sigma->send_global_to_local();

    // Set up A and b
    int i,j,k,Ii,J,Istart,Iend;
    double v;
    
    // Currently, all PETSc parallel matrix formats are partitioned by
    // contiguous chunks of rows across the processors.  Determine which
    // rows of the matrix are locally owned.
    
    MatGetOwnershipRange(linear_sys->A,&Istart,&Iend);
    
    // Set matrix elements for the 3-D, seven-point stencil in parallel.
    // - Each processor needs to insert only elements that it owns
    // locally (but aNY non-local elements will be sent to the
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
        i = Ii/(NX*NY); j = (Ii - i*(NX*NY))/NX; k = Ii - j*NY - i*NX*NY;
        
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
            J = Ii - NX*NY;
            v = -coeffmzhalf;
            MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        if (i<NZ-1)
        {
            J = Ii + NX*NY;
            v = -coeffpzhalf;
            MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        if (j>0)
        {
            J = Ii - NX;
            v = -coeffmyhalf;
            MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        if (j<NY-1)
        {
            J = Ii + NX;
            v = -coeffpyhalf;
            MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        if (k>0)
        {
            J = Ii - 1;
            v = -coeffmxhalf;
            MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        if (k<NX-1)
        {
            J = Ii + 1;
            v = -coeffpxhalf;
            MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        
        // Periodic boundary conditions in y
        if (get_BC_type(Y_BC_TYPE) == DM_BOUNDARY_PERIODIC)
        {
            if (j == 0)
            {
                J = Ii + NX*NY - NX;
                v = -coeffmyhalf;
                MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
            if (j == NY - 1)
            {
                J = Ii - NX*NY + NX;
                v = -coeffpyhalf;
                MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
        }
        // Zero-flux boundary conditions in y. GHOSTED instead of MIRROR in 3D
        else if (get_BC_type(Y_BC_TYPE) == DM_BOUNDARY_GHOSTED)
        {
            if (j == 0)
            {
                J = Ii + NX;
                v = -coeffmyhalf-coeffpyhalf;
                MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
            if (j == NY - 1)
            {
                J = Ii - NX;
                v = -coeffmyhalf-coeffpyhalf;
                MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
        }
        else
        {
            PetscPrintf(PETSC_COMM_WORLD, "lateral BC problem!\n");
        }
        
        // Periodic boundary conditions in x
        if (get_BC_type(X_BC_TYPE) == DM_BOUNDARY_PERIODIC)
        {
            if (k == 0)
            {
                J = Ii + NX - 1;
                v = -coeffmxhalf;
                MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
            if (k == NX - 1)
            {
                J = Ii - NX + 1;
                v = -coeffpxhalf;
                MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
        }
        // Zero-flux boundary conditions in x. GHOSTED instead of MIRROR in 3D
        else if (get_BC_type(X_BC_TYPE) == DM_BOUNDARY_GHOSTED)
        {
            if (k == 0)
            {
                J = Ii + 1;
                v = -coeffmxhalf-coeffpxhalf;
                MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
            if (k == NX - 1)
            {
                J = Ii - 1;
                v = -coeffmxhalf-coeffpxhalf;
                MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
        }
        else
        {
            PetscPrintf(PETSC_COMM_WORLD, "lateral BC problem!\n");
        }
        
        // Derivative boundary conditions in z
        if (phi->bc.upper_BC_type == derivativeBC)
        {
            if (i == 0)
            {
                J = Ii + NX*NY;
                v = -coeffmzhalf-coeffpzhalf;
                MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
        }
        // Derivative boundary conditions in z
        if (phi->bc.lower_BC_type == derivativeBC)
        {
            if (i == NZ-1)
            {
                J = Ii - NX*NY;
                v = -coeffmzhalf-coeffpzhalf;
                MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
        }
        
        // Diagonal term
        v = coeffmxhalf + coeffmyhalf + coeffmzhalf + coeffpxhalf + coeffpyhalf + coeffpzhalf;
        
        MatSetValues(linear_sys->A,1,&Ii,1,&Ii,&v,INSERT_VALUES);
    }
    
    // Assemble matrix, using the 2-step process:
    // MatAssemblyBegin(), MatAssemblyEnd()
    // Computations can be done while messages are in transition
    // by placing code between these two statements.
    
    MatAssemblyBegin(linear_sys->A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(linear_sys->A,MAT_FINAL_ASSEMBLY);
    
    // Set right-hand-side vector.
    VecGetOwnershipRange(linear_sys->b,&Istart,&Iend);
    for (Ii=Istart; Ii<Iend; Ii++)
    {
        i = Ii/(NX*NY); j = (Ii - i*(NX*NY))/NX; k = Ii - j*NX - i*NX*NY;
        
        // Top side
        if (i == 0)
        {
            coeffmzhalf = 0.5*(sigma->local_array[i][j][k] + sigma->local_array[i-1][j][k])/(DELTA_X*DELTA_X);
            v = -source->global_array[i][j][k];
            
            if (phi->bc.upper_BC_type == constantBC)
            {
                v += coeffmzhalf*phi->bc.upper_BC_val;
            }
            else
            {
                v += 2.*DELTA_X*coeffmzhalf*phi->bc.upper_BC_val;
            }
        }
        // Bottom side
        else if (i == NZ - 1)
        {
            coeffpzhalf = 0.5*(sigma->local_array[i][j][k] + sigma->local_array[i+1][j][k]);
            v = -source->global_array[i][j][k];
            
            if (phi->bc.lower_BC_type == constantBC)
            {
                v += coeffpzhalf*phi->bc.lower_BC_val;
            }
            else
            {
                v += 2.*DELTA_X*coeffpzhalf*phi->bc.lower_BC_val;
            }
        }
        // Center points
        else
        {
            v = -source->global_array[i][j][k];
        }
        VecSetValue(linear_sys->b, Ii, v, INSERT_VALUES);
    }
    VecAssemblyBegin(linear_sys->b);
    VecAssemblyEnd(linear_sys->b);
    
    // Actually perform the solve
    KSPSetOperators(linear_sys->ksp,linear_sys->A,linear_sys->A);
    KSPSolve(linear_sys->ksp,linear_sys->b,phi->global_vec);
    
    // Save the solution to a file
    phi->write_to_file(output_file);
    
    return;
}