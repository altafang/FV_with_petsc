#include <petscsys.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <string>
#include <map>
#include "IO_tools.hpp"
#include "nonlocal_field.hpp"
#include "linear_sys.hpp"
#include "poisson_solver_2D.hpp"

void PoissonSolver2D::read_input(const std::string &input_file)
{
    std::map<std::string, std::string> params;
    read_parameters(params, input_file);
    
    unpack(params, "NX", NX);
    unpack(params, "NY", NY);
    
    unpack(params, "DELTA_X", DELTA_X);
    
    unpack(params, "X_BC", X_BC);
    unpack(params, "Y_BC", Y_BC);
}

// Constructor
PoissonSolver2D::PoissonSolver2D(std::string input_file, std::string sigma_file, 
                                 std::string source_file)
{
    // Read input parameters from text file called "input.txt"
    // and unpack the input parameters.
    read_input(input_file);
    
    // Set up the distributed array used for all fields.
    DMDACreate2d(PETSC_COMM_WORLD, X_BC.get_DMBoundaryType(), 
                 Y_BC.get_DMBoundaryType(), DMDA_STENCIL_STAR, 
                 NX, NY, 1, PETSC_DECIDE, 1, 1, PETSC_NULL, PETSC_NULL, &da);
                 
    // Initialize with z boundary conditions
    phi = new NonLocalField<double**>("phi", &da, &X_BC, &Y_BC, NULL, DELTA_X);
    
    // Zero derivative boundary conditions in all directions for sigma
    sigma = new NonLocalField<double**>("sigma", &da, &zeroflux, &zeroflux, NULL, 
                                        DELTA_X);
    
    // Read sigma in from hdf5 file
    sigma->read_from_file(sigma_file);
    sigma->send_global_to_local();
    
    source = new Field<double**>("source", &da);

    // Read source term in from hdf5 file
    source->read_from_file(source_file);

    // Get coords of where this thread is in the global domain.
    int xs, ys, xm, ym;
    DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);

    // Fill in initial values.
    for (int i = ys; i < ys+ym; i++) {
        for (int j = xs; j < xs+xm; j++) {
            phi->global_array[i][j] = 0.;
        }
    }
    
    // Make a LinearSys
    linear_sys = new LinearSys(NX*NY, 5); // For 2D
}

// Destructor
PoissonSolver2D::~PoissonSolver2D()
{
    delete linear_sys;
    delete phi;
    delete sigma;
    delete source;
    DMDestroy(&da);
}

// The code here is adapted from
// http://www.mcs.anl.gov/petsc/petsc-3.6/src/ksp/ksp/examples/tutorials/ex2.c.html
void PoissonSolver2D::run_solver(std::string output_file)
{
    // Make ghost rows available
    sigma->send_global_to_local();

    // Set up A and b
    int i,j,J,Istart,Iend;
    double v;
    
    // Currently, all PETSc parallel matrix formats are partitioned by
    // contiguous chunks of rows across the processors.  Determine which
    // rows of the matrix are locally owned.
    
    MatGetOwnershipRange(linear_sys->A,&Istart,&Iend);
    
    // - Each processor needs to insert only elements that it owns
    // locally (but aNY non-local elements will be sent to the
    // appropriate processor during matrix assembly).
    // - Always specify global rows and columns of matrix entries.
    // 
    // Note: this uses the less common natural ordering that orders first
    // all the unknowns for x = h then for x = 2h etc; Hence you see J = Ii +- n
    // instead of J = I +- m as you might expect. The more standard ordering
    // would first do all variables for y = h, then y = 2h etc.
    
    double coeffpxhalf, coeffpyhalf, coeffmxhalf, coeffmyhalf;
    for (int Ii = Istart; Ii < Iend; ++Ii) {
        i = Ii/NX; j = Ii - i*NX;
        
        // Careful with signs: here, positive on the diagonals.
        // Retrieve coeff values, including neighboring values.
        coeffpxhalf = 0.5*(sigma->local_array[i][j] + sigma->local_array[i][j+1])
                      /(DELTA_X*DELTA_X);
        coeffmxhalf = 0.5*(sigma->local_array[i][j] + sigma->local_array[i][j-1])
                      /(DELTA_X*DELTA_X);
        coeffpyhalf = 0.5*(sigma->local_array[i][j] + sigma->local_array[i+1][j])
                      /(DELTA_X*DELTA_X);
        coeffmyhalf = 0.5*(sigma->local_array[i][j] + sigma->local_array[i-1][j])
                      /(DELTA_X*DELTA_X);
        
        // Fill in the finite-volume stencil
        if (i > 0) {
            J = Ii - NX;
            v = -coeffmyhalf;
            MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        if (i < NY-1) {
            J = Ii + NX;
            v = -coeffpyhalf;
            MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        if (j > 0) {
            J = Ii - 1;
            v = -coeffmxhalf;
            MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
        if (j < NX-1) {
            J = Ii + 1;
            v = -coeffpxhalf;
            MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
        }
                
        // Periodic boundary conditions in y. If one side is periodic, both are.
        if (Y_BC.lower_BC_type == BC_type::periodicBC) {
            if (i == 0) {
                J = Ii + NX*NY - NX;
                v = -coeffmyhalf;
                MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
            if (i == NY - 1) {
                J = Ii - NX*NY + NX;
                v = -coeffpyhalf;
                MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
        }
        
        // Periodic boundary conditions in x
        if (X_BC.lower_BC_type == BC_type::periodicBC) {
            if (j == 0) {
                J = Ii + NX - 1;
                v = -coeffmxhalf;
                MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
            if (j == NX - 1) {
                J = Ii - NX + 1;
                v = -coeffpxhalf;
                MatSetValues(linear_sys->A,1,&Ii,1,&J,&v,INSERT_VALUES);
            }
        }
                        
        // Diagonal term
        v = coeffmxhalf + coeffmyhalf + coeffpxhalf + coeffpyhalf;
                
        // Derivative boundary conditions in y. 
        if (Y_BC.lower_BC_type == BC_type::derivativeBC) {
            if (i == 0) {
                v -= coeffmyhalf;
            }
        }
        if (Y_BC.upper_BC_type == BC_type::derivativeBC) {                    
            if (i == NY - 1) {
                v -= coeffpyhalf;
            }
        }
                
        // Zero-flux boundary conditions in x. 
        if (X_BC.lower_BC_type == BC_type::derivativeBC) {
            if (j == 0) {
                v -= coeffmxhalf;
            }
        }
        if (X_BC.upper_BC_type == BC_type::derivativeBC) {
            if (j == NX - 1) {
                v -= coeffpxhalf;
            }
        }
        
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
    for (int Ii = Istart; Ii < Iend; ++Ii) {
        i = Ii/NX; j = Ii - i*NX;
        
        // Default
        v = -source->global_array[i][j];
        
        // Constant or derivative BCs for y
        if (i == 0) {
            coeffmyhalf = 0.5*(sigma->local_array[i][j] 
                          + sigma->local_array[i-1][j])/(DELTA_X*DELTA_X);
            
            if (Y_BC.lower_BC_type == BC_type::constantBC) {
                v += coeffmyhalf*Y_BC.lower_BC_val;
            } else if (Y_BC.lower_BC_type == BC_type::derivativeBC) {
                v += DELTA_X*coeffmyhalf*Y_BC.lower_BC_val;
            }
        } else if (i == NY - 1) {
            coeffpyhalf = 0.5*(sigma->local_array[i][j] 
                          + sigma->local_array[i+1][j])/(DELTA_X*DELTA_X);
            
            if (Y_BC.upper_BC_type == BC_type::constantBC) {
                v += coeffpyhalf*Y_BC.upper_BC_val;
            } else if (Y_BC.upper_BC_type == BC_type::derivativeBC) {
                v += DELTA_X*coeffpyhalf*Y_BC.upper_BC_val;
            }
        }
        // Constant or derivative BCs for x
        if (j == 0) {
            coeffmxhalf = 0.5*(sigma->local_array[i][j] 
                          + sigma->local_array[i][j-1])/(DELTA_X*DELTA_X);
            
            if (X_BC.lower_BC_type == BC_type::constantBC) {
                v += coeffmxhalf*X_BC.lower_BC_val;
            } else if (X_BC.lower_BC_type == BC_type::derivativeBC) {
                v += DELTA_X*coeffmxhalf*X_BC.lower_BC_val;
            }
        } else if (j == NX - 1) {
            coeffpxhalf = 0.5*(sigma->local_array[i][j] 
                          + sigma->local_array[i][j+1])/(DELTA_X*DELTA_X);
            
            if (X_BC.upper_BC_type == BC_type::constantBC) {
                v += coeffpxhalf*X_BC.upper_BC_val;
            } else if (X_BC.upper_BC_type == BC_type::derivativeBC) {
                v += DELTA_X*coeffpxhalf*X_BC.upper_BC_val;
            }
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
}