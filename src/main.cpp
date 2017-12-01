// Alta Fang, 2017
//
// Solve Poisson's equation using PETSc.
// 
// To run on two processors, execute for example:
//    mpiexec -np 2 ./solve_poisson

#include "poisson_solver_3D.hpp"

int main(int argc,char **args)
{
    PetscInitialize(&argc,&args,NULL,NULL);

    // Create a PoissonSolver
    PoissonSolver3D *solver = new PoissonSolver3D();
    
    // Run the solver
    solver->run_solver();
    
    // Cleanup
    delete solver;
    
    PetscFinalize();
    
    return 0;
}
