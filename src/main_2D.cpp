// Alta Fang, 2017
//
// Solve Poisson's equation using PETSc.
// 
// To run on two processors, execute for example:
//    mpiexec -np 2 ./solve_poisson_2D

#include "poisson_solver_2D.hpp"

int main(int argc,char **args)
{
    PetscInitialize(&argc,&args,NULL,NULL);

    // Create a PoissonSolver
    PoissonSolver2D *solver = new PoissonSolver2D();
    
    // Run the solver
    solver->run_solver();
    
    // Cleanup
    delete solver;
    
    PetscFinalize();
    
    return 0;
}
