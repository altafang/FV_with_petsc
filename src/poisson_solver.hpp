#ifndef POISSON_SOLVER_H
#define POISSON_SOLVER_H

#include "IO_tools.hpp"
#include "nonlocal_field.hpp"
#include "linear_sys.hpp"
#include <string>

class PoissonSolver
{
    public:    
        PoissonSolver(std::string input_file="input.txt");
        ~PoissonSolver();
        void run_solver(std::string output_file="phi.h5");
        
        DM da;
        NonLocalField *phi; // The field we are solving for
        NonLocalField *sigma;
        Field *source;
        
    private:
        LinearSys *linear_sys;
        int NX, NY, NZ;
        double DELTA_X;
        double PHI_UPPER, PHI_LOWER;
        std::string x_BC_type, y_BC_type;
};

#endif