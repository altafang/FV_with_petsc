#ifndef POISSON_SOLVER_3D_H
#define POISSON_SOLVER_3D_H

#include "IO_tools.hpp"
#include "nonlocal_field.hpp"
#include "linear_sys.hpp"
#include <string>

class PoissonSolver3D
{
    public:    
        PoissonSolver3D(std::string input_file="input.txt", 
                        std::string sigma_file="sigma.h5", 
                        std::string source_file="source.h5");
        ~PoissonSolver3D();
        void read_input(const std::string &input_file);
        void run_solver(std::string output_file="phi.h5");
        
        DM da;
        NonLocalField<double***> *phi; // The field we are solving for
        NonLocalField<double***> *sigma;
        Field<double***> *source;
        
    private:
        LinearSys *linear_sys;
        int NX, NY, NZ;
        double DELTA_X;
        BC X_BC, Y_BC, Z_BC; 
};

#endif