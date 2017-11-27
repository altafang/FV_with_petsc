#ifndef POISSON_SOLVER_2D_H
#define POISSON_SOLVER_2D_H

#include "nonlocal_field.hpp"
#include "linear_sys.hpp"
#include <string>

class PoissonSolver2D
{
    public:    
        PoissonSolver2D(std::string input_file="input.txt", 
                        std::string sigma_file="sigma.h5", 
                        std::string source_file="source.h5");
        ~PoissonSolver2D();
        void read_input(const std::string &input_file);
        void run_solver(std::string output_file="phi.h5");
        
        DM da;
        NonLocalField<double**> *phi; // The field we are solving for
        NonLocalField<double**> *sigma;
        Field<double**> *source;
        
    private:
        LinearSys *linear_sys;
        int NX, NY;
        double DELTA_X;
        BC X_BC, Y_BC;
};

#endif