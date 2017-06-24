#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include "nonlocal_field.hpp"

class LinearSolver
{
    public:
        LinearSolver(DM *da, NonLocalField *phi, NonLocalField *sigma, Field *source, double DELTA_X);
        ~LinearSolver();
        void run_solver();
    
    private:
        NonLocalField *phi; // The field we are solving for
        NonLocalField *sigma;
        Field *source;
    
        Mat A;
        KSP ksp;
        PC pc;
        Vec b;
        int nx, ny, nz;
        double DELTA_X;
        DMBoundaryType x_BC_type, y_BC_type; // Lateral boundary condition types
};



#endif
