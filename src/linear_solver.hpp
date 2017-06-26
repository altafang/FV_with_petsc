#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include "nonlocal_field.hpp"

class LinearSolver
{
    public:
        LinearSolver(DM *da, NonLocalField *phi, NonLocalField *sigma, Field *source);
        ~LinearSolver();
        void run_solver(double DELTA_X);
    
    private:
        DM *da;
        NonLocalField *phi; // The field we are solving for
        NonLocalField *sigma;
        Field *source;
    
        Mat A;
        KSP ksp;
        PC pc;
        Vec b;
};



#endif
