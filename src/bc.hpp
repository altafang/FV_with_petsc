#ifndef BC_H
#define BC_H

// Use an enum for boundary condition flags
enum BC_type 
{
    constantBC,
    derivativeBC,
    periodicBC
};

// Upper and lower boundary conditions and values
struct BC 
{
    BC_type lower_BC_type;
    double lower_BC_val;

    BC_type upper_BC_type;
    double upper_BC_val;
    
    // Default is zero-flux
    BC() : lower_BC_type(derivativeBC), lower_BC_val(0), upper_BC_type(derivativeBC), upper_BC_val(0) {}
};

#endif