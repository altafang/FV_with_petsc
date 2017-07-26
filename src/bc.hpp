#ifndef BC_H
#define BC_H

// For simplicity, use strings instead of enum for now
// Use an enum for boundary condition flags
// enum BC_type 
// {
//     constantBC,
//     derivativeBC,
//     periodicBC
// };

// Upper and lower boundary conditions and values
struct BC 
{
    std::string lower_BC_type;
    double lower_BC_val;

    std::string upper_BC_type;
    double upper_BC_val;
    
    // Default is zero-flux
    BC() : lower_BC_type("derivative"), lower_BC_val(0), upper_BC_type("derivative"), upper_BC_val(0) {}
};

#endif