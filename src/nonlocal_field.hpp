#ifndef NONLOCAL_FIELD_H
#define NONLOCAL_FIELD_H

#include "field.hpp"

// Use an enum for boundary condition flags
enum BC_type 
{
    constantBC,
    derivativeBC
};

// Upper and lower boundary conditions and values
struct BC 
{
    BC_type lower_BC_type;
    double lower_BC_val;

    BC_type upper_BC_type;
    double upper_BC_val;
};

template <typename T>
class NonLocalField: public Field<T> // inherits from Field
{
    public:
        NonLocalField(DM *da, BC_type lower_BC_type, double lower_BC_val, \
                        BC_type upper_BC_type, double upper_BC_val);
        ~NonLocalField();
        void send_global_to_local();
    
        Vec local_vec;
        T local_array; // Either 2D or 3D
    
        BC bc;
};

#endif