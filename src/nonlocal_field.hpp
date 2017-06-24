#ifndef NONLOCAL_FIELD_H
#define NONLOCAL_FIELD_H

#include "field.hpp"

// Use an enum for boundary condition flags
enum BC_type {
    constantBC,
    derivativeBC
};

// Upper and lower boundary conditions and values
typedef struct {
    BC_type upper_BC_type;
    double upper_BC_val;
    
    BC_type lower_BC_type;
    double lower_BC_val;
} BC;

class NonLocalField: public Field // inherits from Field
{
    public:
        NonLocalField(DM *da, BC *bc);
        ~NonLocalField();
        void send_global_to_local();
    
        Vec local_vec;
        double ***local_array;
    
        BC *bc; // Pointer to BC
};

#endif