#ifndef NONLOCAL_FIELD_H
#define NONLOCAL_FIELD_H

#include <petscsys.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <string>
#include "bc.hpp"
#include "field.hpp"

template <typename T>
class NonLocalField: public Field<T> // inherits from Field
{
    public:
        NonLocalField(const std::string &name, DM *da, BC *x_bc, BC *y_bc, BC *z_bc, 
                      double DELTA_X);
        ~NonLocalField();
        void send_global_to_local();
    
        Vec local_vec;
        T local_array; // Either 2D or 3D
    
        BC *x_bc;
        BC *y_bc;
        BC *z_bc; // In 2D, z_bc will just be NULL
        
        double DELTA_X;
};

#endif