// A Field does not have boundary conditions, so it only has a global_array

#ifndef FIELD_H
#define FIELD_H

#include <petscsys.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <string>

template <typename T>
class Field
{
    public:
        Field(DM *da);
        virtual ~Field();
        void write_to_file(const std::string &filename);
        void read_from_file(const std::string &filename);
    
        Vec global_vec;
        T global_array; // T is double*** for 3D and double** for 2D
        DM *da; // Pointer to DM
};

#endif