// A Field does not have boundary conditions, so it only has a global_array

#ifndef FIELD_H
#define FIELD_H

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include "petscviewerhdf5.h"
#include <iostream>
#include <string>
#include <iomanip>
#include "IO_tools.hpp"

class Field
{
    public:
        Field(DM *da);
        ~Field();
        void write_to_file(std::string filename);
        void read_from_file(std::string filename);
    
        Vec global_vec;
        double ***global_array;
        DM *da; // Pointer to DM
};

#endif