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

class Field
{
    public:
        Field(DM *da);
        ~Field();
        void write_to_file(std::string filename, int snapshot_counter);
        void read_from_file(std::string filename, int snapshot_counter);
    
        Vec global_vec;
        double ***global_array;
        DM *da; // Pointer to DM
};

// constructor
Field::Field(DM *da): da(da)
{
    DMCreateGlobalVector(*da, &global_vec);
    VecSet(global_vec, 0.); // zero out the Vec to begin with.
    
    DMDAVecGetArray(*da, global_vec, &global_array);
}

// destructor
Field::~Field()
{
    DMDAVecRestoreArray(*da, global_vec, &global_array);
    VecDestroy(&global_vec);
}

// Write field to an hdf5 file
void Field::write_to_file(std::string filename, int snapshot_counter)
{
    PetscViewer viewer;
    
    std::stringstream index_string;
    index_string.str("");
    index_string.clear();
    index_string << std::setfill('0') << std::setw(3) << snapshot_counter; // %03d formatting
    std::string complete_filename = filename + "_" + index_string.str() + ".h5";
    
    PetscViewerHDF5Open(PETSC_COMM_WORLD, complete_filename.c_str(), FILE_MODE_WRITE, &viewer);
    VecView(global_vec, viewer);
    PetscViewerDestroy(&viewer);
    return;
}

// Read field from an hdf5 file
void Field::read_from_file(std::string filename, int snapshot_counter)
{
    PetscViewer viewer;
    
    std::stringstream index_string;
    index_string.str("");
    index_string.clear();
    index_string << std::setfill('0') << std::setw(3) << snapshot_counter; // %03d formatting
    std::string complete_filename = filename + "_" + index_string.str() + ".h5";
    
    PetscViewerHDF5Open(PETSC_COMM_WORLD, complete_filename.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(global_vec, viewer);
    PetscViewerDestroy(&viewer);
    return;
}

#endif