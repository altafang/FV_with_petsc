#include <petscsys.h>
#include <string>
#include <iostream>
#include <iomanip>
#include "petscviewerhdf5.h"
#include "field.hpp"

// constructor
template <typename T>
Field<T>::Field(const std::string &name, DM *da): da(da)
{
    DMCreateGlobalVector(*da, &global_vec);
    VecSet(global_vec, 0.); // zero out the Vec to begin with.
    DMDAVecGetArray(*da, global_vec, &global_array);
    PetscObjectSetName(reinterpret_cast<PetscObject>(global_vec), name.c_str());
}

// destructor
template <typename T>
Field<T>::~Field()
{
    DMDAVecRestoreArray(*da, global_vec, &global_array);
    VecDestroy(&global_vec);
}

// Write field to an hdf5 file
template <typename T>
void Field<T>::write_to_file(const std::string &filename)
{
    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &viewer);
    VecView(global_vec, viewer);
    PetscViewerDestroy(&viewer);
}

// Read field from an hdf5 file
template <typename T>
void Field<T>::read_from_file(const std::string &filename)
{
    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(global_vec, viewer);
    PetscViewerDestroy(&viewer);
}

// Specialization of template to specific types
// 3D fields
template Field<double***>::Field(const std::string &name, DM *da);
template Field<double***>::~Field();
template void Field<double***>::write_to_file(const std::string &filename);
template void Field<double***>::read_from_file(const std::string &filename);
// 2D fields
template Field<double**>::Field(const std::string &name, DM *da);
template Field<double**>::~Field();
template void Field<double**>::write_to_file(const std::string &filename);
template void Field<double**>::read_from_file(const std::string &filename);

