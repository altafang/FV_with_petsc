#include "field.hpp"

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
void Field::write_to_file(std::string filename)
{
    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &viewer);
    VecView(global_vec, viewer);
    PetscViewerDestroy(&viewer);
}

// Read field from an hdf5 file
void Field::read_from_file(std::string filename)
{
    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &viewer);
    VecLoad(global_vec, viewer);
    PetscViewerDestroy(&viewer);
}