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

// constructor
NonLocalField::NonLocalField(DM *da, BC *bc): Field(da), bc(bc)
{
    DMCreateLocalVector(*da, &local_vec);
    DMDAVecGetArray(*da, local_vec, &local_array);
}

// destructor
NonLocalField::~NonLocalField()
{
    DMDAVecRestoreArray(*da, local_vec, &local_array);
    VecDestroy(&local_vec);
}

void NonLocalField::send_global_to_local()
{
    DMGlobalToLocalBegin(*da, global_vec, INSERT_VALUES, local_vec);
    DMGlobalToLocalEnd(*da, global_vec, INSERT_VALUES, local_vec);
    
    // Fill in boundary conditions
    int xs, ys, zs, xm, ym, zm, j, k;
    int nx, ny, nz;
    DMDAGetInfo(*da, NULL, &nx, &ny, &nz, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    DMDAGetCorners(*da, &xs, &ys, &zs, &xm, &ym, &zm);
    if (zs == 0)
    {
        // derivative BC
        if (bc->upper_BC_type == derivativeBC)
        {
            for (j = ys; j < ys+ym; j++)
            {
                for (k = xs; k < xs+xm; k++)
                {
                    local_array[-1][j][k] = global_array[1][j][k] - 2.*bc->upper_BC_val;
                }
            }
        }
        else // const BC
        {
            for (j = ys; j < ys+ym; j++)
            {
                for (k = xs; k < xs+xm; k++)
                {
                    local_array[-1][j][k] = bc->upper_BC_val;
                }
            }
        }
    }
    if (zs + zm == nz)
    {
        // derivative BC
        if (bc->lower_BC_type == derivativeBC)
        {
            for (j = ys; j < ys+ym; j++)
            {
                for (k = xs; k < xs+xm; k++)
                {
                    local_array[nz][j][k] = global_array[nz-2][j][k] + 2.*bc->lower_BC_val;
                }
            }
        }
        else // const BC
        {
            for (j = ys; j < ys+ym; j++)
            {
                for (k = xs; k < xs+xm; k++)
                {
                    local_array[nz][j][k] = bc->lower_BC_val;
                }
            }
        }
    }
    
    // Use zero-derivative boundary conditions.
    // DM_BOUNDARY_MIRROR is not yet implemented in 3D so have to manually fill those cells...
    if (ys == 0)
    {
        for (j = zs; j < zs+zm; j++)
        {
            for (k = xs; k < xs+xm; k++)
            {
                local_array[j][-1][k] = global_array[j][1][k];
            }
        }
    }
    if (ys + ym == ny)
    {
        for (j = zs; j < zs+zm; j++)
        {
            for (k = xs; k < xs+xm; k++)
            {
                local_array[j][ny][k] = global_array[j][ny-2][k];
            }
        }
    }
    
    if (xs == 0)
    {
        for (j = zs; j < zs+zm; j++)
        {
            for (k = ys; k < ys+ym; k++)
            {
                local_array[j][k][-1] = global_array[j][k][1];
            }
        }
    }
    if (xs + xm == nx)
    {
        for (j = zs; j < zs+zm; j++)
        {
            for (k = ys; k < ys+ym; k++)
            {
                local_array[j][k][nx] = global_array[j][k][nx-2];
            }
        }
    }
    
    return;
}

#endif