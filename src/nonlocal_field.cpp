#include "nonlocal_field.hpp"

// constructor
template <typename T>
NonLocalField<T>::NonLocalField(DM *da, BC *x_bc, BC *y_bc, BC *z_bc, double DELTA_X) : Field<T>(da), \
    x_bc(x_bc), y_bc(y_bc), z_bc(z_bc), DELTA_X(DELTA_X)
{
    DMCreateLocalVector(*da, &local_vec);
    DMDAVecGetArray(*da, local_vec, &local_array);
}

// destructor
template <typename T>
NonLocalField<T>::~NonLocalField()
{
    // Notice that you have to use 'this' to access base class variables when using templates.
    // See https://stackoverflow.com/questions/4643074/why-do-i-have-to-access-template-base-class-members-through-the-this-pointer
    DMDAVecRestoreArray(*(this->da), local_vec, &local_array);
    VecDestroy(&local_vec);
}

// Specialization of template to specific types
// 3D fields
template NonLocalField<double***>::NonLocalField(DM *da, BC *x_bc, BC *y_bc, BC *z_bc, double DELTA_X);
template NonLocalField<double***>::~NonLocalField();
// 2D fields
template NonLocalField<double**>::NonLocalField(DM *da, BC *x_bc, BC *y_bc, BC *z_bc, double DELTA_X);
template NonLocalField<double**>::~NonLocalField();

// Implementation specific to 3D
template<>
void NonLocalField<double***>::send_global_to_local()
{
    DMGlobalToLocalBegin(*(this->da), global_vec, INSERT_VALUES, local_vec);
    DMGlobalToLocalEnd(*(this->da), global_vec, INSERT_VALUES, local_vec);
    
    // Fill in boundary conditions
    int xs, ys, zs, xm, ym, zm;
    int nx, ny, nz;
    DMDAGetInfo(*(this->da), NULL, &nx, &ny, &nz, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    DMDAGetCorners(*(this->da), &xs, &ys, &zs, &xm, &ym, &zm);
    // z BCs
    if (zs == 0)
    {
        // derivative BC
        if (z_bc->lower_BC_type == derivativeBC)
        {
            for (int j = ys; j < ys+ym; ++j)
            {
                for (int k = xs; k < xs+xm; ++k)
                {
                    local_array[-1][j][k] = global_array[1][j][k] - 2.*z_bc->lower_BC_val;
                }
            }
        }
        else if (z_bc->lower_BC_type == constantBC) // const BC
        {
            for (int j = ys; j < ys+ym; ++j)
            {
                for (int k = xs; k < xs+xm; ++k)
                {
                    local_array[-1][j][k] = z_bc->lower_BC_val;
                }
            }
        }
    }
    if (zs + zm == nz)
    {
        // derivative BC
        if (z_bc->upper_BC_type == derivativeBC)
        {
            for (int j = ys; j < ys+ym; ++j)
            {
                for (int k = xs; k < xs+xm; ++k)
                {
                    local_array[nz][j][k] = global_array[nz-2][j][k] + 2.*DELTA_X*z_bc->upper_BC_val;
                }
            }
        }
        else if (z_bc->upper_BC_type == constantBC) // const BC
        {
            for (int j = ys; j < ys+ym; ++j)
            {
                for (int k = xs; k < xs+xm; ++k)
                {
                    local_array[nz][j][k] = z_bc->upper_BC_val;
                }
            }
        }
    }
    
    if (ys == 0)
    {
        if (y_bc->lower_BC_type == derivativeBC)
        {
            for (int j = zs; j < zs+zm; ++j)
            {
                for (int k = xs; k < xs+xm; ++k)
                {
                    local_array[j][-1][k] = global_array[j][1][k] - 2.*DELTA_X*y_bc->lower_BC_val;
                }
            }
        }
        else if (y_bc->lower_BC_type == constantBC)
        {
            for (int j = zs; j < zs+zm; ++j)
            {
                for (int k = xs; k < xs+xm; ++k)
                {
                    local_array[j][-1][k] = y_bc->lower_BC_val;
                }
            }
        }
    }
    if (ys + ym == ny)
    {
        if (y_bc->upper_BC_type == derivativeBC)
        {
            for (int j = zs; j < zs+zm; ++j)
            {
                for (int k = xs; k < xs+xm; ++k)
                {
                    local_array[j][ny][k] = global_array[j][ny-2][k] + 2.*DELTA_X*y_bc->upper_BC_val;
                }
            }
        }
        else if (y_bc->upper_BC_type == constantBC)
        {
            for (int j = zs; j < zs+zm; ++j)
            {
                for (int k = xs; k < xs+xm; ++k)
                {
                    local_array[j][ny][k] = y_bc->upper_BC_val;
                }
            }
        }
    }
    
    if (xs == 0)
    {
        if (x_bc->lower_BC_type == derivativeBC)
        {
            for (int j = zs; j < zs+zm; ++j)
            {
                for (int k = ys; k < ys+ym; ++k)
                {
                    local_array[j][k][-1] = global_array[j][k][1] - 2.*DELTA_X*x_bc->lower_BC_val;
                }
            }
        }
        else if (x_bc->lower_BC_type == constantBC)
        {
            for (int j = zs; j < zs+zm; ++j)
            {
                for (int k = ys; k < ys+ym; ++k)
                {
                    local_array[j][k][-1] = x_bc->lower_BC_val;
                }
            }
        }
    }
    if (xs + xm == nx)
    {
        if (x_bc->upper_BC_type == derivativeBC)
        {
            for (int j = zs; j < zs+zm; ++j)
            {
                for (int k = ys; k < ys+ym; ++k)
                {
                    local_array[j][k][nx] = global_array[j][k][nx-2] + 2.*DELTA_X*x_bc->upper_BC_val;
                }
            }
        }
        else if (x_bc->upper_BC_type == constantBC)
        {
            for (int j = zs; j < zs+zm; ++j)
            {
                for (int k = ys; k < ys+ym; ++k)
                {
                    local_array[j][k][nx] = x_bc->upper_BC_val;
                }
            }
        }
    }
    
    return;
}

// Implementation specific to 2D
template<>
void NonLocalField<double**>::send_global_to_local()
{
    DMGlobalToLocalBegin(*(this->da), global_vec, INSERT_VALUES, local_vec);
    DMGlobalToLocalEnd(*(this->da), global_vec, INSERT_VALUES, local_vec);
    
    /* Fill in upper and lower boundary conditions */
    int xs, ys, xm, ym, nx, ny;
    DMDAGetInfo(*(this->da), NULL, &nx, &ny, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    DMDAGetCorners(*(this->da), &xs, &ys, NULL, &xm, &ym, NULL);
        
    if (ys == 0)
    {
        if (y_bc->lower_BC_type == derivativeBC)
        {
            for (int k = xs; k < xs+xm; ++k)
            {
                local_array[-1][k] = global_array[1][k] - 2.*DELTA_X*y_bc->lower_BC_val;
            }
        }
        else if (y_bc->lower_BC_type == constantBC)
        {
            for (int k = xs; k < xs+xm; ++k)
            {
                local_array[-1][k] = y_bc->lower_BC_val;
            }
        }
    }
    if (ys + ym == ny)
    {
        if (y_bc->upper_BC_type == derivativeBC)
        {
            for (int k = xs; k < xs+xm; ++k)
            {
                local_array[ny][k] = global_array[ny-2][k] + 2.*DELTA_X*y_bc->upper_BC_val;
            }
        }
        else if (y_bc->upper_BC_type == constantBC)
        {
            for (int k = xs; k < xs+xm; ++k)
            {
                local_array[ny][k] = y_bc->upper_BC_val;
            }
        }
    }
    
    if (xs == 0)
    {
        if (x_bc->lower_BC_type == derivativeBC)
        {
            for (int k = ys; k < ys+ym; ++k)
            {
                local_array[k][-1] = global_array[k][1] - 2.*DELTA_X*x_bc->lower_BC_val;
            }
        }
        else if (x_bc->lower_BC_type == constantBC)
        {
            for (int k = ys; k < ys+ym; ++k)
            {
                local_array[k][-1] = x_bc->lower_BC_val;
            }
        }
    }
    if (xs + xm == nx)
    {
        if (x_bc->upper_BC_type == derivativeBC)
        {
            for (int k = ys; k < ys+ym; ++k)
            {
                local_array[k][nx] = global_array[k][nx-2] + 2.*DELTA_X*x_bc->upper_BC_val;
            }
        }
        else if (x_bc->upper_BC_type == constantBC)
        {
            for (int k = ys; k < ys+ym; ++k)
            {
                local_array[k][nx] = x_bc->upper_BC_val;
            }
        }
    }
    return;
}




