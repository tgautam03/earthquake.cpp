#include <iostream>

#include "./include/Wave3D.hpp"

int main(int argc, char const *argv[])
{
    // Space discretization
    float x_start = 0; float x_end = 10; unsigned long nx = 201;
    float y_start = 0; float y_end = 10; unsigned long ny = 201;
    float z_start = 0; float z_end = 10; unsigned long nz = 201;

    // Time discretization
    float t_start = 0; float t_end = 3; unsigned long nt = 3001;

    // Required resources
    double sizeof_float = sizeof(float);
    std::cout << "Need atleast: " << (4*nx*ny*nz*sizeof_float + nx*nz*nt*sizeof_float + nx*ny*nt*sizeof_float) / 1073741824.0 << " GB \n \n";

    // Medium velocity
    Tensor3dFP32 c = Tensor3dFP32(nx,ny,nz);
    float c_air = 2.0f;
    float c_ground = 5.0f;
    int interface_loc = ny / 4;
    for (unsigned long iz = 0; iz < nz; iz++)
    {
        for (unsigned long iy = 0; iy < ny; iy++)
        {
            for (unsigned long ix = 0; ix < nx; ix++)
            {
                if (iy < interface_loc)
                    c(ix,iy,iz) = c_air;
                else
                    c(ix,iy,iz) = c_ground;
            }
        }
    }

    Wave3D simulation = Wave3D(x_start, x_end, nx, 
                                y_start, y_end, ny, 
                                z_start, z_end, nz, 
                                t_start, t_end, nt);

    // Define source function
    simulation.add_src_dGaussian(0.1f, 5.0f, x_end/2, y_end/1.25, z_end/2);
    
    // Solve
    simulation.solve_slices(c, 0.5, z_end/2);

    return 0;
}
