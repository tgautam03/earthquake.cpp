#include <iostream>

#include "./include/Wave1D.hpp"

int main(int argc, char const *argv[])
{
    // Space discretization
    float y_start = 0; float y_end = 20; unsigned int ny = 2001;

    // Time discretization
    float t_start = 0; float t_end = 10; unsigned int nt = 10001;

    // Medium velocity
    float c_air = 2.0f;
    float c_ground = 5.0f;
    int interface_loc = ny / 4;
    float *c = new float[ny];
    for (unsigned int i = 0; i < ny; i++)
    {
        if (i < interface_loc)
            c[i] = c_air;
        else
            c[i] = c_ground;
    }
        

    // Required resources
    double sizeof_float = sizeof(float);
    std::cout << "Need atleast: " << (ny*nt*sizeof_float + ny*sizeof_float) / 1073741824.0 << " GB \n \n";

    // Define 1D wave simulation settings 
    Wave1D simulation = Wave1D(y_start, y_end, ny, t_start, t_end, nt);

    // Define source function
    simulation.add_src_dGaussian(0.15f, 5.0f, (y_end-y_start)/2.0f);
    
    // Solve using 3 point stencil with top:0 & bottom:0 boundary conditions
    simulation.solve_3pt_00(c);
    simulation.solve_3pt_da(c, 0.75);
    simulation.solve_5pt_da(c, 0.75);
    simulation.solve_7pt_da(c, 0.75);

    // Free dynamic array
    delete[] c;

    return 0;
}
