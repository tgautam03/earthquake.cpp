#include "../include/Wave1D.hpp"

Wave1D::Wave1D(float ys_, float ye_, unsigned long ny_, 
            float ts_, float te_, float nt_): ys(ys_), ye(ye_), ny(ny_), ts(ts_), te(te_), nt(nt_),
                    dy((ye - ys) / static_cast<float>(ny-1)), 
                    dt((te - ts) / static_cast<float>(nt-1))
{
    // Defining y and t grid
    y = new float[ny];
    linspace_1d(y, ny, ys, dy);

    t = new float[nt];
    linspace_1d(t, nt, ts, dt);
}

/*=================================================
================ Source function ==================
=================================================*/
void Wave1D::add_src_dGaussian(float t0, float f0, float src_loc)
{
    // Source location on the grid point
    i_src = static_cast<unsigned long>(src_loc / dy) + 1;
    
    // Defining source function
    src = new float[nt];
    for (int i = 0; i < nt; i++)
    {
        float t_i = t[i];
        src[i] = -8 * (t_i - t0) * f0 * exp(-1 * 4*f0*(t_i - t0) * 4*f0*(t_i - t0));
    }
}

/*=================================================
================ 3 point stencil ==================
=================================================*/
void Wave1D::solve_3pt_00(float *c)
{
    std::cout << "Simulating using 3 point stencil; top boundary: 0; bottom boundary: 0 \n";

    // Checking CFL condition
    float max_c = 0;
    float cur_c;
    for (int i = 0; i < ny; i++)
    {
        cur_c = c[i];
        if (cur_c > max_c)
            max_c = cur_c;
    }
    assert((max_c * dt/dy) < 1 && "CFL condition not satisfied!");

    // Initialise Solution matrix
    MatrixFP32 u = MatrixFP32(ny, nt);

    auto start = std::chrono::high_resolution_clock::now();
    // Loop over time
    for (unsigned long it = 1; it < nt-1; it++)
    {
        // Loop over space
        for (unsigned long iy = 1; iy < ny-1; iy++)
        {
            // 2nd derivative wrt x
            float d2u_dx2 = (u(iy+1,it) - 2*u(iy,it) + u(iy-1,it)) / (dy*dy);

            // Updating solution
            if (iy == i_src)
                u(iy,it+1) = (c[iy]*dt * c[iy]*dt) * d2u_dx2 + 2*u(iy,it) - u(iy,it-1) + dt*dt * src[it] / dy;
            else
                u(iy,it+1) = (c[iy]*dt * c[iy]*dt) * d2u_dx2 + 2*u(iy,it) - u(iy,it-1);
        }

        // Top: Zero boundary
        u(0,it+1) = 0; 
        
        // Bottom: Zero boundary
        u(ny-1,it+1) = 0;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Simulation took: " << duration.count() << " microseconds" << std::endl;

    // Saving simulation details
    // Saving solution as bin
    std::cout << "Saving results and grid details to: data/1d/ \n";
    u.save("data/1d/sol_3pt_00.bin");
    // Saving grid as bin
    save_vec(y, ny, "data/1d/y_3pt_00.bin");
    save_vec(t, nt, "data/1d/t_3pt_00.bin");
    save_vec(c, ny, "data/1d/c_3pt_00.bin");
    std::cout << "Done \n \n";
}

void Wave1D::solve_3pt_da(float *c, float alpha)
{
    std::cout << "Simulating using 3 point stencil; top boundary: damping; bottom boundary: absorbing \n";

    // Checking CFL condition
    float max_c = 0;
    float cur_c;
    for (int i = 0; i < ny; i++)
    {
        cur_c = c[i];
        if (cur_c > max_c)
            max_c = cur_c;
    }
    assert((max_c * dt/dy) < 1 && "CFL condition not satisfied!");

    // Initialise Solution matrix
    MatrixFP32 u = MatrixFP32(ny, nt);

    auto start = std::chrono::high_resolution_clock::now();
    // Loop over time
    for (unsigned long it = 1; it < nt-1; it++)
    {
        // Loop over space
        for (unsigned long iy = 1; iy < ny-1; iy++)
        {
            // 2nd derivative wrt x
            float d2u_dx2 = (u(iy+1,it) - 2*u(iy,it) + u(iy-1,it)) / (dy*dy);

            // Updating solution
            if (iy == i_src)
                u(iy,it+1) = (c[iy]*dt * c[iy]*dt) * d2u_dx2 + 2*u(iy,it) - u(iy,it-1) + dt*dt * src[it] / dy;
            else
                u(iy,it+1) = (c[iy]*dt * c[iy]*dt) * d2u_dx2 + 2*u(iy,it) - u(iy,it-1);
        }

        // Top: Damping boundary
        u(0,it+1) = u(0,it) + (c[0]*dt/dy)*(u(1,it) - u(0,it)) - alpha*(u(0,it) - u(0,it-1));
        
        // Bottom: Absorbing boundary
        u(ny-1,it+1) = u(ny-2,it) + (c[ny-1]*dt-dy)/(c[ny-1]*dt+dy) * (u(ny-2,it+1)-u(ny-1,it));

        /*
        // Top: Absorbing boundary
        u(0,it+1) = u(1,it) + (c[0]*dt-dy)/(c[0]*dt+dy) * (u(1,it+1)-u(0,it));
        
        // Bottom: Damping boundary
        u(ny-1,it+1) = u(ny-1,it) + (c[ny-1]*dt/dy)*(u(ny-2,it) - u(ny-1,it)) - alpha*(u(ny-1,it) - u(ny-1,it-1));
        */
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Simulation took: " << duration.count() << " microseconds" << std::endl;

    // Saving simulation details
    // Saving solution as bin
    std::cout << "Saving results and grid details to: data/1d/ \n";
    u.save("data/1d/sol_3pt_da.bin");
    // Saving grid as bin
    save_vec(y, ny, "data/1d/y_3pt_da.bin");
    save_vec(t, nt, "data/1d/t_3pt_da.bin");
    save_vec(c, ny, "data/1d/c_3pt_da.bin");
    std::cout << "Done \n \n";
}

/*=================================================
================ 5 point stencil ==================
=================================================*/
void Wave1D::solve_5pt_da(float *c, float alpha)
{
    std::cout << "Simulating using 5 point stencil; top boundary: damping; bottom boundary: absorbing \n";

    // Checking CFL condition
    float max_c = 0;
    float cur_c;
    for (int i = 0; i < ny; i++)
    {
        cur_c = c[i];
        if (cur_c > max_c)
            max_c = cur_c;
    }
    assert((max_c * dt/dy) < 1 && "CFL condition not satisfied!");

    // Initialise Solution matrix
    MatrixFP32 u = MatrixFP32(ny, nt);

    auto start = std::chrono::high_resolution_clock::now();
    // Loop over time
    for (unsigned long it = 1; it < nt-1; it++)
    {
        // Loop over space
        for (unsigned long iy = 2; iy < ny-2; iy++)
        {
            // 2nd derivative wrt x (5 point approximation)
            float d2u_dx2 = (-1/12*u(iy+2,it) + 4/3*u(iy+1,it) - 5/2*u(iy,it) + 4/3*u(iy-1,it) - 1/12*u(iy-2,it))/(dy*dy);

            // Updating solution
            if (iy == i_src)
                u(iy,it+1) = (c[iy]*dt * c[iy]*dt) * d2u_dx2 + 2*u(iy,it) - u(iy,it-1) + dt*dt * src[it] / dy;
            else
                u(iy,it+1) = (c[iy]*dt * c[iy]*dt) * d2u_dx2 + 2*u(iy,it) - u(iy,it-1);
        }

        // Top: Damping boundary
        u(1,it+1) = u(1,it) + (c[1]*dt/dy)*(u(2,it) - u(1,it)) - alpha*(u(1,it) - u(1,it-1));
        u(0,it+1) = u(0,it) + (c[0]*dt/dy)*(u(1,it) - u(0,it)) - alpha*(u(0,it) - u(0,it-1));
        
        // Bottom: Absorbing boundary
        u(ny-2,it+1) = u(ny-3,it) + (c[ny-2]*dt-dy)/(c[ny-2]*dt+dy) * (u(ny-3,it+1)-u(ny-2,it));
        u(ny-1,it+1) = u(ny-2,it) + (c[ny-1]*dt-dy)/(c[ny-1]*dt+dy) * (u(ny-2,it+1)-u(ny-1,it));
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Simulation took: " << duration.count() << " microseconds" << std::endl;

    // Saving simulation details
    // Saving solution as bin
    std::cout << "Saving results and grid details to: data/1d/ \n";
    u.save("data/1d/sol_5pt_da.bin");
    // Saving grid as bin
    save_vec(y, ny, "data/1d/y_5pt_da.bin");
    save_vec(t, nt, "data/1d/t_5pt_da.bin");
    save_vec(c, ny, "data/1d/c_5pt_da.bin");
    std::cout << "Done \n \n";
}

/*=================================================
================ 7 point stencil ==================
=================================================*/
void Wave1D::solve_7pt_da(float *c, float alpha)
{
    std::cout << "Simulating using 7 point stencil; top boundary: absorbing; bottom boundary: absorbing \n";

    // Checking CFL condition
    float max_c = 0;
    float cur_c;
    for (int i = 0; i < ny; i++)
    {
        cur_c = c[i];
        if (cur_c > max_c)
            max_c = cur_c;
    }
    assert((max_c * dt/dy) < 1 && "CFL condition not satisfied!");

    // Initialise Solution matrix
    MatrixFP32 u = MatrixFP32(ny, nt);

    auto start = std::chrono::high_resolution_clock::now();
    // Loop over time
    for (unsigned long it = 1; it < nt-1; it++)
    {
        // Loop over space
        for (unsigned long iy = 3; iy < ny-3; iy++)
        {
            // 2nd derivative wrt x
            float d2u_dx2 = (2*u(iy-3,it) - 27*u(iy-2,it) + 270*u(iy-1,it) - 490*u(iy,it) + 270*u(iy+1,it) - 27*u(iy+2,it) + 2*u(iy+3,it))/(180*dy*dy);

            // Updating solution
            if (iy == i_src)
                u(iy,it+1) = (c[iy]*dt * c[iy]*dt) * d2u_dx2 + 2*u(iy,it) - u(iy,it-1) + dt*dt * src[it] / dy;
            else
                u(iy,it+1) = (c[iy]*dt * c[iy]*dt) * d2u_dx2 + 2*u(iy,it) - u(iy,it-1);
        }

        // Top: Damping boundary condition
        u(2,it+1) = u(2,it) + (c[2]*dt/dy)*(u(3,it) - u(2,it)) - alpha*(u(2,it) - u(2,it-1));
        u(1,it+1) = u(1,it) + (c[1]*dt/dy)*(u(2,it) - u(1,it)) - alpha*(u(1,it) - u(1,it-1));
        u(0,it+1) = u(0,it) + (c[0]*dt/dy)*(u(1,it) - u(0,it)) - alpha*(u(0,it) - u(0,it-1));

        
        // Bottom: Absorbing boundary condition
        u(ny-3,it+1) = u(ny-4,it) + (c[ny-3]*dt-dy)/(c[ny-3]*dt+dy) * (u(ny-4,it+1)-u(ny-3,it));
        u(ny-2,it+1) = u(ny-3,it) + (c[ny-2]*dt-dy)/(c[ny-2]*dt+dy) * (u(ny-3,it+1)-u(ny-2,it));
        u(ny-1,it+1) = u(ny-2,it) + (c[ny-1]*dt-dy)/(c[ny-1]*dt+dy) * (u(ny-2,it+1)-u(ny-1,it));
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Simulation took: " << duration.count() << " microseconds" << std::endl;

    // Saving simulation details
    // Saving solution as bin
    std::cout << "Saving results and grid details to: data/1d/ \n";
    u.save("data/1d/sol_7pt_da.bin");
    // Saving grid as bin
    save_vec(y, ny, "data/1d/y_7pt_da.bin");
    save_vec(t, nt, "data/1d/t_7pt_da.bin");
    save_vec(c, ny, "data/1d/c_7pt_da.bin");
    std::cout << "Done \n \n";
}

Wave1D::~Wave1D()
{
    delete[] y;
    delete[] t;
    delete[] src;
}