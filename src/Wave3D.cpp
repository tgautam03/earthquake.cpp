#include "../include/Wave3D.hpp"

Wave3D::Wave3D(float xs_, float xe_, unsigned long nx_,
    float ys_, float ye_, unsigned long ny_,
    float zs_, float ze_, unsigned long nz_,
    float ts_, float te_, float nt_): xs(xs_), xe(xe_), nx(nx_),
                ys(ys_), ye(ye_), ny(ny_),
                zs(zs_), ze(ze_), nz(nz_),
                ts(ts_), te(te_), nt(nt_),
                dx((xe - xs) / static_cast<float>(nx-1)),
                dy((ye - ys) / static_cast<float>(ny-1)),
                dz((ze - zs) / static_cast<float>(nz-1)),
                dt((te - ts) / static_cast<float>(nt-1))
{
    // Defining x,y,z and t grid
    x = new float[nx];
    linspace_1d(x, nx, xs, dx);

    y = new float[ny];
    linspace_1d(y, ny, ys, dy);

    z = new float[nz];
    linspace_1d(z, nz, zs, dz);

    t = new float[nt];
    linspace_1d(t, nt, ts, dt);
}

void Wave3D::add_src_dGaussian(float t0, float f0, float src_xloc, float src_yloc, float src_zloc)
{
    // Source location on the grid point
    ix_src = static_cast<unsigned long>(src_xloc / dx) + 1;
    iy_src = static_cast<unsigned long>(src_yloc / dy) + 1;
    iz_src = static_cast<unsigned long>(src_zloc / dz) + 1;

    // Defining source function
    src = new float[nt];
    for (int i = 0; i < nt; i++)
    {
        float t_i = t[i];
        src[i] = -8 * (t_i - t0) * f0 * exp(-1 * 4*f0*(t_i - t0) * 4*f0*(t_i - t0));
    }
}

void Wave3D::solve(Tensor3dFP32 &c, float alpha)
{
    std::cout << "Simulating using 7 point stencil \n";

    // Checking CFL condition
    float max_c = 0;
    float cur_c;
    for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
            for (int iz = 0; iz < nz; iz++) {
                cur_c = c(ix,iy,iz);
                if (cur_c > max_c)
                    max_c = cur_c;
            }
        }
    }
    assert((max_c * dt * sqrt((1/(dx*dx)) + (1/(dy*dy)) + (1/(dz*dz)))) < 1 && "CFL condition not satisfied!");

    // Initialise Solution matrix
    Tensor4dFP32 u = Tensor4dFP32(nx, ny, nz, nt);

    // Loop over time
    auto start = std::chrono::high_resolution_clock::now();
    for (unsigned long it = 1; it < nt-1; it++) {
        // Loop over space
        for (unsigned long iz = 3; iz < nz-3; iz++) {
            for (unsigned long iy = 3; iy < ny-3; iy++) {
                for (unsigned long ix = 3; ix < nx-3; ix++) {
                    // 2nd derivative wrt x
                    float d2u_dx2 = (2*u(ix-3,iy,iz,it) - 27*u(ix-2,iy,iz,it) + 270*u(ix-1,iy,iz,it) - 490*u(ix,iy,iz,it) + 270*u(ix+1,iy,iz,it) - 27*u(ix+2,iy,iz,it) + 2*u(ix+3,iy,iz,it))/(180*dx*dx);
                    float d2u_dy2 = (2*u(ix,iy-3,iz,it) - 27*u(ix,iy-2,iz,it) + 270*u(ix,iy-1,iz,it) - 490*u(ix,iy,iz,it) + 270*u(ix,iy+1,iz,it) - 27*u(ix,iy+2,iz,it) + 2*u(ix,iy+3,iz,it))/(180*dy*dy);
                    float d2u_dz2 = (2*u(ix,iy,iz-3,it) - 27*u(ix,iy,iz-2,it) + 270*u(ix,iy,iz-1,it) - 490*u(ix,iy,iz,it) + 270*u(ix,iy,iz+1,it) - 27*u(ix,iy,iz+2,it) + 2*u(ix,iy,iz+3,it))/(180*dz*dz);

                    // Updating solution
                    if ((ix == ix_src) && (iy == iy_src) && (iz == iz_src))
                        u(ix,iy,iz,it+1) = (c(ix,iy,iz)*dt * c(ix,iy,iz)*dt) * (d2u_dx2 + d2u_dy2 + d2u_dz2) + 2*u(ix,iy,iz,it) - u(ix,iy,iz,it-1) + dt*dt * src[it] / (dx*dy*dz);
                    else
                        u(ix,iy,iz,it+1) = (c(ix,iy,iz)*dt * c(ix,iy,iz)*dt) * (d2u_dx2 + d2u_dy2 + d2u_dz2) + 2*u(ix,iy,iz,it) - u(ix,iy,iz,it-1);
                }
            }
        }

        // Boundaries
        for (unsigned long iy = 0; iy < nx; iy++) {
            for (unsigned long ix = 0; ix < ny; ix++) {
                // Front face: Absorbing
                u(ix,iy,2,it+1) = u(ix,iy,3,it) + (c(ix,iy,2)*dt-dy)/(c(ix,iy,2)*dt+dy) * (u(ix,iy,3,it+1)-u(ix,iy,2,it));
                u(ix,iy,1,it+1) = u(ix,iy,2,it) + (c(ix,iy,1)*dt-dy)/(c(ix,iy,1)*dt+dy) * (u(ix,iy,2,it+1)-u(ix,iy,1,it));
                u(ix,iy,0,it+1) = u(ix,iy,1,it) + (c(ix,iy,0)*dt-dy)/(c(ix,iy,0)*dt+dy) * (u(ix,iy,1,it+1)-u(ix,iy,0,it));

                // Back face: Absorbing
                u(ix,iy,nz-3,it+1) = u(ix,iy,nz-4,it) + (c(ix,iy,nz-3)*dt-dy)/(c(ix,iy,nz-3)*dt+dy) * (u(ix,iy,nz-4,it+1)-u(ix,iy,nz-3,it));
                u(ix,iy,nz-2,it+1) = u(ix,iy,nz-3,it) + (c(ix,iy,nz-2)*dt-dy)/(c(ix,iy,nz-2)*dt+dy) * (u(ix,iy,nz-3,it+1)-u(ix,iy,nz-2,it));
                u(ix,iy,nz-1,it+1) = u(ix,iy,nz-2,it) + (c(ix,iy,nz-1)*dt-dy)/(c(ix,iy,nz-1)*dt+dy) * (u(ix,iy,nz-2,it+1)-u(ix,iy,nz-1,it));
            }
        }

        for (unsigned long iz = 0; iz < nz; iz++) {
            for (unsigned long iy = 0; iy < ny; iy++) {
                // Side face: Absorbing
                u(2,iy,iz,it+1) = u(3,iy,iz,it) + (c(2,iy,iz)*dt-dy)/(c(2,iy,iz)*dt+dy) * (u(3,iy,iz,it+1)-u(2,iy,iz,it));
                u(1,iy,iz,it+1) = u(2,iy,iz,it) + (c(1,iy,iz)*dt-dy)/(c(1,iy,iz)*dt+dy) * (u(2,iy,iz,it+1)-u(1,iy,iz,it));
                u(0,iy,iz,it+1) = u(1,iy,iz,it) + (c(0,iy,iz)*dt-dy)/(c(0,iy,iz)*dt+dy) * (u(1,iy,iz,it+1)-u(0,iy,iz,it));

                // Side face: Absorbing
                u(nx-3,iy,iz,it+1) = u(nx-4,iy,iz,it) + (c(nx-3,iy,iz)*dt-dy)/(c(nx-3,iy,iz)*dt+dy) * (u(nx-4,iy,iz,it+1)-u(nx-3,iy,iz,it));
                u(nx-2,iy,iz,it+1) = u(nx-3,iy,iz,it) + (c(nx-2,iy,iz)*dt-dy)/(c(nx-2,iy,iz)*dt+dy) * (u(nx-3,iy,iz,it+1)-u(nx-2,iy,iz,it));
                u(nx-1,iy,iz,it+1) = u(nx-2,iy,iz,it) + (c(nx-1,iy,iz)*dt-dy)/(c(nx-1,iy,iz)*dt+dy) * (u(nx-2,iy,iz,it+1)-u(nx-1,iy,iz,it));
            }
        }

        for (unsigned long iz = 0; iz < nz; iz++) {
            for (unsigned long ix = 0; ix < nx; ix++) {
                // Top face: Damping
                u(ix,2,iz,it+1) = u(ix,2,iz,it) + (c(ix,2,iz)*dt/dy)*(u(ix,3,iz,it) - u(ix,2,iz,it)) - alpha*(u(ix,2,iz,it) - u(ix,2,iz,it-1));
                u(ix,1,iz,it+1) = u(ix,1,iz,it) + (c(ix,1,iz)*dt/dy)*(u(ix,2,iz,it) - u(ix,1,iz,it)) - alpha*(u(ix,1,iz,it) - u(ix,1,iz,it-1));
                u(ix,0,iz,it+1) = u(ix,0,iz,it) + (c(ix,0,iz)*dt/dy)*(u(ix,1,iz,it) - u(ix,0,iz,it)) - alpha*(u(ix,0,iz,it) - u(ix,0,iz,it-1));


                // Bottom face: Absorbing
                u(ix,ny-3,iz,it+1) = u(ix,ny-4,iz,it) + (c(ix,ny-3,iz)*dt-dy)/(c(ix,ny-3,iz)*dt+dy) * (u(ix,ny-4,iz,it+1)-u(ix,ny-3,iz,it));
                u(ix,ny-2,iz,it+1) = u(ix,ny-3,iz,it) + (c(ix,ny-2,iz)*dt-dy)/(c(ix,ny-2,iz)*dt+dy) * (u(ix,ny-3,iz,it+1)-u(ix,ny-2,iz,it));
                u(ix,ny-1,iz,it+1) = u(ix,ny-2,iz,it) + (c(ix,ny-1,iz)*dt-dy)/(c(ix,ny-1,iz)*dt+dy) * (u(ix,ny-2,iz,it+1)-u(ix,ny-1,iz,it));
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "Simulation took: " << duration.count() << " seconds" << std::endl;

    // Saving simulation details
    std::cout << "Saving results and grid details to: data/3d/ \n";
    // Saving solution as bin
    u.save("data/3d/sol.bin");
    // Saving grid as bin
    save_vec(x, nx, "data/3d/x.bin");
    save_vec(y, ny, "data/3d/y.bin");
    save_vec(z, nz, "data/3d/z.bin");
    save_vec(t, nt, "data/3d/t.bin");
    c.save("data/3d/c.bin");
    std::cout << "Done \n";
}

void Wave3D::solve_slices(Tensor3dFP32 &c, float alpha, float where_z)
{
    std::cout << "Simulating using 7 point stencil \n";

    // Checking CFL condition
    float max_c = 0;
    float cur_c;
    for (int iy = 0; iy < ny; iy++)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            for (int iz = 0; iz < nz; iz++)
            {
                cur_c = c(ix,iy,iz);
                if (cur_c > max_c)
                    max_c = cur_c;
            }
        }
    }
    assert((max_c * dt * sqrt((1/(dx*dx)) + (1/(dy*dy)) + (1/(dz*dz)))) < 1 && "CFL condition not satisfied!");

    // Initialise slices
    Tensor3dFP32 u_surface = Tensor3dFP32(nx, nz, nt);
    
    int z_loc = static_cast<unsigned long>(where_z / dz) + 1;
    Tensor3dFP32 u_xy = Tensor3dFP32(nx, ny, nt);

    // Initialise Solution matrix
    Tensor3dFP32 u_prev = Tensor3dFP32(nx, ny, nz);
    Tensor3dFP32 u_now = Tensor3dFP32(nx, ny, nz);
    Tensor3dFP32 u_next = Tensor3dFP32(nx, ny, nz);


    // Loop over time
    auto start = std::chrono::high_resolution_clock::now();
    for (unsigned long it = 1; it < nt-1; it++)
    {
        // Loop over space
        for (unsigned long iz = 3; iz < nz-3; iz++)
        {
            for (unsigned long iy = 3; iy < ny-3; iy++)
            {
                for (unsigned long ix = 3; ix < nx-3; ix++)
                {
                    // 2nd derivative wrt x
                    float d2u_dx2 = (2*u_now(ix-3,iy,iz) - 27*u_now(ix-2,iy,iz) + 270*u_now(ix-1,iy,iz) - 490*u_now(ix,iy,iz) + 270*u_now(ix+1,iy,iz) - 27*u_now(ix+2,iy,iz) + 2*u_now(ix+3,iy,iz))/(180*dx*dx);
                    float d2u_dy2 = (2*u_now(ix,iy-3,iz) - 27*u_now(ix,iy-2,iz) + 270*u_now(ix,iy-1,iz) - 490*u_now(ix,iy,iz) + 270*u_now(ix,iy+1,iz) - 27*u_now(ix,iy+2,iz) + 2*u_now(ix,iy+3,iz))/(180*dy*dy);
                    float d2u_dz2 = (2*u_now(ix,iy,iz-3) - 27*u_now(ix,iy,iz-2) + 270*u_now(ix,iy,iz-1) - 490*u_now(ix,iy,iz) + 270*u_now(ix,iy,iz+1) - 27*u_now(ix,iy,iz+2) + 2*u_now(ix,iy,iz+3))/(180*dz*dz);

                    // Updating solution
                    if ((ix == ix_src) && (iy == iy_src) && (iz == iz_src))
                        u_next(ix,iy,iz) = (c(ix,iy,iz)*dt * c(ix,iy,iz)*dt) * (d2u_dx2 + d2u_dy2 + d2u_dz2) + 2*u_now(ix,iy,iz) - u_prev(ix,iy,iz) + dt*dt * src[it] / (dx*dy*dz);
                    else
                        u_next(ix,iy,iz) = (c(ix,iy,iz)*dt * c(ix,iy,iz)*dt) * (d2u_dx2 + d2u_dy2 + d2u_dz2) + 2*u_now(ix,iy,iz) - u_prev(ix,iy,iz);
                }
            }
        }
        // ==================================================================== //
        // ============================== Boundaries ========================== //
        // ==================================================================== //
        for (unsigned long iy = 0; iy < nx; iy++)
        {
            for (unsigned long ix = 0; ix < ny; ix++)
            {
                // Front face: Absorbing
                u_next(ix,iy,2) = u_now(ix,iy,3) + (c(ix,iy,2)*dt-dy)/(c(ix,iy,2)*dt+dy) * (u_next(ix,iy,3)-u_now(ix,iy,2));
                u_next(ix,iy,1) = u_now(ix,iy,2) + (c(ix,iy,1)*dt-dy)/(c(ix,iy,1)*dt+dy) * (u_next(ix,iy,2)-u_now(ix,iy,1));
                u_next(ix,iy,0) = u_now(ix,iy,1) + (c(ix,iy,0)*dt-dy)/(c(ix,iy,0)*dt+dy) * (u_next(ix,iy,1)-u_now(ix,iy,0));

                // Back face: Absorbing
                u_next(ix,iy,nz-3) = u_now(ix,iy,nz-4) + (c(ix,iy,nz-3)*dt-dy)/(c(ix,iy,nz-3)*dt+dy) * (u_next(ix,iy,nz-4)-u_now(ix,iy,nz-3));
                u_next(ix,iy,nz-2) = u_now(ix,iy,nz-3) + (c(ix,iy,nz-2)*dt-dy)/(c(ix,iy,nz-2)*dt+dy) * (u_next(ix,iy,nz-3)-u_now(ix,iy,nz-2));
                u_next(ix,iy,nz-1) = u_now(ix,iy,nz-2) + (c(ix,iy,nz-1)*dt-dy)/(c(ix,iy,nz-1)*dt+dy) * (u_next(ix,iy,nz-2)-u_now(ix,iy,nz-1));
            }
        }

        for (unsigned long iz = 0; iz < nz; iz++)
        {
            for (unsigned long iy = 0; iy < ny; iy++)
            {
                // Side face: Absorbing
                u_next(2,iy,iz) = u_now(3,iy,iz) + (c(2,iy,iz)*dt-dy)/(c(2,iy,iz)*dt+dy) * (u_next(3,iy,iz)-u_now(2,iy,iz));
                u_next(1,iy,iz) = u_now(2,iy,iz) + (c(1,iy,iz)*dt-dy)/(c(1,iy,iz)*dt+dy) * (u_next(2,iy,iz)-u_now(1,iy,iz));
                u_next(0,iy,iz) = u_now(1,iy,iz) + (c(0,iy,iz)*dt-dy)/(c(0,iy,iz)*dt+dy) * (u_next(1,iy,iz)-u_now(0,iy,iz));

                // Side face: Absorbing
                u_next(nx-3,iy,iz) = u_now(nx-4,iy,iz) + (c(nx-3,iy,iz)*dt-dy)/(c(nx-3,iy,iz)*dt+dy) * (u_next(nx-4,iy,iz)-u_now(nx-3,iy,iz));
                u_next(nx-2,iy,iz) = u_now(nx-3,iy,iz) + (c(nx-2,iy,iz)*dt-dy)/(c(nx-2,iy,iz)*dt+dy) * (u_next(nx-3,iy,iz)-u_now(nx-2,iy,iz));
                u_next(nx-1,iy,iz) = u_now(nx-2,iy,iz) + (c(nx-1,iy,iz)*dt-dy)/(c(nx-1,iy,iz)*dt+dy) * (u_next(nx-2,iy,iz)-u_now(nx-1,iy,iz));
            }
        }

        for (unsigned long iz = 0; iz < nz; iz++)
        {
            for (unsigned long ix = 0; ix < nx; ix++)
            {
                // Top face: Damping
                u_next(ix,2,iz) = u_now(ix,2,iz) + (c(ix,2,iz)*dt/dy)*(u_now(ix,3,iz) - u_now(ix,2,iz)) - alpha*(u_now(ix,2,iz) - u_prev(ix,2,iz));
                u_next(ix,1,iz) = u_now(ix,1,iz) + (c(ix,1,iz)*dt/dy)*(u_now(ix,2,iz) - u_now(ix,1,iz)) - alpha*(u_now(ix,1,iz) - u_prev(ix,1,iz));
                u_next(ix,0,iz) = u_now(ix,0,iz) + (c(ix,0,iz)*dt/dy)*(u_now(ix,1,iz) - u_now(ix,0,iz)) - alpha*(u_now(ix,0,iz) - u_prev(ix,0,iz));


                // Bottom face: Absorbing
                u_next(ix,ny-3,iz) = u_now(ix,ny-4,iz) + (c(ix,ny-3,iz)*dt-dy)/(c(ix,ny-3,iz)*dt+dy) * (u_next(ix,ny-4,iz)-u_now(ix,ny-3,iz));
                u_next(ix,ny-2,iz) = u_now(ix,ny-3,iz) + (c(ix,ny-2,iz)*dt-dy)/(c(ix,ny-2,iz)*dt+dy) * (u_next(ix,ny-3,iz)-u_now(ix,ny-2,iz));
                u_next(ix,ny-1,iz) = u_now(ix,ny-2,iz) + (c(ix,ny-1,iz)*dt-dy)/(c(ix,ny-1,iz)*dt+dy) * (u_next(ix,ny-2,iz)-u_now(ix,ny-1,iz));
            }
        }

        // ==================================================================== //
        // ===================== Now -> Prev; Future -> Now =================== //
        // ==================================================================== //
        for (unsigned long iz = 0; iz < nz; iz++)
        {
            for (unsigned long iy = 0; iy < ny; iy++)
            {
                for (unsigned long ix = 0; ix < nx; ix++)
                {
                    u_prev(ix, iy, iz) = u_now(ix, iy, iz);
                    u_now(ix, iy, iz) = u_next(ix, iy, iz);
                }
            }
        }

        // ==================================================================== //
        // ============================ Saving Slices ========================= //
        // ==================================================================== //
        // Vertical slice
        for (unsigned long iy = 0; iy < nx; iy++)
        {
            for (unsigned long ix = 0; ix < ny; ix++)
                u_xy(ix, iy, it+1) = u_now(ix, iy, z_loc);
        }

        // Surface slice
        for (unsigned long iz = 0; iz < nz; iz++)
        {
            for (unsigned long ix = 0; ix < nx; ix++)
                u_surface(ix, iz, it+1) = u_now(ix, 0, iz);
        }


    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "Simulation took: " << duration.count() << " seconds" << std::endl;

    // Saving simulation details
    std::cout << "Saving results and grid details to: data/3d/ \n";
    // Saving solution as bin
    u_xy.save("data/3d/sol_vert.bin");
    u_surface.save("data/3d/sol_surf.bin");
    // Saving grid as bin
    save_vec(x, nx, "data/3d/x.bin");
    save_vec(y, ny, "data/3d/y.bin");
    save_vec(z, nz, "data/3d/z.bin");
    save_vec(t, nt, "data/3d/t.bin");
    c.save("data/3d/c.bin");
    std::cout << "Done \n";
}

Wave3D::~Wave3D()
{
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] t;
    delete[] src;
}