#include <iostream>
#include <vector>         // For using std::vector to store pockets
#include <random>         // For random number generation
#include <cmath>          // For mathematical operations like square roo

#include "./include/Wave3D.hpp"

int main(int argc, char const *argv[])
{
    // Space discretization
    float x_start = 0; float x_end = 10; unsigned long nx = 201;
    float y_start = 0; float y_end = 10; unsigned long ny = 201;
    float z_start = 0; float z_end = 10; unsigned long nz = 201;

    // Time discretization
    float t_start = 0; float t_end = 5; unsigned long nt = 5001;

    // Required resources
    double sizeof_float = sizeof(float);
    std::cout << "Need atleast: " << (4*nx*ny*nz*sizeof_float + nx*nz*nt*sizeof_float + nx*ny*nt*sizeof_float) / 1073741824.0 << " GB \n \n";

    // Medium velocity with randomized pockets
    Tensor3dFP32 c = Tensor3dFP32(nx, ny, nz);
    const float c_air = 5.0f;
    const float c_ground = 10.0f;
    const float c_oil = 3.0f;    // Typical oil velocity
    const float c_gas = 2.5f;    // Typical gas velocity
    const int interface_loc = ny / 4;

    // Random number setup
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    // Generate random pocket parameters
    const int num_oil_pockets = 8;
    const int num_gas_pockets = 5;
    const float max_pocket_size = 0.2f;  // Max 20% of domain size

    struct Pocket {
        float cx, cy, cz;  // Center coordinates
        float rx, ry, rz;  // Radii
    };

    std::vector<Pocket> oil_pockets;
    std::vector<Pocket> gas_pockets;

    // Create oil pockets
    for (int i = 0; i < num_oil_pockets; ++i) {
        oil_pockets.push_back({
            dist(gen) * nx,                  // X-center
            interface_loc + dist(gen) * (ny - interface_loc), // Y-center (ground layer)
            dist(gen) * nz,                  // Z-center
            max_pocket_size * nx * dist(gen), // X-radius
            max_pocket_size * ny * dist(gen), // Y-radius
            max_pocket_size * nz * dist(gen)  // Z-radius
        });
    }

    // Create gas pockets (smaller and fewer)
    for (int i = 0; i < num_gas_pockets; ++i) {
        gas_pockets.push_back({
            dist(gen) * nx,
            interface_loc + dist(gen) * (ny - interface_loc),
            dist(gen) * nz,
            0.5f * max_pocket_size * nx * dist(gen),
            0.5f * max_pocket_size * ny * dist(gen),
            0.5f * max_pocket_size * nz * dist(gen)
        });
    }

    // Fill velocity model
    for (unsigned long iz = 0; iz < nz; iz++) {
        for (unsigned long iy = 0; iy < ny; iy++) {
            for (unsigned long ix = 0; ix < nx; ix++) {
                if (iy < interface_loc) {
                    c(ix, iy, iz) = c_air;
                } else {
                    // Default to ground velocity
                    c(ix, iy, iz) = c_ground;
                    
                    // Check gas pockets first (higher priority)
                    for (const auto& p : gas_pockets) {
                        const float dx = (ix - p.cx)/p.rx;
                        const float dy = (iy - p.cy)/p.ry;
                        const float dz = (iz - p.cz)/p.rz;
                        if (dx*dx + dy*dy + dz*dz <= 1.0f) {
                            c(ix, iy, iz) = c_gas;
                            break;
                        }
                    }
                    
                    // Check oil pockets if not in gas pocket
                    if (c(ix, iy, iz) == c_ground) {
                        for (const auto& p : oil_pockets) {
                            const float dx = (ix - p.cx)/p.rx;
                            const float dy = (iy - p.cy)/p.ry;
                            const float dz = (iz - p.cz)/p.rz;
                            if (dx*dx + dy*dy + dz*dz <= 1.0f) {
                                c(ix, iy, iz) = c_oil;
                                break;
                            }
                        }
                    }
                }
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
