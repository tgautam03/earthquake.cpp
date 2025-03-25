#include "../include/utils.hpp"

// [n_points] points from [start] that are [grid spacing] apart
void linspace_1d(float *arr, int n_points, float start, float grid_spacing)
{
    for (int i = 0; i < n_points; i++)
        arr[i] = start + grid_spacing*i;
}

// Printing 1D array [arr]
void print_vec(float *arr, int n_points)
{
    for (int i = 0; i < n_points; i++)
        std::cout << FIXED_FLOAT(arr[i]) << " ";
    std::cout << "\n";
}

// Saving ptr data as raw binary
void save_vec(float *ptr, int len, std::string DIR_LOC)
{
    std::ofstream outfile_result(DIR_LOC, std::ios::binary);
    outfile_result.write(reinterpret_cast<char*>(ptr), len * sizeof(float));
}