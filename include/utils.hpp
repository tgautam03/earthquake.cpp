#ifndef UTILS
#define UTILS

#include <iostream> // For printing
#include <iomanip> // Setting decimal places for printing
#include <fstream> // Saving data as bin

// Macro to print 2 decimal places
#define FIXED_FLOAT(val) std::fixed << std::setprecision(2) << (val)

void linspace_1d(float *arr, int n_points, float start, float grid_spacing);

void print_vec(float *arr, int n_points);

// Saving ptr data as raw binary
void save_vec(float *ptr, int len, std::string DIR_LOC);

#endif