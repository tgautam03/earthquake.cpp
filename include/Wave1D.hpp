#ifndef WAVE1D
#define WAVE1D

#include <cmath>
#include <chrono>
#include <cassert>

#include "MatrixFP32.hpp"
#include "utils.hpp"

class Wave1D
{
private:
    // Space discretization
    float *y; const float ys; const float ye; unsigned long ny; const float dy;

    // Time discretization
    float *t; const float ts; const float te; unsigned long nt; const float dt;

    // Source function
    float *src; unsigned long i_src;

public:
    // Default constructor
    Wave1D(float ys_, float ye_, unsigned long ny_, float ts_, float te_, float nt_);

    // Add source
    void add_src_dGaussian(float t0, float f0, float src_loc); // Derivative of Gaussian source

    // Solve
    void solve_3pt_00(float *c);    // 3 point; Top boundary: zero; Bottom boundary: zero
    void solve_3pt_da(float *c, float alpha);    // 3 point; Top boundary: damping; Bottom boundary: absorbing
    void solve_5pt_da(float *c, float alpha);    // 5 point; Top boundary: damping; Bottom boundary: absorbing
    void solve_7pt_da(float *c, float alpha);    // 7 point; Top boundary: damping; Bottom boundary: absorbing

    ~Wave1D();
};

#endif