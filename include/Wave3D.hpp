#ifndef WAVE3D
#define WAVE3D

#include <cmath>
#include <cassert>
#include <chrono>

#include "Tensor4dFP32.hpp"
#include "Tensor3dFP32.hpp"
#include "utils.hpp"

class Wave3D
{
private:
    // Space discretization
    float *x; const float xs; const float xe; unsigned long nx; const float dx;
    float *y; const float ys; const float ye; unsigned long ny; const float dy;
    float *z; const float zs; const float ze; unsigned long nz; const float dz;

    // Time discretization
    float *t; const float ts; const float te; unsigned long nt; const float dt;

    // Source function
    float *src; unsigned long ix_src; unsigned long iy_src; unsigned long iz_src;

public:
    // Default constructor
    Wave3D(float xs_, float xe_, unsigned long nx_, 
            float ys_, float ye_, unsigned long ny_, 
            float zs_, float ze_, unsigned long nz_, 
            float ts_, float te_, float nt_);

    // Add source
    void add_src_dGaussian(float t0, float f0, float src_xloc, float src_yloc, float src_zloc); // Derivative of Gaussian source

    // Solve
    void solve(Tensor3dFP32 &c, float alpha);
    void solve_slices(Tensor3dFP32 &c, float alpha, float where_z);

    ~Wave3D();
};

#endif