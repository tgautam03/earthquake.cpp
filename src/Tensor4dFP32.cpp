#include "../include/Tensor4dFP32.hpp"

// Constructor
Tensor4dFP32::Tensor4dFP32(unsigned long nx, unsigned long ny, unsigned long nz, unsigned long nt): _nx(nx), _ny(ny), _nz(nz), _nt(nt)
{
    // Memory allocation
    _ptr = new float[nx*ny*nz*nt];
    
    // Zero initialization
    for (unsigned long i = 0; i < nx*ny*nz*nt; i++)
        _ptr[i] = 0.0f;
}

// Get value at (i,j)
float Tensor4dFP32::operator()(unsigned long ix, unsigned long iy, unsigned long iz, unsigned long it) const
{
    return _ptr[it*_nx*_ny*_nz + iz*_nx*_ny + iy*_nx + ix];
}

// Set value at (i,j)
float& Tensor4dFP32::operator()(unsigned long ix, unsigned long iy, unsigned long iz, unsigned long it)
{
    return _ptr[it*_nx*_ny*_nz + iz*_nx*_ny + iy*_nx + ix];
}

// Saving ptr data as raw binary
void Tensor4dFP32::save(std::string SAVE_LOC)
{
    std::ofstream outfile_result(SAVE_LOC, std::ios::binary);
    outfile_result.write(reinterpret_cast<char*>(_ptr), _nx*_ny*_nz*_nt * sizeof(float));
}

// Destructor
Tensor4dFP32::~Tensor4dFP32()
{
    delete[] _ptr;
}