#include "../include/Tensor3dFP32.hpp"

// Constructor
Tensor3dFP32::Tensor3dFP32(unsigned long nx, unsigned long ny, unsigned long nz): _nx(nx), _ny(ny), _nz(nz)
{
    // Memory allocation
    _ptr = new float[nx*ny*nz];
    
    // Zero initialization
    for (unsigned long i = 0; i < nx*ny*nz; i++)
        _ptr[i] = 0.0f;
}

// Get value at (i,j)
float Tensor3dFP32::operator()(unsigned long ix, unsigned long iy, unsigned long iz) const
{
    return _ptr[iz*_nx*_ny + iy*_nx + ix];
}

// Set value at (i,j)
float& Tensor3dFP32::operator()(unsigned long ix, unsigned long iy, unsigned long iz)
{
    return _ptr[iz*_nx*_ny + iy*_nx + ix];
}

// Saving ptr data as raw binary
void Tensor3dFP32::save(std::string SAVE_LOC)
{
    std::ofstream outfile_result(SAVE_LOC, std::ios::binary);
    outfile_result.write(reinterpret_cast<char*>(_ptr), _nx*_ny*_nz * sizeof(float));
}

// Destructor
Tensor3dFP32::~Tensor3dFP32()
{
    delete[] _ptr;
}