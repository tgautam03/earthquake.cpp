#include "../include/MatrixFP32.hpp"

// Constructor
MatrixFP32::MatrixFP32(unsigned long ny, unsigned long nt): _ny(ny), _nt(nt)
{
    // Memory allocation
    _ptr = new float[ny*nt];
    
    // Zero initialization
    for (unsigned long i = 0; i < ny*nt; i++)
        _ptr[i] = 0.0f;
}

// Get value at (i,j)
float MatrixFP32::operator()(unsigned long iy, unsigned long it) const
{
    return _ptr[it*_ny + iy];
}

// Set value at (i,j)
float& MatrixFP32::operator()(unsigned long iy, unsigned long it)
{
    return _ptr[it*_ny + iy];
}

// Saving ptr data as raw binary
void MatrixFP32::save(std::string SAVE_LOC)
{
    std::ofstream outfile_result(SAVE_LOC, std::ios::binary);
    outfile_result.write(reinterpret_cast<char*>(_ptr), _ny*_nt * sizeof(float));
}

// Destructor
MatrixFP32::~MatrixFP32()
{
    delete[] _ptr;
}