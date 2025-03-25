#ifndef TENSOR4DFP32
#define TENSOR4DFP32

#include <string>   // Data saving location
#include <fstream>  // Data saving
#include <iostream> // Printing

class Tensor4dFP32
{
private:
    // Data Pointer
    float *_ptr;
public:
    // Matrix Dimensions
    const unsigned long _nx;
    const unsigned long _ny;
    const unsigned long _nz;
    const unsigned long _nt;
    
    // Initialize dynamic array to zeros
    Tensor4dFP32(unsigned long nx, unsigned long ny, unsigned long nz, unsigned long nt);

    // Get value at ptr(i,j)
    float operator()(unsigned long ix, unsigned long iy, unsigned long iz, unsigned long it) const;

    // Set value at ptr(i,j)
    float& operator()(unsigned long ix, unsigned long iy, unsigned long iz, unsigned long it);

    // Save ptr data as raw binary
    void save(std::string SAVE_LOC);

    // Destructor
    ~Tensor4dFP32();
};

#endif