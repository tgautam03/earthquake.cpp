#ifndef TENSOR3DFP32
#define TENSOR3DFP32

#include <string>   // Data saving location
#include <fstream>  // Data saving
#include <iostream> // Printing

class Tensor3dFP32
{
private:
    // Data Pointer
    float *_ptr;
public:
    // Matrix Dimensions
    const unsigned long _nx;
    const unsigned long _ny;
    const unsigned long _nz;
    
    // Initialize dynamic array to zeros
    Tensor3dFP32(unsigned long nx, unsigned long ny, unsigned long nz);

    // Get value at ptr(i,j)
    float operator()(unsigned long ix, unsigned long iy, unsigned long iz) const;

    // Set value at ptr(i,j)
    float& operator()(unsigned long ix, unsigned long iy, unsigned long iz);

    // Save ptr data as raw binary
    void save(std::string SAVE_LOC);

    // Destructor
    ~Tensor3dFP32();
};

#endif