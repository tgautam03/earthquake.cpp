#ifndef MATRIXFP32
#define MATRIXFP32

#include <string>   // Data saving location
#include <fstream>  // Data saving
#include <iostream> // Printing

class MatrixFP32
{
private:
    // Data Pointer
    float *_ptr;
public:
    // Matrix Dimensions
    const unsigned long _ny;
    const unsigned long _nt;
    
    // Initialize dynamic array to zeros
    MatrixFP32(unsigned long ny, unsigned long nt);

    // Get value at ptr(i,j)
    float operator()(unsigned long iy, unsigned long it) const;

    // Set value at ptr(i,j)
    float& operator()(unsigned long iy, unsigned long it);

    // Save ptr data as raw binary
    void save(std::string SAVE_LOC);

    // Destructor
    ~MatrixFP32();
};

#endif