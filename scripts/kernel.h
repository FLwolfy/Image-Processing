#pragma once

#include <vector>

class Kernel
{
public:
    Kernel(std::vector<float> initValues);
    Kernel(std::initializer_list<float> initValues);
    Kernel(std::vector<float> initValues, int size);
    Kernel(std::initializer_list<float> initValues, int size);

    // Conversion Operators
    operator std::vector<float>() const { return values; }
    float& operator[](size_t index) { return values[index]; }
    const float& operator[](size_t index) const { return values[index]; }

    // Kernel Generating Functions
    static Kernel Sobel(int size, bool horizontal);
    static Kernel Laplacian(int size);
    static Kernel Mean(int size);
    static Kernel Gaussian(float sigma, int size);

public:
    std::vector<float> values;
    int size;
};
