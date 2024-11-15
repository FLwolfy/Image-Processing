#pragma once

#include <vector>
#include <string>

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

    // Regular Kernel
    static Kernel Mean(int size);
    static Kernel Gaussian(float sigma, int size);

    // Edge Detection Kernel
    static Kernel Sobel(int size, bool horizontal);
    static Kernel Laplacian(int size);

    // Morphological Kernel
    static std::vector<Kernel> Pattern(const std::string& type, bool conditional);

public:
    std::vector<float> values;
    int size;

private:
    static std::vector<Kernel> GenerateAllPatterns(const std::vector<std::vector<float>>& patterns, int size);
};
