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
    operator std::vector<float>() const { return m_values; }
    float& operator[](size_t index) { return m_values[index]; }
    const float& operator[](size_t index) const { return m_values[index]; }

    // Regular Kernel
    static Kernel Mean(int size);
    static Kernel Gaussian(float sigma, int size);

    // Edge Detection Kernel
    static Kernel Sobel(int size, bool horizontal);
    static Kernel Laplacian(int size);

    // Morphological Kernel
    enum PatternType { SHRINK, THIN, SKELETONIZE, EROSION };
    static std::vector<Kernel> Patterns(PatternType type, bool conditional);

    // Dithering Kernel
    static Kernel BayerIndex(int size);
    static Kernel BayerThreshold(int size);
    static Kernel FloydSteinberg();

    // Texture analysis Kernel
    static std::vector<Kernel> LawsFilters(int size);

    // Bit Quads Kernel
    static std::vector<Kernel> BitQuads(char pattern);

public:
    std::vector<float> m_values;
    int m_size;

private:
    static std::vector<Kernel> GenerateAllPatterns(const std::vector<std::vector<float>>& patterns, int size);
};
