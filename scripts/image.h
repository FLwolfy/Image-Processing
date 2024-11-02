#pragma once

#include <rawIO.h>
#include <process.h>

#include <stdio.h>
#include <stdlib.h>
#include <vector>

class Image
{
public:
    Image(int m_width, int m_height, int m_bytesPerPixel)
        : m_width(m_width), m_height(m_height), m_bytesPerPixel(m_bytesPerPixel) 
        { m_data = std::vector<unsigned char>(m_width * m_height * m_bytesPerPixel); }; 
    ~Image() = default;

    // IO functions
    inline void Load(const char* inputPath) { InputRaw(m_data.data(), inputPath, m_width, m_height, m_bytesPerPixel); }
    inline void Save(const char* outputPath) { OutputRaw(m_data.data(), outputPath, m_width, m_height, m_bytesPerPixel); }

    // Histogram functions
    std::vector<unsigned int> GetHist() const;
    std::vector<unsigned int> GetCumulativeHist() const;

    // Image processing functions
    static Image GrayScale(const Image& img);
    static Image Negative(const Image& img);
    static Image WaterMark(const Image& img, const Image& watermark, int offsetX, int offsetY, float filterWhiteThreshold, float blendRate);
    static Image LinearScale(const Image& img, int min, int max);
    static Image Enhance(const Image& img, int binSize);

public:
    std::vector<unsigned char> m_data;
    int m_width;
    int m_height;
    int m_bytesPerPixel;
};

