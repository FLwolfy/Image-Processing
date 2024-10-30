#pragma once

#include <stdio.h>
#include <stdlib.h>

#include <rawIO.h>
#include <process.h>

class Image
{
public:
    Image(int m_width, int m_height, int m_bytesPerPixel)
        : m_width(m_width), m_height(m_height), m_bytesPerPixel(m_bytesPerPixel) {};
    ~Image() { free(m_data); };

    // IO functions
    inline void Load(const char* input_file_path) { free(m_data); m_data = InputRaw(input_file_path, m_width, m_height, m_bytesPerPixel); }
    inline void Save(const char* output_file_path) { OutputRaw(output_file_path, m_data, m_width, m_height, m_bytesPerPixel); };

    // Image processing functions
    static Image& GrayScale(Image& img);
    static Image& WaterMark(Image& img, Image& watermark, int offsetX, int offsetY, float filterWhiteThreshold, float blendRate);

private:
    inline void SetData(unsigned char* data) { m_data = data; }
    
public:
    unsigned char* m_data = nullptr;
    int m_width;
    int m_height;
    int m_bytesPerPixel;
};

