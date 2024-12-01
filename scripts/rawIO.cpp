#include <rawIO.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <string>

void InputRaw(unsigned char* data, const char* inputPath, int width, int height, int bytesPerPixel) 
{
    FILE* file = fopen(inputPath, "rb");
    if (file == nullptr) 
    {
        throw std::runtime_error(std::string("Cannot open file: ") + inputPath + "\nFailed to load image");
    }

    int dataSize = width * height * bytesPerPixel;
    fread(data, sizeof(unsigned char), dataSize, file);
    fclose(file);
}

void OutputRaw(const unsigned char* data, const char* outputPath, int width, int height, int bytesPerPixel) 
{
    FILE* file = fopen(outputPath, "wb");
    if (file == nullptr) 
    {
        throw std::runtime_error(std::string("Cannot open file: ") + outputPath + "\nFailed to load image");
    }

    int dataSize = width * height * bytesPerPixel;
    fwrite(data, sizeof(unsigned char), dataSize, file);
    fclose(file);
}
