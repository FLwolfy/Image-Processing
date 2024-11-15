#include <rawIO.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <string>

void InputRaw(unsigned char* data, const char* inputPath, int width, int height, int bytesPerPixel) 
{
    FILE* file = nullptr;
    int dataSize = width * height * bytesPerPixel;

    errno_t err = fopen_s(&file, inputPath, "rb");
    if (err != 0) 
    {
        throw std::runtime_error(std::string("Cannot open file: ") + inputPath + "\nFailed to load image");
    }

    fread(data, sizeof(unsigned char), dataSize, file);
    fclose(file);
}

void OutputRaw(const unsigned char* data, const char* outputPath, int width, int height, int bytesPerPixel) 
{
    FILE* file = nullptr;
    int dataSize = width * height * bytesPerPixel;

    errno_t err = fopen_s(&file, outputPath, "wb");
    if (err != 0) 
    {
        throw std::runtime_error(std::string("Cannot open file: ") + outputPath + "\nFailed to load image");
    }

    fwrite(data, sizeof(unsigned char), dataSize, file);
    fclose(file);
}
