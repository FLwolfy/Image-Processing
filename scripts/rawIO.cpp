#include <rawIO.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>

void InputRaw(unsigned char* data, const char* inputPath, int width, int height, int bytesPerPixel) 
{
    FILE* file;
    int dataSize = width * height * bytesPerPixel;

    if (!(file = fopen(inputPath, "rb"))) 
    {
        throw std::runtime_error(std::string("Cannot open file: ") + inputPath + "\nFailed to load image");
    }

    fread(data, sizeof(unsigned char), dataSize, file);
    fclose(file);
}

void OutputRaw(unsigned char* data, const char* outputPath, int width, int height, int bytesPerPixel) 
{
    FILE* file;
    int dataSize = width * height * bytesPerPixel;

    if (!(file = fopen(outputPath, "wb"))) 
    {
        throw std::runtime_error(std::string("Cannot open file: ") + outputPath + "\nFailed to load image");
    }

    fwrite(data, sizeof(unsigned char), dataSize, file);
    fclose(file);
}