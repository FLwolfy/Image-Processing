#include <rawIO.h>

#include <stdio.h>
#include <stdlib.h>

unsigned char* InputRaw(const char* input_file_path, int width, int height, int bytesPerPixel) 
{
    FILE* file;
    int dataSize = width * height * bytesPerPixel;

    unsigned char* data = (unsigned char*)malloc(dataSize);
    if (!data) 
    {
        printf("\n Memory allocation failed\n");
        exit(1);
    }

    if (!(file = fopen(input_file_path, "rb"))) 
    {
        printf("\n Cannot open file: %s\n", input_file_path);
        free(data);
        exit(1);
    }

    fread(data, sizeof(unsigned char), dataSize, file);
    fclose(file);

    return data;
}

void OutputRaw(const char* output_file_path, unsigned char* data, int width, int height, int bytesPerPixel) 
{
    FILE* file;
    int dataSize = width * height * bytesPerPixel;

    if (!(file = fopen(output_file_path, "wb"))) 
    {
        printf("\n Cannot open file: %s\n", output_file_path);
        exit(1);
    }

    fwrite(data, sizeof(unsigned char), dataSize, file);
    fclose(file);
}