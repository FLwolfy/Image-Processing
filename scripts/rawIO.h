#pragma once

unsigned char* InputRaw(const char* input_file_path, int width, int height, int bytesPerPixel);
void OutputRaw(const char* output_file_path, unsigned char* data, int width, int height, int bytesPerPixel);
