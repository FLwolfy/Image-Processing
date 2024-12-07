#pragma once

#include <kernel.h>

///////////////////////// Math Functions /////////////////////////

unsigned char Max(const unsigned char* array, int length);

unsigned char Min(const unsigned char* array, int length);

unsigned char Maximin(const unsigned char* array, int length);

unsigned char Minimax(const unsigned char* array, int length);

unsigned char Median(const unsigned char* array, int length, bool pseudo);

///////////////////////// Array Functions /////////////////////////

std::vector<unsigned char> Sort(
    const unsigned char* array,
    int length
);

std::vector<unsigned char> Hash(
    const unsigned char* data,
    int width, int height,
    int bytesPerPixel,
    int channel,
    unsigned long long seed
);

std::vector<unsigned char> Convolve(
    const unsigned char* data,
    int width, int height,
    int bytesPerPixel,
    int channel,
    const Kernel& kernel
);

std::vector<float> ConvolvePrecise(
    const unsigned char* data,
    int width, int height,
    int bytesPerPixel,
    int channel,
    const Kernel& kernel
);

std::vector<bool> Mask(
    const unsigned char* data,
    int width, int height,
    int bytesPerPixel,
    int channel,
    const std::vector<Kernel>& kernels
);

std::vector<bool> MaskBool(
    const std::vector<bool>& data,
    int width, int height,
    const std::vector<Kernel>& kernels
);

std::vector<int> kMeansClustering(
    const std::vector<std::vector<float>>& data,
    int numOfClusters
);