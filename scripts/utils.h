#pragma once

#include <kernel.h>

///////////////////////// Math Functions /////////////////////////

unsigned char Max(const unsigned char* array, int length);

unsigned char Min(const unsigned char* array, int length);

unsigned char Maximin(const unsigned char* array, int length);

unsigned char Minimax(const unsigned char* array, int length);

unsigned char Mean(const unsigned char* array, int length);

unsigned char Median(const unsigned char* array, int length, bool pseudo);

float Variance(const unsigned char* array, int length);

float MSE(const unsigned char* data1, const unsigned char* data2, int length);

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

int MaskCount(
    const unsigned char* data,
    int width, int height,
    int bytesPerPixel,
    int channel,
    const Kernel& kernel
);

///////////////////////// K-Means Clustering /////////////////////////

std::vector<int> KMEANSClustering(
    const std::vector<std::vector<float>>& featureMatrix,
    int K,
    int max_iterations
);