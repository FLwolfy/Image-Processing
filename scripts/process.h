#pragma once

#include <vector>
#include <string>
#include <unordered_map>

///////////////////////// Regular processing functions /////////////////////////

std::vector<unsigned char> SeparateChannel(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel, 
    int channel
);

std::vector<unsigned char> ToGrayScale(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel
);

std::vector<unsigned char> AddWatermark(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel, 
    const unsigned char* watermark,
    int watermarkWidth, 
    int watermarkHeight, 
    int offsetX, 
    int offsetY,
    unsigned char threshold,
    float blendRate
);

std::vector<unsigned char> ToNegative(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel
);

///////////////////////// Enhancement functions /////////////////////////

std::vector<unsigned char> ToLinearScale(
    const unsigned char* data, 
    int width, int height,
    int bytesPerPixel,
    int channel,
    int min, int max
);

std::vector<unsigned char> EqualizeHistogram(
    const unsigned char* data,
    const unsigned int* cumulativeHist,
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int binSize
);

///////////////////////// Noise Removal functions /////////////////////////

std::vector<unsigned char> MeanFilter(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int windowSize
);

std::vector<unsigned char> MedianFilter(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int windowSize,
    bool pseudo
);

std::vector<unsigned char> GaussianFilter(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int windowSize,
    float STD
);

std::vector<unsigned char> BilateralFilter(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int windowSize,
    float spaceSTD,
    float colorSTD
);

///////////////////////// Edge Detection functions /////////////////////////

std::vector<unsigned char> ToSobelEdge(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int windowSize,
    const std::string& suppressedMethod,
    const std::string& thresholdMethod,
    const std::unordered_map<std::string, float>& thresholds
);

std::vector<unsigned char> ToLaplacianEdge(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int windowSize,
    float threshold
);

///////////////////////// Morphological functions /////////////////////////

std::vector<unsigned char> Morpho(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    const std::string& type,
    int iterations
);

std::vector<unsigned char> OpenClose(
    const unsigned char* data,
    int width, int height, 
    int bytesPerPixel,
    int channel,
    bool open
);
