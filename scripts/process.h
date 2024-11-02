#pragma once

#include <vector>

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
    float filterWhiteThreshold,
    float blendRate
);

std::vector<unsigned char> ToNegative(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel
);

std::vector<unsigned char> ToLinearScale(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int min, int max
);

std::vector<unsigned char> HistoryEqualizedEnhance(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int binSize,
    const unsigned int* cumulativeHist
);
