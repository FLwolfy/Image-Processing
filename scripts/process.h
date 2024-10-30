#pragma once

#include <vector>

std::vector<unsigned char> ToGrayScale(
    const std::vector<unsigned char>& data, 
    int width, int height, 
    int bytesPerPixel
);

std::vector<unsigned char> AddWatermark(
    const std::vector<unsigned char>& data, 
    int width, int height, 
    int bytesPerPixel, 
    const std::vector<unsigned char>& watermark, 
    int watermarkWidth, 
    int watermarkHeight, 
    int offsetX, 
    int offsetY,
    float filterWhiteThreshold,
    float blendRate
);

std::vector<unsigned char> ToNegative(
    const std::vector<unsigned char>& data, 
    int width, int height, 
    int bytesPerPixel
);

std::vector<unsigned char> ToLinearScale(
    const std::vector<unsigned char>& data, 
    int width, int height, 
    int bytesPerPixel,
    int min, int max
);
