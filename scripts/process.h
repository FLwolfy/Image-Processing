#pragma once

unsigned char* ToGrayScale(
    unsigned char* data, 
    int width, int height, 
    int bytesPerPixel
);

unsigned char* AddWatermark(
    unsigned char* data, 
    int width, int height, 
    int bytesPerPixel, 
    unsigned char* watermark, 
    int watermarkWidth, 
    int watermarkHeight, 
    int offsetX, 
    int offsetY,
    float filterWhiteThreshold,
    float blendRate
);

