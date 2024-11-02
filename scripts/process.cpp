#include <process.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <iostream>

std::vector<unsigned char> ToGrayScale(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel
) {
    std::vector<unsigned char> grayData(width * height);

    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;

            if (bytesPerPixel < 3) 
            {
                grayData[y * width + x] = data[index];
            } 
            else 
            {
                grayData[y * width + x] = (unsigned char)(0.299 * data[index] + 0.587 * data[index + 1] + 0.114 * data[index + 2]);
            }
        }
    }

    return grayData;
}

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
) {
    std::vector<unsigned char> watermarkedData(width * height * bytesPerPixel);

    blendRate = blendRate > 1 ? 1 : blendRate < 0 ? 0 : blendRate;
    int channel = bytesPerPixel < 3 ? 1 : 3;

    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            int watermarkIndex = ((y - offsetY) * watermarkWidth + (x - offsetX)) * bytesPerPixel;

            if (y >= offsetY && y < offsetY + watermarkHeight && x >= offsetX && x < offsetX + watermarkWidth) 
            {
                // Calculate if the watermark pixel is white (needed to be filtered)
                bool isWhite = false;
                if (filterWhiteThreshold > 0) 
                {
                    isWhite = true;
                    for (int i = 0; i < channel; i++) 
                    {
                        if (watermark[watermarkIndex + i] < (unsigned char)filterWhiteThreshold) 
                        {
                            isWhite = false;
                            break;
                        }
                    }
                }

                // Blend the watermark pixel with the original pixel
                if (!isWhite)
                {
                    for (int i = 0; i < channel; i++)
                    {
                        watermarkedData[index + i] = (unsigned char)(data[index + i] * (1 - blendRate) + watermark[watermarkIndex + i] * blendRate);
                    }
                }

                // If the watermark pixel is white, keep the original pixel
                else 
                {
                    for (int i = 0; i < channel; i++)
                    {
                        watermarkedData[index + i] = data[index + i];
                    }
                }
            }
            else 
            {
                for (int i = 0; i < bytesPerPixel; i++)
                {
                    watermarkedData[index + i] = data[index + i];
                }
            }
        }
    }

    return watermarkedData;
}

std::vector<unsigned char> ToNegative(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel
) {
    std::vector<unsigned char> negativeData(width * height * bytesPerPixel);

    int channel = bytesPerPixel < 3 ? 1 : 3;

    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            for (int i = 0; i < channel; i++)
            {
                negativeData[index + i] = 255 - data[index + i];
            }
        }
    }

    return negativeData;
}

std::vector<unsigned char> ToLinearScale(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int min, int max
) {
    std::vector<unsigned char> linearScaledData(width * height * bytesPerPixel);
    int channel = bytesPerPixel < 3 ? 1 : 3;

    // Store the min and max value of each channel
    std::vector<int> channelMin(channel, 255);
    std::vector<int> channelMax(channel, 0);

    // Calculate the min and max value of each channel
    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            for (int i = 0; i < channel; i++) 
            {
                channelMin[i] = std::min(channelMin[i], static_cast<int>(data[index + i]));
                channelMax[i] = std::max(channelMax[i], static_cast<int>(data[index + i]));
            }
        }
    }

    // Linear scale each channel
    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            for (int i = 0; i < channel; i++) 
            {
                if (channelMax[i] != channelMin[i]) 
                {
                    linearScaledData[index + i] = static_cast<unsigned char>
                    (
                        min + (data[index + i] - channelMin[i]) * (max - min) / (channelMax[i] - channelMin[i])
                    );
                } 
                else 
                {
                    linearScaledData[index + i] = static_cast<unsigned char>(min);
                }
            }
        }
    }

    return linearScaledData;
}

std::vector<unsigned char> HistoryEqualizedEnhance(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int binSize,
    const unsigned int* cumulativeHist
) {
    std::vector<unsigned char> enhancedData(width * height * bytesPerPixel);
    int channel = bytesPerPixel < 3 ? 1 : 3;
    binSize = binSize < 1 ? 1 : (binSize > 256 ? 256 : binSize);
    int binNum = width * height / binSize;
    binNum = binNum < 1 ? 1 : binNum;

    // Calculate the bin color
    std::vector<int> binColor((binSize + 1) * channel, 0);
    for (int i = 0; i < channel; i++) 
    {
        int flag = 0;

        for (int j = 0; j < binSize; j++)
        {
            for (int k = flag; k < 256; k++) 
            {
                if (cumulativeHist[i * 256 + k] >= binNum * j) 
                {
                    binColor[i * (binSize + 1) + j] = k;
                    flag = k;
                    break;
                }
            }
        }
    }

    // Calculate the equalized image
    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            for (int i = 0; i < channel; i++) 
            {
                int binIndex = cumulativeHist[data[index + i]] / binNum;
                binIndex = binIndex < binSize ? binIndex : binSize - 1;

                int toMin = 255 * binIndex / binSize;
                int toMax = 255 * (binIndex + 1) / binSize;
                
                int min = binColor[i * (binSize + 1) + binIndex];
                int max = binColor[i * (binSize + 1) + binIndex + 1];

                if (min == max) 
                {
                    enhancedData[index + i] = data[index + i];
                }
                else 
                {
                    enhancedData[index + i] = toMin + (data[index + i] - min) * (toMax - toMin) / (max - min);
                }
            }
        }
    }

    return enhancedData;
}
