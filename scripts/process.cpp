#include <process.h>

#include <utils.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <iostream>

///////////////////////// Regular processing functions /////////////////////////

std::vector<unsigned char> SeparateChannel(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel, 
    int channel
) {
    std::vector<unsigned char> separatedData(width * height);

    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            separatedData[y * width + x] = data[index + channel];
        }
    }

    return separatedData;
}

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

///////////////////////// Enhancement functions /////////////////////////

std::vector<unsigned char> ToLinearScale(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int min, int max
) {
    std::vector<unsigned char> linearScaledData(width * height * bytesPerPixel);

    channel = channel < 0 ? 0 : (channel >= bytesPerPixel ? bytesPerPixel - 1 : channel);

    // Calculate the min and max of the channel
    unsigned char channelMin = 255;
    unsigned char channelMax = 0;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int index = (y * width + x) * bytesPerPixel;
            channelMin = std::min(channelMin, data[index + channel]);
            channelMax = std::max(channelMax, data[index + channel]);
        }
    }

    // Linear scale the channel
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int index = (y * width + x) * bytesPerPixel;
            if (channelMax != channelMin) {
                linearScaledData[index + channel] = static_cast<unsigned char>(
                    min + (data[index + channel] - channelMin) * (max - min) / (channelMax - channelMin)
                );
            } else {
                linearScaledData[index + channel] = static_cast<unsigned char>(min);
            }

            for (int i = 0; i < bytesPerPixel; i++) {
                if (i != channel) {
                    linearScaledData[index + i] = data[index + i];
                }
            }
        }
    }

    return linearScaledData;
}

std::vector<unsigned char> EqualizeHistogram(
    const unsigned char* data, 
    const unsigned int* cumulativeHist,
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int binSize
) {
    std::vector<unsigned char> enhancedData(width * height * bytesPerPixel);

    channel = channel < 0 ? 0 : (channel >= bytesPerPixel ? bytesPerPixel - 1 : channel);

    binSize = binSize < 1 ? 1 : (binSize > 256 ? 256 : binSize);
    int binNum = width * height / binSize;
    binNum = binNum < 1 ? 1 : binNum;

    // Calculate bin color
    std::vector<int> binColor(binSize + 1, 0);
    int flag = 0;
    for (int j = 0; j < binSize; j++) {
        for (int k = flag; k < 256; k++) {
            if (cumulativeHist[channel * 256 + k] >= binNum * j) {
                binColor[j] = k;
                flag = k;
                break;
            }
        }
    }


    // Calculate the enhanced data
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int index = (y * width + x) * bytesPerPixel;

            for (int i = 0; i < bytesPerPixel; i++) {
                if (i == channel) {  
                    int binIndex = cumulativeHist[data[index + i]] / binNum;
                    binIndex = binIndex < binSize ? binIndex : binSize - 1;

                    int toMin = 255 * binIndex / binSize;
                    int toMax = 255 * (binIndex + 1) / binSize;
                    
                    int min = binColor[binIndex];
                    int max = binColor[binIndex + 1];

                    if (min == max) {
                        enhancedData[index + i] = (toMin + toMax) / 2;
                    } else {
                        enhancedData[index + i] = toMin + (data[index + i] - min) * (toMax - toMin) / (max - min);
                    }
                } else {
                    // Keep the other channels
                    enhancedData[index + i] = data[index + i];
                }
            }
        }
    }

    return enhancedData;
}


///////////////////////// Noise Removal functions /////////////////////////

std::vector<unsigned char> MeanFilter(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int windowSize
) {
    std::vector<unsigned char> denoisedData(width * height * bytesPerPixel);

    channel = channel < 0 ? 0 : (channel >= bytesPerPixel ? bytesPerPixel - 1 : channel);

    // Make sure the window size is odd and not too big
    windowSize = windowSize < 1 ? 1 : (windowSize > 100 ? 100 : windowSize);
    windowSize = windowSize % 2 == 0 ? windowSize + 1 : windowSize;

    int halfSize = windowSize / 2;
    int windowArea = windowSize * windowSize;

    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            int sum = 0;

            // Calculate the sum
            for (int i = -halfSize; i <= halfSize; i++) 
            {
                for (int j = -halfSize; j <= halfSize; j++) 
                {
                    // If out of bound, use the nearest pixel
                    int curX = std::min(std::max(x + j, 0), width - 1);
                    int curY = std::min(std::max(y + i, 0), height - 1);
                    int curIndex = (curY * width + curX) * bytesPerPixel;

                    sum += data[curIndex + channel];
                }
            }

            // Calculate the mean
            for (int i = 0; i < bytesPerPixel; i++) 
            {
                if (i == channel) 
                {
                    denoisedData[index + channel] = sum / windowArea;
                } 
                else 
                {
                    // Keep the other channels
                    denoisedData[index + i] = data[index + i];
                }       
            }
        }
    }

    return denoisedData;
}

std::vector<unsigned char> MedianFilter(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int windowSize,
    bool pseudo
) {
    std::vector<unsigned char> denoisedData(width * height * bytesPerPixel);

    channel = channel < 0 ? 0 : (channel >= bytesPerPixel ? bytesPerPixel - 1 : channel);

    // Make sure the window size is odd and not too big
    windowSize = windowSize < 1 ? 1 : (windowSize > 100 ? 100 : windowSize);
    windowSize = windowSize % 2 == 0 ? windowSize + 1 : windowSize;

    int halfSize = windowSize / 2;
    int windowArea = windowSize * windowSize;

    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            std::vector<unsigned char> window(windowArea, 0);

            // Calculate the window
            for (int i = -halfSize; i <= halfSize; i++) 
            {
                for (int j = -halfSize; j <= halfSize; j++) 
                {
                    // If out of bound, use the nearest pixel
                    int curX = std::min(std::max(x + j, 0), width - 1);
                    int curY = std::min(std::max(y + i, 0), height - 1);
                    int curIndex = (curY * width + curX) * bytesPerPixel;

                    window[(i + halfSize) * windowSize + (j + halfSize)] = data[curIndex + channel];
                }
            }

            // Calculate the median
            for (int i = 0; i < bytesPerPixel; i++) 
            {
                if (i == channel) 
                {
                    // Use the pseudo median
                    if (pseudo) 
                    {
                        denoisedData[index + channel] = Maximin(window) / 2 + Minimax(window) / 2;
                    }

                    // Use the real median
                    else 
                    {
                        window = Sort(window);
                        denoisedData[index + channel] = window[windowArea / 2];
                    }
                } 
                else 
                {
                    // Keep the other channels
                    denoisedData[index + i] = data[index + i];
                }
            }
        }
    }

    return denoisedData;
}

std::vector<unsigned char> GaussianFilter(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int windowSize,
    float STD
) {
    std::vector<unsigned char> denoisedData(width * height * bytesPerPixel);

    channel = channel < 0 ? 0 : (channel >= bytesPerPixel ? bytesPerPixel - 1 : channel);

    // Make sure the window size is odd and not too big
    windowSize = windowSize < 1 ? 1 : (windowSize > 100 ? 100 : windowSize);
    windowSize = windowSize % 2 == 0 ? windowSize + 1 : windowSize;

    int halfSize = windowSize / 2;
    int windowArea = windowSize * windowSize;

    // Calculate the Gaussian kernel
    std::vector<float> kernel(windowArea, 0);
    float sum = 0;
    for (int i = -halfSize; i <= halfSize; i++) 
    {
        for (int j = -halfSize; j <= halfSize; j++) 
        {
            kernel[(i + halfSize) * windowSize + (j + halfSize)] = exp(-(i * i + j * j) / (2 * STD * STD));
            sum += kernel[(i + halfSize) * windowSize + (j + halfSize)];
        }
    }

    // Normalize the kernel
    for (int i = 0; i < windowArea; i++) 
    {
        kernel[i] /= sum;
    }

    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            float sum = 0;

            // Calculate the sum
            for (int i = -halfSize; i <= halfSize; i++) 
            {
                for (int j = -halfSize; j <= halfSize; j++) 
                {
                    // If out of bound, use the nearest pixel
                    int curX = std::min(std::max(x + j, 0), width - 1);
                    int curY = std::min(std::max(y + i, 0), height - 1);
                    int curIndex = (curY * width + curX) * bytesPerPixel;

                    sum += data[curIndex + channel] * kernel[(i + halfSize) * windowSize + (j + halfSize)];
                }
            }

            // Calculate the mean
            for (int i = 0; i < bytesPerPixel; i++) 
            {
                if (i == channel) 
                {
                    denoisedData[index + channel] = static_cast<unsigned char>(sum);
                } 
                else 
                {
                    // Keep the other channels
                    denoisedData[index + i] = data[index + i];
                }
            }
        }
    }

    return denoisedData;
}

std::vector<unsigned char> BilateralFilter(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int windowSize,
    float spaceSTD,
    float colorSTD
) {
    std::vector<unsigned char> filteredData(width * height * bytesPerPixel);

    channel = channel < 0 ? 0 : (channel >= bytesPerPixel ? bytesPerPixel - 1 : channel);

    // Make sure the window size is odd and not too big
    windowSize = windowSize < 1 ? 1 : (windowSize > 100 ? 100 : windowSize);
    windowSize = windowSize % 2 == 0 ? windowSize + 1 : windowSize;

    int halfSize = windowSize / 2;
    float spaceCoeff = -0.5 / (spaceSTD * spaceSTD);
    float colorCoeff = -0.5 / (colorSTD * colorSTD);

    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            float sum = 0.0;
            float weightSum = 0.0;

            for (int i = -halfSize; i <= halfSize; i++) 
            {
                for (int j = -halfSize; j <= halfSize; j++) 
                {
                    // If out of bound, use the nearest pixel
                    int curX = std::min(std::max(x + j, 0), width - 1);
                    int curY = std::min(std::max(y + i, 0), height - 1);
                    int curIndex = (curY * width + curX) * bytesPerPixel;

                    float spaceDist = i * i + j * j;
                    float colorDist = data[index + channel] - data[curIndex + channel];

                    float weight = std::exp(spaceDist * spaceCoeff + colorDist * colorCoeff);
                    sum += data[curIndex + channel] * weight;
                    weightSum += weight;
                }
            }

            for (int i = 0; i < bytesPerPixel; i++) 
            {
                if (i == channel) 
                {
                    filteredData[index + channel] = static_cast<unsigned char>(sum / weightSum);
                } 
                else 
                {
                    // Keep the other channels
                    filteredData[index + i] = data[index + i];
                }
            }     
        }
    }

    return filteredData;
}

