#include <process.h>

#include <utils.h>

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <iostream>
#include <string>
#include <random>

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
    unsigned char noise,
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
                if (noise > 0) 
                {
                    isWhite = true;
                    for (int i = 0; i < channel; i++) 
                    {
                        if (watermark[watermarkIndex + i] < noise) 
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
    int bytesPerPixel,
    int channel
) {
    std::vector<unsigned char> negativeData(width * height * bytesPerPixel);

    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            for (int i = 0; i < bytesPerPixel; i++)
            {
                if (i == channel)
                {
                    negativeData[index + i] = 255 - data[index + i];
                } 
                else 
                {
                    negativeData[index + i] = data[index + i];
                }
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
    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            if (channelMax != channelMin) 
            {
                linearScaledData[index + channel] = static_cast<unsigned char>(
                    min + (data[index + channel] - channelMin) * (max - min) / (channelMax - channelMin)
                );
            } else 
            {
                linearScaledData[index + channel] = static_cast<unsigned char>(min);
            }

            for (int i = 0; i < bytesPerPixel; i++) 
            {
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
            if (cumulativeHist[channel * 256 + k] >= (unsigned int)(binNum * j)) {
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
    channel = channel < 0 ? 0 : (channel >= bytesPerPixel ? bytesPerPixel - 1 : channel);

    // Make sure the window size is odd and not too big
    windowSize = windowSize < 1 ? 1 : (windowSize > 100 ? 100 : windowSize);
    windowSize = windowSize % 2 == 0 ? windowSize + 1 : windowSize;

    return Convolve(
        data, 
        width, height, 
        bytesPerPixel, 
        channel, 
        Kernel::Mean(windowSize)
    );
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
                    denoisedData[index + channel] = Median(window.data(), windowArea, pseudo);
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
    channel = channel < 0 ? 0 : (channel >= bytesPerPixel ? bytesPerPixel - 1 : channel);

    // Make sure the window size is odd and not too big
    windowSize = windowSize < 1 ? 1 : (windowSize > 100 ? 100 : windowSize);
    windowSize = windowSize % 2 == 0 ? windowSize + 1 : windowSize;

    return Convolve(
        data, 
        width, height, 
        bytesPerPixel, 
        channel, 
        Kernel::Gaussian(STD, windowSize)
    );
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
    float spaceCoeff = -0.5f / (spaceSTD * spaceSTD);
    float colorCoeff = -0.5f / (colorSTD * colorSTD);

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

                    float spaceDist = float(i * i + j * j);
                    float colorDist = float(data[index + channel] - data[curIndex + channel]);

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
) {
    std::vector<unsigned char> edgeData(width * height * bytesPerPixel);

    channel = channel < 0 ? 0 : (channel >= bytesPerPixel ? bytesPerPixel - 1 : channel);

    // Make sure the window size is odd and not too big
    windowSize = windowSize < 1 ? 1 : (windowSize > 100 ? 100 : windowSize);
    windowSize = windowSize % 2 == 0 ? windowSize + 1 : windowSize;

    Kernel kernel = Kernel::Sobel(windowSize, true);
    std::vector<float> gradientX = ConvolvePrecise(data, width, height, bytesPerPixel, channel, kernel);

    kernel = Kernel::Sobel(windowSize, false);
    std::vector<float> gradientY = ConvolvePrecise(data, width, height, bytesPerPixel, channel, kernel);

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            float gx = gradientX[index + channel];
            float gy = gradientY[index + channel];
            float gradient = std::sqrt(gx * gx + gy * gy);
            float angle = abs(std::atan2(gy, gx) * 180 / 3.14159265358979323846f);

            // Apply Suppressing
            if (suppressedMethod == "none") {}
            else if (suppressedMethod == "non-maximum") 
            {
                int angleIndex;
                if (angle <= 22.5 || angle > 157.5) angleIndex = 0;  // 0°
                else if (angle > 22.5 && angle <= 67.5) angleIndex = 1;  // 45°
                else if (angle > 67.5 && angle <= 112.5) angleIndex = 2;  // 90°
                else angleIndex = 3;  // 135°

                int windowHalfSize = windowSize / 2;
                std::vector<float> g(windowHalfSize * 2, 0);

                switch (angleIndex)
                {
                    case 0:  // 0°
                        for (int i = -windowHalfSize; i <= windowHalfSize; i++) 
                        {
                            if (i == 0 || x + i < 0 || x + i >= width) { continue; }
                            int curIndex = (y * width + x + i) * bytesPerPixel;
                            g.push_back(std::abs(gradientX[curIndex + channel]));
                        }
                        break;
                    case 1:  // 45°
                        for (int i = -windowHalfSize; i <= windowHalfSize; i++) 
                        {
                            if (i == 0 || x + i < 0 || x + i >= width || y + i < 0 || y + i >= height) { continue; }
                            int curIndex = ((y + i) * width + x + i) * bytesPerPixel;
                            g.push_back(std::abs(gradientX[curIndex + channel]));
                        }
                        break;
                    case 2:  // 90°
                        for (int i = -windowHalfSize; i <= windowHalfSize; i++) 
                        {
                            if (i == 0 || y + i < 0 || y + i >= height) { continue; }
                            int curIndex = ((y + i) * width + x) * bytesPerPixel;
                            g.push_back(std::abs(gradientY[curIndex + channel]));
                        }
                        break;
                    case 3:  // 135°
                        for (int i = -windowHalfSize; i <= windowHalfSize; i++) 
                        {
                            if (i == 0 || x + i < 0 || x + i >= width || y - i < 0 || y - i >= height) { continue; }
                            int curIndex = ((y - i) * width + x + i) * bytesPerPixel;
                            g.push_back(std::abs(gradientX[curIndex + channel]));
                        }
                        break;
                }

                if (!g.empty() && gradient < *std::max_element(g.begin(), g.end())) 
                {
                    gradient = 0;
                }
            }
            else 
            {
                throw std::invalid_argument("Invalid suppressed method");
            }

            // Apply thresholding
            unsigned char result;
            if (thresholdMethod == "auto") 
            {
                float threshold = 0.15f * Max(data, width * height * bytesPerPixel);
                if (gradient < threshold) 
                {
                    result = 0;
                } 
                else 
                {
                    result = 255;
                }
            }
            else if (thresholdMethod == "manual") 
            {
                float threshold = thresholds.at("threshold");
                if (gradient < threshold) 
                {
                    result = 0;
                } 
                else 
                {
                    result = 255;
                }
            }
            else if (thresholdMethod == "hysteresis")
            {
                if (thresholds.find("low_threshold") == thresholds.end() || thresholds.find("high_threshold") == thresholds.end()) 
                {
                    throw std::invalid_argument("Missing low_threshold or high_threshold");
                }

                float lowThreshold = thresholds.at("low_threshold");
                float highThreshold = thresholds.at("high_threshold");

                if (gradient < lowThreshold) 
                {
                    result = 0;
                } 
                else if (gradient > highThreshold) 
                {
                    result = 255;
                } 
                else 
                {
                    // Check the 8 neighbors
                    bool isEdge = false;
                    for (int i = -1; i <= 1; i++) 
                    {
                        for (int j = -1; j <= 1; j++) 
                        {
                            if (i == 0 && j == 0) { continue; }
                            if (x + j < 0 || x + j >= width || y + i < 0 || y + i >= height) { continue; }

                            int curX = std::min(std::max(x + j, 0), width - 1);
                            int curY = std::min(std::max(y + i, 0), height - 1);
                            int curIndex = (curY * width + curX) * bytesPerPixel;

                            if (gradient > abs(gradientX[curIndex + channel]) && gradient > abs(gradientY[curIndex + channel])) 
                            {
                                isEdge = true;
                                break;
                            }
                        }
                    }

                    // If connected to a strong edge, keep the edge
                    if (isEdge) 
                    {
                        result = 255;
                    }
                    // Else, remove the edge
                    else 
                    {
                        result = 0;
                    }
                }
            }
            else 
            {
                throw std::invalid_argument("Invalid threshold method");
            }

            // Apply the result to the edge data
            for (int i = 0; i < bytesPerPixel; i++) 
            {  
                if (i == channel)
                {
                    edgeData[index + i] = result;
                }
                else
                {
                    // Keep the other channels
                    edgeData[index + i] = data[index + i];
                }
            }
        }
    }

    return edgeData;
}

std::vector<unsigned char> ToLaplacianEdge(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int windowSize,
    float threshold
) {
    std::vector<unsigned char> edgeData(width * height * bytesPerPixel);

    channel = channel < 0 ? 0 : (channel >= bytesPerPixel ? bytesPerPixel - 1 : channel);

    // Make sure the window size is odd and not too big
    windowSize = windowSize < 1 ? 1 : (windowSize > 100 ? 100 : windowSize);
    windowSize = windowSize % 2 == 0 ? windowSize + 1 : windowSize;

    Kernel kernel = Kernel::Laplacian(windowSize);
    std::vector<float> laplacian = ConvolvePrecise(data, width, height, bytesPerPixel, channel, kernel);

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            float gradient = std::abs(laplacian[index + channel]);

            // Apply thresholding
            unsigned char result;
            if (gradient <= threshold) 
            {
                result = 0;
            } 
            else 
            {
                result = 255;
            } 

            // Apply the result to the edge data
            for (int i = 0; i < bytesPerPixel; i++) 
            {  
                if (i == channel)
                {
                    edgeData[index + i] = result;
                }
                else
                {
                    // Keep the other channels
                    edgeData[index + i] = data[index + i];
                }
            }
        }
    }

    return edgeData;
}

///////////////////////// Morphological functions /////////////////////////

std::vector<unsigned char> Morpho(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int type,
    int iterations
) {
    std::vector<unsigned char> morphedData(data, data + width * height * bytesPerPixel);
    Kernel::PatternType patternType;

    channel = channel < 0 ? 0 : (channel >= bytesPerPixel ? bytesPerPixel - 1 : channel);
    if (type < 0 || type > 4) 
    {
        throw std::invalid_argument("Invalid morphological type");
    }
    else if (type == 4)
    {
        // Dilation is the same as erosion with the negative image
        patternType = Kernel::PatternType::EROSION;
        morphedData = ToNegative(data, width, height, bytesPerPixel, channel);
    }
    else
    {
        patternType = (Kernel::PatternType)type;
    }

    // Apply the algorithm
    std::vector<Kernel> conditional = Kernel::Patterns(patternType, true);
    std::vector<Kernel> unconditional = Kernel::Patterns(patternType, false);

    for (int iter = 0; iter < iterations; iter++) 
    {
        std::vector<bool> masks = Mask(morphedData.data(), width, height, bytesPerPixel, channel, conditional);
        std::vector<bool> preserved = MaskBool(masks, width, height, unconditional);

        // Apply the result to the morphed data
        for (int y = 0; y < height; y++) 
        {
            for (int x = 0; x < width; x++) 
            {
                int index = (y * width + x) * bytesPerPixel;
                for (int i = 0; i < bytesPerPixel; i++) 
                {
                    if (i == channel) 
                    {
                        // G = X ∩ [!M ∪ P]
                        if (type >= 3)
                        {
                            morphedData[index + i] *= masks[y * width + x];
                        }
                        else
                        {
                            morphedData[index + i] *= (!masks[y * width + x] || preserved[y * width + x]);
                        }     
                    } 
                }
            }
        }
    }

    if (type == 4) 
    {
        // Dilation is the same as erosion with the negative image
        morphedData = ToNegative(morphedData.data(), width, height, bytesPerPixel, channel);
    }
    
    return morphedData;
}

///////////////////////// Digital Halftoning functions /////////////////////////

std::vector<unsigned char> FixedDithering(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    unsigned char threshold
) {
    std::vector<unsigned char> ditheredData(width * height * bytesPerPixel);

    channel = channel < 0 ? 0 : (channel >= bytesPerPixel ? bytesPerPixel - 1 : channel);

    // Apply the fixed dithering algorithm
    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            for (int i = 0; i < bytesPerPixel; i++) 
            {
                if (i == channel) 
                {
                    ditheredData[index + i] = data[index + i] > threshold ? 255 : 0;
                } 
                else 
                {
                    // Keep the other channels
                    ditheredData[index + i] = data[index + i];
                }
            }
        }
    }

    return ditheredData;
}

std::vector<unsigned char> RandomDithering(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    bool localHash,
    unsigned long long seed
) {
    std::vector<unsigned char> ditheredData(width * height * bytesPerPixel);
    std::vector<unsigned char> hashedData(width * height);

    channel = channel < 0 ? 0 : (channel >= bytesPerPixel ? bytesPerPixel - 1 : channel);

    if (localHash) 
    {
        hashedData = Hash(data, width, height, bytesPerPixel, channel, seed);
    }
    else
    {
        std::mt19937 gen(seed);
        std::uniform_int_distribution<> dis(0, 255);
        for (int i = 0; i < width * height; ++i) 
        {
            hashedData[i] = dis(gen);
        }
    }


    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            for (int i = 0; i < bytesPerPixel; i++) 
            {
                if (i == channel) 
                {
                    ditheredData[index + i] = data[index + i] > hashedData[y * width + x] ? 255 : 0;
                } 
                else 
                {
                    // Keep the other channels
                    ditheredData[index + i] = data[index + i];
                }
            }
        }
    }

    return ditheredData;
}

std::vector<unsigned char> ClusterDithering(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int clusterSize
) {
    std::vector<unsigned char> ditheredData(width * height * bytesPerPixel);

    channel = channel < 0 ? 0 : (channel >= bytesPerPixel ? bytesPerPixel - 1 : channel);

    // Apply the cluster dithering algorithm
    for (int y = 0; y < height; y += clusterSize) 
    {
        for (int x = 0; x < width; x += clusterSize) 
        {
            int clusterSum = 0;
            for (int i = 0; i < clusterSize; i++) 
            {
                for (int j = 0; j < clusterSize; j++) 
                {
                    int curX = std::min(std::max(x + j, 0), width - 1);
                    int curY = std::min(std::max(y + i, 0), height - 1);
                    int curIndex = (curY * width + curX) * bytesPerPixel;

                    clusterSum += data[curIndex + channel];
                }
            }

            int clusterAvg = clusterSum / (clusterSize * clusterSize);
            for (int i = 0; i < clusterSize; i++) 
            {
                for (int j = 0; j < clusterSize; j++) 
                {
                    int curX = std::min(std::max(x + j, 0), width - 1);
                    int curY = std::min(std::max(y + i, 0), height - 1);
                    int curIndex = (curY * width + curX) * bytesPerPixel;

                    for (int k = 0; k < bytesPerPixel; k++) 
                    {
                        if (k == channel) 
                        {
                            ditheredData[curIndex + k] = data[curIndex + k] > clusterAvg ? 255 : 0;
                        } 
                        else 
                        {
                            // Keep the other channels
                            ditheredData[curIndex + k] = data[curIndex + k];
                        }
                    }
                }
            }
        }
    }

    return ditheredData;
}

std::vector<unsigned char> BayerDithering(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int bayerSize,
    int numOfLevels
) {
    if (numOfLevels < 2) 
    {
        throw std::invalid_argument("numOfLevels must be at least 2.");
    }

    std::vector<unsigned char> ditheredData(width * height * bytesPerPixel);
    Kernel bayerMatrix = Kernel::BayerThreshold(bayerSize);

    std::vector<unsigned char> levels(numOfLevels);
    for (int i = 0; i < numOfLevels; ++i) 
    {
        levels[i] = static_cast<unsigned char>(i * 255 / (numOfLevels - 1));
    }

    channel = channel < 0 ? 0 : (channel >= bytesPerPixel ? bytesPerPixel - 1 : channel);

    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;

            for (int i = 0; i < bytesPerPixel; i++) 
            {
                if (i == channel) 
                {
                    // Normalize Bayer matrix threshold
                    float bayerThreshold = bayerMatrix[(y % bayerSize) * bayerSize + (x % bayerSize)] * (255 / (numOfLevels - 1));

                    // Map pixel to the nearest level
                    unsigned char pixelValue = data[index + i];
                    unsigned char newLevel = 0;
                    for (int j = 0; j < numOfLevels - 1; ++j) 
                    {
                        if (pixelValue >= levels[j] && pixelValue <= levels[j + 1]) 
                        {
                            newLevel = (pixelValue - levels[j] > bayerThreshold) ? levels[j + 1] : levels[j];
                            break;
                        }
                    }

                    ditheredData[index + i] = newLevel;
                } 
                else 
                {
                    ditheredData[index + i] = data[index + i];
                }
            }
        }
    }

    return ditheredData;
}

std::vector<unsigned char> FloydSteinbergEDD(
    const unsigned char* data,
    int width, int height,
    int bytesPerPixel,
    int channel,
    const std::string& ditherMethod,
    int param,
    bool serpentine
) {
    std::vector<unsigned char> errorDiffusedData(data, data + width * height * bytesPerPixel);
    Kernel FSMatrix = Kernel::FloydSteinberg();
    Kernel BayerMatrix = Kernel::BayerThreshold(2);
    if (param != 2 && ditherMethod == "bayer") BayerMatrix = Kernel::BayerThreshold(param);

    channel = channel < 0 ? 0 : (channel >= bytesPerPixel ? bytesPerPixel - 1 : channel);

    // Apply the Floyd-Steinberg error diffusion algorithm with serpentine scanning
    for (int y = 0; y < height; y++) {
        int startX, endX, step;
        if (serpentine) {
            // Serpentine scanning
            startX = (y % 2 == 0) ? 0 : width - 1;
            endX = (y % 2 == 0) ? width : -1;
            step = (y % 2 == 0) ? 1 : -1;
        } else {
            // Standard scanning
            startX = 0;
            endX = width;
            step = 1;
        }

        for (int x = startX; x != endX; x += step) {
            int index = (y * width + x) * bytesPerPixel;
            unsigned char oldPixel = errorDiffusedData[index + channel];
            unsigned char newPixel;
            if (ditherMethod == "fixed") newPixel = oldPixel > (unsigned char)param ? 255 : 0;
            else if (ditherMethod == "bayer") newPixel = oldPixel > BayerMatrix[(y % 2) * 2 + (x % 2)] * 255 ? 255 : 0;
            else throw std::invalid_argument("Invalid dither method");
            int error = oldPixel - newPixel;

            // Update the pixel with the new value
            errorDiffusedData[index + channel] = newPixel;

            // Diffuse the error to the neighboring pixels
            if (x + 1 < width) 
            {
                int result = errorDiffusedData[index + bytesPerPixel + channel] + error * FSMatrix[5];
                errorDiffusedData[index + bytesPerPixel + channel] = std::min(255, std::max(0, result));
            }
            if (x - 1 >= 0 && y + 1 < height) 
            {
                int result = errorDiffusedData[(y + 1) * width * bytesPerPixel + (x - 1) * bytesPerPixel + channel] + error * FSMatrix[6];
                errorDiffusedData[(y + 1) * width * bytesPerPixel + (x - 1) * bytesPerPixel + channel] = std::min(255, std::max(0, result));
            }
            if (y + 1 < height) 
            {
                int result = errorDiffusedData[(y + 1) * width * bytesPerPixel + x * bytesPerPixel + channel] + error * FSMatrix[7];
                errorDiffusedData[(y + 1) * width * bytesPerPixel + x * bytesPerPixel + channel] = std::min(255, std::max(0, result));
            }
            if (x + 1 < width && y + 1 < height) 
            {
                int result = errorDiffusedData[(y + 1) * width * bytesPerPixel + (x + 1) * bytesPerPixel + channel] + error * FSMatrix[8];
                errorDiffusedData[(y + 1) * width * bytesPerPixel + (x + 1) * bytesPerPixel + channel] = std::min(255, std::max(0, result));
            }
        }
    }

    return errorDiffusedData;
}