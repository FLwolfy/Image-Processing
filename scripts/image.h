#pragma once

#include <rawIO.h>
#include <process.h>

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <vector>

class Image
{
public:
    Image(int m_width, int m_height, int m_bytesPerPixel)
        : m_width(m_width), m_height(m_height), m_bytesPerPixel(m_bytesPerPixel) 
        { m_data = std::vector<unsigned char>(m_width * m_height * m_bytesPerPixel); }; 
    ~Image() = default;

    // IO functions
    inline void Load(const char* inputPath) { InputRaw(m_data.data(), inputPath, m_width, m_height, m_bytesPerPixel); }
    inline void Save(const char* outputPath) { OutputRaw(m_data.data(), outputPath, m_width, m_height, m_bytesPerPixel); }

    // Histogram functions
    std::vector<unsigned int> GetHist() const;
    std::vector<unsigned int> GetCumulativeHist() const;

    // Regular processing functions
    static Image ChannelSeparate(const Image& img, int channel);
    static Image GrayScale(const Image& img);
    static Image Negative(const Image& img, int channel = 0);
    static Image WaterMark(const Image& img, const Image& watermark, int offsetX, int offsetY, unsigned char threshold, float blendRate);
    
    // Enhancement functions
    static Image LinearScale(const Image& img, int channel, int min = 0, int max = 255);
    static Image HistEqualize(const Image& img, int channel, int binSize);

    // Noise Removal functions
    static Image MeanDenoise(const Image& img, int channel, int windowSize);
    static Image MedianDenoise(const Image& img, int channel, int windowSize, bool pseudo);
    static Image GaussianDenoise(const Image& img, int channel, int windowSize, float STD);
    static Image BilateralDenoise(const Image& img, int channel, int windowSize, float spaceSTD, float colorSTD);

    // Edge Detection functions
    static Image SobelEdge(const Image& img, int channel, int windowSize, const std::string& suppressedMethod = "none", const std::string& thresholdMethod = "auto", const std::unordered_map<std::string, float>& thresholds = {});
    static Image LaplacianEdge(const Image& img, int channel, int windowSize, float noise);
    
    // Morphological functions
    static Image Shrink(const Image& img, int channel, int iterations = 8);
    static Image Thin(const Image& img, int channel, int iterations = 8);
    static Image Skeletonize(const Image& img, int channel, int iterations = 8);
    static Image Erode(const Image& img, int channel, int iterations = 8);
    static Image Dilate(const Image& img, int channel, int iterations = 8);
    static Image Open(const Image& img, int channel);
    static Image Close(const Image& img, int channel);

    // Digital Halftoning functions
    static Image FixedDither(const Image& img, int channel, unsigned char threshold);
    static Image RandomDither(const Image& img, int channel, bool localHash = false, unsigned long long seed = 0);
    static Image BayerDither(const Image& img, int channel, int windowSize, int numOfLevels = 2);
    static Image ClusterDither(const Image& img, int channel, int clusterSize);
    static Image FSEDDither(const Image& img, int channel, const std::string& method, int param, bool serpentine = false);

public:
    std::vector<unsigned char> m_data;
    int m_width;
    int m_height;
    int m_bytesPerPixel;
};

