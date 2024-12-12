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

std::vector<std::vector<std::pair<int, int>>> FindContours(
    const unsigned char* edges,
    int width, int height,
    int bytesPerPixel,
    int minLength = 30
);

std::vector<unsigned char> ToContours(
    const unsigned char* edges,
    int width, int height,
    int bytesPerPixel,
    int minLength = 30
);

///////////////////////// Morphological functions /////////////////////////

std::vector<unsigned char> Morpho(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int type, // 0: shrink, 1: thin, 2: skeletonize, 3: erosion, 4: dilation
    int iterations
);

///////////////////////// Digital Halftoning functions /////////////////////////

std::vector<unsigned char> FixedDithering(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    unsigned char threshold
);

std::vector<unsigned char> RandomDithering(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    bool localHash,
    unsigned long long seed
);

std::vector<unsigned char> ClusterDithering(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int clusterSize
);

std::vector<unsigned char> BayerDithering(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int bayerSize,
    int numOfLevels
);

std::vector<unsigned char> FloydSteinbergEDD(
    const unsigned char* data,
    int width, int height,
    int bytesPerPixel,
    int channel,
    const std::string& ditherMethod,
    int param, // fixed, then input "threshold"; bayer, then input "size"
    bool serpentine
);

///////////////////////// Geometric Modification functions /////////////////////////

std::vector<unsigned char> Rotating(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    float angle, // clockwise,
    int interpolation // 0: nearest, 1: bilinear
);

std::vector<unsigned char> Scaling(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    float scaleX, 
    float scaleY,
    int interpolation // 0: nearest, 1: bilinear
);

std::vector<unsigned char> Translating(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    float offsetX, 
    float offsetY,
    int interpolation // 0: nearest, 1: bilinear
);

std::vector<unsigned char> Shearing(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    float shearX,
    float shearY,
    int interpolation // 0: nearest, 1: bilinear
);

std::vector<unsigned char> SquareToCircleWarp(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel
);

std::vector<unsigned char> CircleToSquareWarp(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel
);

std::vector<unsigned char> ToPerspectiveWarp(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    const std::vector<std::pair<int, int>>& dstPoints,
    const std::vector<std::pair<int, int>>& srcPoints  
);

///////////////////////// Texture Analysis functions /////////////////////////////

std::vector<float> TextureFeatureExtract(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int filterSize
);

std::vector<unsigned char> TextureSegmentation(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int channel,
    int filterSize,
    int patchSize,
    int numOfClusters,
    int numOfIterations
);

///////////////////////// Feature Extraction functions /////////////////////////////

std::vector<std::tuple<const unsigned char*, int, int>> SegmentImage(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int minArea = 9
);

float AreaRate(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel
);

float PerimeterRate(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel
);

float EulerNumber(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    bool connectivity4 // 4 or 8
);

float Symmetry(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel
);

float SpatialMoment(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel,
    int p, 
    int q
);

std::pair<float, float> Centroid(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel
);

float Circularity(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel
);

float StatisticalMean(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel
);

float StatisticalVariance(
    const unsigned char* data, 
    int width, int height, 
    int bytesPerPixel
);

