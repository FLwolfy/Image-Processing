#include <image.h>
#include <utils.h>
#include <stdexcept>
#include <cmath>


///////////////////////////////////////////////////////////////
///////////////////////// PUBLIC APIS /////////////////////////
///////////////////////////////////////////////////////////////

////////////// Regular processing functions //////////////

Image Image::ChannelSeparate(const Image& img, int channel)
{
    std::vector<unsigned char> separatedData = SeparateChannel(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel);

    Image separatedImage = Image(img.m_width, img.m_height, 1);
    separatedImage.m_data = separatedData;

    return separatedImage;
}

Image Image::GrayScale(const Image &img)
{
    std::vector<unsigned char> grayData = ToGrayScale(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel);
    
    Image grayImage = Image(img.m_width, img.m_height, 1);
    grayImage.m_data = grayData;

    return grayImage;
}

Image Image::Negative(const Image& img, int channel)
{
    std::vector<unsigned char> negativeData = ToNegative(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel);

    Image negativeImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    negativeImage.m_data = negativeData;

    return negativeImage;
}

Image Image::WaterMark(const Image& img, const Image& watermark, int offsetX, int offsetY, unsigned char threshold, float blendRate)
{
    std::vector<unsigned char> watermarkedData = AddWatermark(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, watermark.m_data.data(), watermark.m_width, watermark.m_height, offsetX, offsetY, threshold, blendRate);

    Image watermarkedImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    watermarkedImage.m_data = watermarkedData;

    return watermarkedImage;
}

////////////// Enhancement functions //////////////

Image Image::LinearScale(const Image& img, int channel, int min, int max)
{
    std::vector<unsigned char> linearScaledData = ToLinearScale(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, min, max);

    Image linearScaledImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    linearScaledImage.m_data = linearScaledData;

    return linearScaledImage;
}

Image Image::HistEqualize(const Image& img, int channel, int binSize)
{
    std::vector<unsigned int> cumulativeHistData = img.GetCumulativeHist();
    std::vector<unsigned char> enhancedData = EqualizeHistogram(img.m_data.data(), cumulativeHistData.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, binSize);

    Image enhancedImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    enhancedImage.m_data = enhancedData;

    return enhancedImage;
}

////////////// Noise Removal functions //////////////

Image Image::MeanDenoise(const Image& img, int channel, int windowSize)
{
    std::vector<unsigned char> denoisedData = MeanFilter(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, windowSize);

    Image denoisedImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    denoisedImage.m_data = denoisedData;

    return denoisedImage;
}

Image Image::MedianDenoise(const Image& img, int channel, int windowSize, bool pseudo)
{
    std::vector<unsigned char> denoisedData = MedianFilter(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, windowSize, pseudo);

    Image denoisedImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    denoisedImage.m_data = denoisedData;

    return denoisedImage;
}

Image Image::GaussianDenoise(const Image& img, int channel, int windowSize, float STD)
{
    std::vector<unsigned char> denoisedData = GaussianFilter(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, windowSize, STD);

    Image denoisedImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    denoisedImage.m_data = denoisedData;

    return denoisedImage;
}

Image Image::BilateralDenoise(const Image& img, int channel, int windowSize, float spaceSTD, float colorSTD)
{
    std::vector<unsigned char> denoisedData = BilateralFilter(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, windowSize, spaceSTD, colorSTD);

    Image denoisedImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    denoisedImage.m_data = denoisedData;

    return denoisedImage;
}

////////////// Edge Detection functions //////////////

Image Image::SobelEdge(const Image& img, int channel, int windowSize, const std::string& suppressedMethod, const std::string& thresholdMethod, const std::unordered_map<std::string, float>& thresholds)
{
    std::vector<unsigned char> edgeData = ToSobelEdge(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, windowSize, suppressedMethod, thresholdMethod, thresholds);

    Image edgeImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    edgeImage.m_data = edgeData;

    return edgeImage;
}

Image Image::LaplacianEdge(const Image& img, int channel, int windowSize, float noise)
{
    std::vector<unsigned char> edgeData = ToLaplacianEdge(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, windowSize, noise);

    Image edgeImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    edgeImage.m_data = edgeData;

    return edgeImage;
}

std::vector<std::vector<std::pair<int, int>>> Image::GetEdgeContours(int minLength)
{
    return FindContours(m_data.data(), m_width, m_height, m_bytesPerPixel, minLength);
}


////////////// Morphological functions //////////////

Image Image::Shrink(const Image& img, int channel, int iterations)
{
    std::vector<unsigned char> morphData = Morpho(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, 0, iterations);

    Image morphImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    morphImage.m_data = morphData;

    return morphImage;
}

Image Image::Thin(const Image& img, int channel, int iterations)
{
    std::vector<unsigned char> morphData = Morpho(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, 1, iterations);

    Image morphImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    morphImage.m_data = morphData;

    return morphImage;
}

Image Image::Skeletonize(const Image& img, int channel, int iterations)
{
    std::vector<unsigned char> morphData = Morpho(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, 2, iterations);

    Image morphImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    morphImage.m_data = morphData;

    return morphImage;
}

Image Image::Erode(const Image& img, int channel, int iterations)
{
    std::vector<unsigned char> morphData = Morpho(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, 3, iterations);

    Image morphImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    morphImage.m_data = morphData;

    return morphImage;
}

Image Image::Dilate(const Image& img, int channel, int iterations)
{
    std::vector<unsigned char> morphData = Morpho(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, 4, iterations);

    Image morphImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    morphImage.m_data = morphData;

    return morphImage;
}

Image Image::Open(const Image& img, int channel)
{
    std::vector<unsigned char> erodedData = Morpho(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, 3, 1);
    std::vector<unsigned char> openedData = Morpho(erodedData.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, 4, 1);

    Image morphImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    morphImage.m_data = openedData;

    return morphImage;
}

Image Image::Close(const Image& img, int channel)
{
    std::vector<unsigned char> dilatedData = Morpho(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, 4, 1);
    std::vector<unsigned char> closedData = Morpho(dilatedData.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, 3, 1);

    Image morphImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    morphImage.m_data = closedData;

    return morphImage;
}

///////////// Digital Halftoning functions /////////////

Image Image::FixedDither(const Image& img, int channel, unsigned char threshold)
{
    std::vector<unsigned char> ditheredData = FixedDithering(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, threshold);

    Image ditheredImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    ditheredImage.m_data = ditheredData;

    return ditheredImage;
}

Image Image::RandomDither(const Image& img, int channel, bool localHash, unsigned long long seed)
{
    std::vector<unsigned char> ditheredData = RandomDithering(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, localHash, seed);

    Image ditheredImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    ditheredImage.m_data = ditheredData;

    return ditheredImage;
}

Image Image::BayerDither(const Image& img, int channel, int windowSize, int numOfLevels)
{
    std::vector<unsigned char> ditheredData = BayerDithering(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, windowSize, numOfLevels);

    Image ditheredImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    ditheredImage.m_data = ditheredData;

    return ditheredImage;
}

Image Image::ClusterDither(const Image& img, int channel, int clusterSize)
{
    std::vector<unsigned char> ditheredData = ClusterDithering(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, clusterSize);

    Image ditheredImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    ditheredImage.m_data = ditheredData;

    return ditheredImage;
}

Image Image::FSEDDither(const Image& img, int channel, const std::string& method, int param, bool serpentine)
{
    std::vector<unsigned char> errorDiffusedData = FloydSteinbergEDD(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, method, param, serpentine);

    Image errorDiffusedImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    errorDiffusedImage.m_data = errorDiffusedData;

    return errorDiffusedImage;
}

///////////// Geometric Modification functions /////////////

Image Image::Rotate(const Image& img, float angle, const std::string& interpolateMethod)
{
    std::vector<unsigned char> rotatedData;

    if (interpolateMethod == "nearest")
    {
        rotatedData = Rotating(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, angle, 0);
    }
    else if (interpolateMethod == "bilinear")
    {
        rotatedData = Rotating(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, angle, 1);
    }
    else
    {
        throw std::invalid_argument("Invalid interpolation method");
    }

    float radian = angle * 3.14159265358979323846f / 180.0f;
    float cosAngle = std::cos(radian);
    float sinAngle = std::sin(radian);

    int newWidth = static_cast<int>(std::abs(img.m_width * cosAngle) + std::abs(img.m_height * sinAngle));
    int newHeight = static_cast<int>(std::abs(img.m_width * sinAngle) + std::abs(img.m_height * cosAngle));
    
    Image rotatedImage = Image(newWidth, newHeight, img.m_bytesPerPixel);
    rotatedImage.m_data = rotatedData;

    return rotatedImage;
}

Image Image::Scale(const Image& img, float scaleX, float scaleY, const std::string& interpolateMethod)
{
    std::vector<unsigned char> scaledData;

    if (interpolateMethod == "nearest")
    {
        scaledData = Scaling(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, scaleX, scaleY, 0);
    }
    else if (interpolateMethod == "bilinear")
    {
        scaledData = Scaling(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, scaleX, scaleY, 1);
    }
    else
    {
        throw std::invalid_argument("Invalid interpolation method");
    }

    int newWidth = static_cast<int>(img.m_width * scaleX);
    int newHeight = static_cast<int>(img.m_height * scaleY);
    
    Image scaledImage = Image(newWidth, newHeight, img.m_bytesPerPixel);
    scaledImage.m_data = scaledData;

    return scaledImage;
}

Image Image::Translate(const Image& img, float offsetX, float offsetY, const std::string& interpolateMethod)
{
    std::vector<unsigned char> translatedData;

    if (interpolateMethod == "nearest")
    {
        translatedData = Translating(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, offsetX, offsetY, 0);
    }
    else if (interpolateMethod == "bilinear")
    {
        translatedData = Translating(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, offsetX, offsetY, 1);
    }
    else
    {
        throw std::invalid_argument("Invalid interpolation method");
    }

    int newWidth = static_cast<int>(img.m_width + std::abs(offsetX));
    int newHeight = static_cast<int>(img.m_height + std::abs(offsetY));

    Image translatedImage = Image(newWidth, newHeight, img.m_bytesPerPixel);
    translatedImage.m_data = translatedData;

    return translatedImage;
}

Image Image::Shear(const Image& img, float shearX, float shearY, const std::string& interpolateMethod)
{
    std::vector<unsigned char> shearedData;

    if (interpolateMethod == "nearest")
    {
        shearedData = Shearing(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, shearX, shearY, 0);
    }
    else if (interpolateMethod == "bilinear")
    {
        shearedData = Shearing(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, shearX, shearY, 1);
    }
    else
    {
        throw std::invalid_argument("Invalid interpolation method");
    }

    int newWidth = static_cast<int>(img.m_width + std::abs(shearX * img.m_height));
    int newHeight = static_cast<int>(img.m_height + std::abs(shearY * img.m_width));

    Image shearedImage = Image(newWidth, newHeight, img.m_bytesPerPixel);
    shearedImage.m_data = shearedData;

    return shearedImage;
}

Image Image::CircleWarp(const Image& img, bool inverse)
{
    std::vector<unsigned char> warpedData;

    if (inverse) 
    {
        warpedData = CircleToSquareWarp(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel);
    }
    else
    {
        warpedData = SquareToCircleWarp(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel);
    }
    
    Image warpedImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    warpedImage.m_data = warpedData;

    return warpedImage;
}

Image Image::PerspectiveWarp(const Image& img, const std::vector<std::pair<int, int>>& dstPoints, const std::vector<std::pair<int, int>>& srcPoints)
{
    std::vector<unsigned char> warpedData;

    if (srcPoints.empty())
    {
        std::vector<std::pair<int, int>> defaultSrcPoints = { {0, 0}, {img.m_width - 1, 0}, {img.m_width - 1, img.m_height - 1}, {0, img.m_height - 1} };
        warpedData = ToPerspectiveWarp(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, dstPoints, defaultSrcPoints);
    }
    else
    {
        warpedData = ToPerspectiveWarp(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, dstPoints, srcPoints);
    }

    Image warpedImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    warpedImage.m_data = warpedData;

    return warpedImage;
}

///////////// Texture Analysis functions /////////////

std::vector<int> Image::TextureCluster(const std::vector<Image>& imgs, int filterSize, int numOfClusters, int numOfIterations)
{
    std::vector<std::vector<float>> featureMatrix;

    for (const Image& img : imgs)
    {
        std::vector<float> feature = TextureFeatureExtract(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, filterSize);
        featureMatrix.push_back(feature);
    }

    std::vector<int> clusterResult = KMEANSClustering(featureMatrix, numOfClusters, numOfIterations);
    return clusterResult;
}

Image Image::TextureSegment(const Image& img, int channel, int filterSize, int patchSize, int numOfClusters, int numOfIterations)
{
    std::vector<unsigned char> segmentedData = TextureSegmentation(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, channel, filterSize, patchSize, numOfClusters, numOfIterations);

    Image segmentedImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    segmentedImage.m_data = segmentedData;

    return segmentedImage;
}

///////////// Feature Extraction functions /////////////

std::vector<Image> Image::Segment(const Image& img, int minArea)
{
    std::vector<std::tuple<const unsigned char*, int, int>> segments = SegmentImage(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, minArea);

    std::vector<Image> segmentedImages;
    for (const std::tuple<const unsigned char*, int, int>& segment : segments)
    {
        const unsigned char* data;
        int width, height;
        std::tie(data, width, height) = segment;

        Image segmentImage = Image(width, height, img.m_bytesPerPixel);
        segmentImage.m_data.assign(data, data + (width * height * img.m_bytesPerPixel));
        segmentedImages.push_back(segmentImage);
    }

    return segmentedImages;
}

float Image::GetAspectRatio()
{
    return static_cast<float>(m_width) / m_height;
}

float Image::GetAreaRate()
{
    return AreaRate(m_data.data(), m_width, m_height, m_bytesPerPixel);
}

float Image::GetPerimeterRate()
{
    return PerimeterRate(m_data.data(), m_width, m_height, m_bytesPerPixel);
}

float Image::GetEulerNumber(bool connectivity4)
{
    return EulerNumber(m_data.data(), m_width, m_height, m_bytesPerPixel, connectivity4);
}

float Image::GetSpatialMoment(int p, int q)
{
    return SpatialMoment(m_data.data(), m_width, m_height, m_bytesPerPixel, p, q);
}

std::pair<float, float> Image::GetCentroid()
{
    return Centroid(m_data.data(), m_width, m_height, m_bytesPerPixel);
}

float Image::GetSymmetry()
{
    return Symmetry(m_data.data(), m_width, m_height, m_bytesPerPixel);
}

float Image::GetCircularity()
{
    return Circularity(m_data.data(), m_width, m_height, m_bytesPerPixel);
}

//////////////////////////////////////////////////////////////////
///////////////////////// Histogram APIS /////////////////////////
//////////////////////////////////////////////////////////////////

std::vector<unsigned int> Image::GetHist() const
{
    std::vector<unsigned int> histData(256 * m_bytesPerPixel, 0);

    for (int y = 0; y < m_height; y++) 
    {
        for (int x = 0; x < m_width; x++) 
        {
            int index = (y * m_width + x) * m_bytesPerPixel;
            for (int i = 0; i < m_bytesPerPixel; i++) 
            {
                histData[i * 256 + m_data[index + i]]++;
            }
        }
    }

    return histData;
}

std::vector<unsigned int> Image::GetCumulativeHist() const
{
    std::vector<unsigned int> histData = GetHist();
    std::vector<unsigned int> cumulativeHistData(256 * m_bytesPerPixel, 0);

    for (int channel = 0; channel < m_bytesPerPixel; channel++) 
    {
        int offset = channel * 256;
        cumulativeHistData[offset] = histData[offset];

        for (int i = 1; i < 256; i++) 
        {
            cumulativeHistData[offset + i] = cumulativeHistData[offset + i - 1] + histData[offset + i];
        }
    }

    return cumulativeHistData;
}