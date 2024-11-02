#include <image.h>

///////////////////////////////////////////////////////////////
///////////////////////// PUBLIC APIS /////////////////////////
///////////////////////////////////////////////////////////////

////////////// Regular processing functions //////////////

Image Image::GrayScale(const Image &img)
{
    std::vector<unsigned char> grayData = ToGrayScale(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel);
    
    Image grayImage = Image(img.m_width, img.m_height, 1);
    grayImage.m_data = grayData;

    return grayImage;
}

Image Image::Negative(const Image& img)
{
    std::vector<unsigned char> negativeData = ToNegative(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel);

    Image negativeImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    negativeImage.m_data = negativeData;

    return negativeImage;
}

Image Image::WaterMark(const Image& img, const Image& watermark, int offsetX, int offsetY, float filterWhiteThreshold, float blendRate)
{
    std::vector<unsigned char> watermarkedData = AddWatermark(img.m_data.data(), img.m_width, img.m_height, img.m_bytesPerPixel, watermark.m_data.data(), watermark.m_width, watermark.m_height, offsetX, offsetY, filterWhiteThreshold, blendRate);

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

////////////////////////////////////////////////////////////////
///////////////////////// PRIVATE APIS /////////////////////////
////////////////////////////////////////////////////////////////

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
        cumulativeHistData[offset] = histData[offset]; // 初始化第一个 bin

        for (int i = 1; i < 256; i++) 
        {
            cumulativeHistData[offset + i] = cumulativeHistData[offset + i - 1] + histData[offset + i];
        }
    }

    return cumulativeHistData;
}