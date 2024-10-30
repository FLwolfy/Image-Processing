#include <image.h>

Image Image::GrayScale(const Image& img)
{
    std::vector<unsigned char> grayData = ToGrayScale(img.m_data, img.m_width, img.m_height, img.m_bytesPerPixel);
    
    Image grayImage = Image(img.m_width, img.m_height, 1);
    grayImage.m_data = grayData;

    return grayImage;
}

Image Image::Negative(const Image& img)
{
    std::vector<unsigned char> negativeData = ToNegative(img.m_data, img.m_width, img.m_height, img.m_bytesPerPixel);

    Image negativeImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    negativeImage.m_data = negativeData;

    return negativeImage;
}

Image Image::WaterMark(const Image& img, const Image& watermark, int offsetX, int offsetY, float filterWhiteThreshold, float blendRate)
{
    std::vector<unsigned char> watermarkedData = AddWatermark(img.m_data, img.m_width, img.m_height, img.m_bytesPerPixel, watermark.m_data, watermark.m_width, watermark.m_height, offsetX, offsetY, filterWhiteThreshold, blendRate);

    Image watermarkedImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    watermarkedImage.m_data = watermarkedData;

    return watermarkedImage;
}

Image Image::LinearScale(const Image& img, int min, int max)
{
    std::vector<unsigned char> linearScaledData = ToLinearScale(img.m_data, img.m_width, img.m_height, img.m_bytesPerPixel, min, max);

    Image linearScaledImage = Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    linearScaledImage.m_data = linearScaledData;

    return linearScaledImage;
}