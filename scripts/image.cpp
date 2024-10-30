#include <image.h>

Image& Image::GrayScale(Image& img)
{
    unsigned char* grayData = ToGrayScale(img.m_data, img.m_width, img.m_height, img.m_bytesPerPixel);
    
    Image* grayImage = new Image(img.m_width, img.m_height, 1);
    grayImage->SetData(grayData);

    return *grayImage;
}

Image& Image::WaterMark(Image& img, Image& watermark, int offsetX, int offsetY, float filterWhiteThreshold, float blendRate)
{
    unsigned char* watermarkedData = AddWatermark(img.m_data, img.m_width, img.m_height, img.m_bytesPerPixel, watermark.m_data, watermark.m_width, watermark.m_height, offsetX, offsetY, filterWhiteThreshold, blendRate);

    Image* watermarkedImage = new Image(img.m_width, img.m_height, img.m_bytesPerPixel);
    watermarkedImage->SetData(watermarkedData);

    return *watermarkedImage;
}