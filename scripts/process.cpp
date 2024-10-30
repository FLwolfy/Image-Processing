#include <process.h>

#include <stdio.h>
#include <stdlib.h>

unsigned char* ToGrayScale(
    unsigned char* data, 
    int width, int height, 
    int bytesPerPixel
) {
    unsigned char* grayData = (unsigned char*)malloc(width * height);
    if (!grayData) 
    {
        printf("\n Memory allocation failed\n");
        exit(1);
    }

    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            grayData[y * width + x] = (unsigned char)(0.299 * data[index] + 0.587 * data[index + 1] + 0.114 * data[index + 2]);
        }
    }

    return grayData;
}

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
) {
    unsigned char* watermarkedData = (unsigned char*)malloc(width * height * bytesPerPixel);
    if (!watermarkedData) 
    {
        printf("\n Memory allocation failed\n");
        exit(1);
    }

    blendRate = blendRate > 1 ? 1 : blendRate < 0 ? 0 : blendRate;

    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            int watermarkIndex = ((y - offsetY) * watermarkWidth + (x - offsetX)) * bytesPerPixel;

            if (y >= offsetY && y < offsetY + watermarkHeight && x >= offsetX && x < offsetX + watermarkWidth) 
            {
                // Color Image
                if (bytesPerPixel >= 3) 
                {
                    bool isWhite = filterWhiteThreshold > 0 &&
                        watermark[watermarkIndex] >= (unsigned char)filterWhiteThreshold &&
                        watermark[watermarkIndex + 1] >= (unsigned char)filterWhiteThreshold &&
                        watermark[watermarkIndex + 2] >= (unsigned char)filterWhiteThreshold;

                    if (!isWhite)
                    {
                        watermarkedData[index] = (unsigned char)(data[index] * (1 - blendRate) + watermark[watermarkIndex] * blendRate);
                        watermarkedData[index + 1] = (unsigned char)(data[index + 1] * (1 - blendRate) + watermark[watermarkIndex + 1] * blendRate);
                        watermarkedData[index + 2] = (unsigned char)(data[index + 2] * (1 - blendRate) + watermark[watermarkIndex + 2] * blendRate);
                    } 
                    else 
                    {
                        watermarkedData[index] = data[index];
                        watermarkedData[index + 1] = data[index + 1];
                        watermarkedData[index + 2] = data[index + 2];
                    }
                }

                // Gray Scale Image
                else if (bytesPerPixel == 1) 
                {
                    bool isWhite = filterWhiteThreshold > 0 && watermark[watermarkIndex] >= (unsigned char)filterWhiteThreshold;

                    if (!isWhite)
                    {
                        watermarkedData[index] = (unsigned char)(data[index] * (1 - blendRate) + watermark[watermarkIndex] * blendRate);
                    }
                    else 
                    {
                        watermarkedData[index] = data[index];
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
