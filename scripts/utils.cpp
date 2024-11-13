#include <utils.h>
#include <functional>

///////////////////////// Non-Destructive Functions /////////////////////////

unsigned char Max(const unsigned char* array, int length) 
{
    unsigned char max = 0;

    for (int i = 0; i < length; ++i) 
    {
        max = array[i] > max ? array[i] : max;
    }

    return max;
}

unsigned char Min(const unsigned char* array, int length) 
{
    unsigned char min = 255;

    for (int i = 0; i < length; ++i) 
    {
        min = array[i] < min ? array[i] : min;
    }

    return min;
}

unsigned char Maximin(const unsigned char* array, int length) 
{
    int half = (length + 1) / 2;
    int step = length - half + 1;
    std::vector<unsigned char> minArray(step);
    std::vector<unsigned char> tempBuffer(half);

    for (int i = 0; i < step; ++i) 
    {
        std::copy(array + i, array + i + half, tempBuffer.begin());
        minArray[i] = Min(tempBuffer.data(), half);
    }

    return Max(minArray.data(), step);
}

unsigned char Minimax(const unsigned char* array, int length) 
{
    int half = (length + 1) / 2;
    int step = length - half + 1;
    std::vector<unsigned char> maxArray(step);
    std::vector<unsigned char> tempBuffer(half);

    for (int i = 0; i < step; ++i) 
    {
        std::copy(array + i, array + i + half, tempBuffer.begin());
        maxArray[i] = Max(tempBuffer.data(), half);
    }

    return Min(maxArray.data(), step);
}

unsigned char Median(const unsigned char* array, int length, bool pseudo) 
{
    if (pseudo) 
    {
        return Maximin(array, length) / 2 + Minimax(array, length) / 2;
    }
    else
    {
        std::vector<unsigned char> sortedArray = Sort(array, length);
        return sortedArray[length / 2];
    }
}

///////////////////////// Array Functions /////////////////////////

std::vector<unsigned char> Sort(
    const unsigned char* array,
    int length
) {
    std::vector<unsigned char> sortedArray(array, array + length);

    auto partition = [](std::vector<unsigned char>& arr, int low, int high) -> int 
    {
        unsigned char pivot = arr[high];
        int i = low - 1;

        for (int j = low; j < high; ++j) 
        {
            if (arr[j] < pivot) 
            {
                ++i;
                std::swap(arr[i], arr[j]);
            }
        }

        std::swap(arr[i + 1], arr[high]);
        return i + 1;
    };

    std::function<void(int, int)> quickSort = [&](int low, int high) 
    {
        if (low < high) 
        {
            int pi = partition(sortedArray, low, high);
            quickSort(low, pi - 1);
            quickSort(pi + 1, high);
        }
    };

    if (length > 0)
    {
        quickSort(0, length - 1);
    }

    return sortedArray;
}

std::vector<unsigned char> Convolve(
    const unsigned char* data,
    int width, int height,
    int bytesPerPixel,
    int channel,
    const Kernel& kernel
) {
    int kernelSize = kernel.size;
    int halfSize = kernelSize / 2;

    std::vector<unsigned char> convolvedData(width * height * bytesPerPixel);

    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            float sum = 0;

            // Calculate the sum
            for (int i = -halfSize; i <= halfSize; ++i) 
            {
                for (int j = -halfSize; j <= halfSize; ++j) 
                {
                    // If out of bound, use the nearest pixel
                    int curX = std::min(std::max(x + j, 0), width - 1);
                    int curY = std::min(std::max(y + i, 0), height - 1);
                    int curIndex = (curY * width + curX) * bytesPerPixel;

                    sum += data[curIndex + channel] * kernel[(i + halfSize) * kernelSize + (j + halfSize)];
                }
            }

            // Assign the convolved pixel
            for (int i = 0; i < bytesPerPixel; ++i) 
            {
                if (channel == i)
                {
                    convolvedData[index + i] = static_cast<unsigned char>(sum);
                }
                else
                {
                    convolvedData[index + i] = data[index + i];
                }
            }
        }
    }

    return convolvedData;
}

std::vector<float> ConvolvePrecise(
    const unsigned char* data,
    int width, int height,
    int bytesPerPixel,
    int channel,
    const Kernel& kernel
) {
    int kernelSize = kernel.size;
    int halfSize = kernelSize / 2;

    std::vector<float> convolvedData(width * height * bytesPerPixel);

    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            float sum = 0;

            // Calculate the sum
            for (int i = -halfSize; i <= halfSize; ++i) 
            {
                for (int j = -halfSize; j <= halfSize; ++j) 
                {
                    // If out of bound, use the nearest pixel
                    int curX = std::min(std::max(x + j, 0), width - 1);
                    int curY = std::min(std::max(y + i, 0), height - 1);
                    int curIndex = (curY * width + curX) * bytesPerPixel;

                    sum += data[curIndex + channel] * kernel[(i + halfSize) * kernelSize + (j + halfSize)];
                }
            }

            // Assign the convolved pixel
            for (int i = 0; i < bytesPerPixel; ++i) 
            {
                if (channel == i)
                {
                    convolvedData[index + i] = sum;
                }
                else
                {
                    convolvedData[index + i] = data[index + i];
                }
            }
        }
    }

    return convolvedData;
}
