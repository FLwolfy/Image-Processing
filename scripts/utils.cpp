#include <utils.h>
#include <functional>

unsigned char Max(const std::vector<unsigned char>& array) 
{
    char max = 0;

    for (unsigned char value : array) 
    {
        max = value > max ? value : max;
    }

    return max;
}

unsigned char Min(const std::vector<unsigned char>& array) 
{
    unsigned char min = 255;

    for (unsigned char value : array) 
    {
        min = value < min ? value : min;
    }

    return min;
}

unsigned char Maximin(const std::vector<unsigned char>& array) 
{
    int length = array.size();
    int half = (length + 1) / 2;
    int step = length - half + 1;
    std::vector<unsigned char> minArray(step);
    std::vector<unsigned char> tempBuffer(half);

    for (int i = 0; i < step; ++i) 
    {
        std::copy(array.begin() + i, array.begin() + i + half, tempBuffer.begin());
        minArray[i] = Min(tempBuffer);
    }

    return Max(minArray);
}

unsigned char Minimax(const std::vector<unsigned char>& array) 
{
    int length = array.size();
    int half = (length + 1) / 2;
    int step = length - half + 1;
    std::vector<unsigned char> maxArray(step);
    std::vector<unsigned char> tempBuffer(half);

    for (int i = 0; i < step; ++i) 
    {
        std::copy(array.begin() + i, array.begin() + i + half, tempBuffer.begin());
        maxArray[i] = Max(tempBuffer);
    }

    return Min(maxArray);
}

std::vector<unsigned char> Sort(const std::vector<unsigned char>& array) 
{
    std::vector<unsigned char> copy = array;

    auto partition = [](std::vector<unsigned char>& vec, int low, int high) -> int 
    {
        unsigned char pivot = vec[high];
        int i = low - 1;

        for (int j = low; j < high; ++j) 
        {
            if (vec[j] < pivot) 
            {
                ++i;
                std::swap(vec[i], vec[j]);
            }
        }

        std::swap(vec[i + 1], vec[high]);
        return i + 1;
    };

    std::function<void(int, int)> quickSort = [&](int low, int high) 
    {
        if (low < high) 
        {
            int pi = partition(copy, low, high);
            quickSort(low, pi - 1);
            quickSort(pi + 1, high);
        }
    };

    if (!array.empty())
    {
        quickSort(0, array.size() - 1);
    }

    return copy;
}