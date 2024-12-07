#include <utils.h>
#include <functional>
#include <iostream>
#include <random>

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

std::vector<unsigned char> Hash(
    const unsigned char* data,
    int width, int height,
    int bytesPerPixel,
    int channel,
    unsigned long long seed
) {
    std::vector<unsigned char> hashedData(width * height);

    // A prime number multiplier for mixing
    const unsigned long long prime = 0x100000001B3;

    // Loop through every pixel's channel
    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            int idx = (y * width + x) * bytesPerPixel + channel;
            unsigned char pixelValue = data[idx];

            unsigned long long hash = seed;  // Initialize hash with seed value

            // Mix pixel value with current hash using bitwise operations
            hash ^= pixelValue;  // XOR with the pixel value to mix
            hash *= prime;       // Multiply with a prime to spread out bits
            
            // Add coordinates (x, y) to further disperse the hash
            hash ^= (x + y * width);  // XOR with the coordinate value to include spatial information
            
            // Additional randomizing transformations
            hash += (hash << 21);  // Left shift and add to introduce non-linearity
            hash ^= (hash >> 35);  // Right shift and XOR to break patterns
            hash *= prime;         // Again multiply to disperse bits
            
            // Mask the hash to prevent overflow
            hash ^= (hash >> 33);
            hash *= prime;

            // Store the resulting hash (keep it within byte range)
            hashedData[y * width + x] = static_cast<unsigned char>(hash & 0xFF);
        }
    }

    return hashedData;
}

std::vector<unsigned char> Convolve(
    const unsigned char* data,
    int width, int height,
    int bytesPerPixel,
    int channel,
    const Kernel& kernel
) {
    int kernelSize = kernel.m_size;
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
    int kernelSize = kernel.m_size;
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

std::vector<bool> Mask(
    const unsigned char* data,
    int width, int height,
    int bytesPerPixel,
    int channel,
    const std::vector<Kernel>& kernels
) {
    std::vector<bool> mask(width * height, false);

    for (const Kernel& kernel : kernels) 
    {
        int kernelSize = kernel.m_size;
        int halfSize = kernelSize / 2;

        for (int y = 0; y < height; y++) 
        {
            for (int x = 0; x < width; x++) 
            {
                bool match = true;

                for (int i = -halfSize; i <= halfSize; i++) 
                {
                    for (int j = -halfSize; j <= halfSize; j++) 
                    {
                        // Handle boundary cases
                        int curX = std::min(std::max(x + j, 0), width - 1);
                        int curY = std::min(std::max(y + i, 0), height - 1);
                        int curIndex = (curY * width + curX) * bytesPerPixel;

                        // Convert pixel value to binary
                        int pixelValue = data[curIndex + channel] == 0 ? 0 : 1;

                        if (pixelValue != kernel[(i + halfSize) * kernelSize + (j + halfSize)]) 
                        {
                            match = false;
                            break;
                        }
                    }
                    if (!match) break;
                }

                if (match) 
                {
                    mask[y * width + x] = true;
                }
            }
        }
    }

    return mask;
}

std::vector<bool> MaskBool(
    const std::vector<bool>& data,
    int width, int height,
    const std::vector<Kernel>& kernels
) {
    std::vector<bool> mask(width * height, false);

    for (const Kernel& kernel : kernels) 
    {
        int kernelSize = kernel.m_size;
        int halfSize = kernelSize / 2;

        for (int y = 0; y < height; y++) 
        {
            for (int x = 0; x < width; x++) 
            {
                bool match = true;

                for (int i = -halfSize; i <= halfSize; i++) 
                {
                    for (int j = -halfSize; j <= halfSize; j++) 
                    {
                        int curX = std::min(std::max(x + j, 0), width - 1);
                        int curY = std::min(std::max(y + i, 0), height - 1);
                        int curIndex = curY * width + curX;

                        if (data[curIndex] != kernel[(i + halfSize) * kernelSize + (j + halfSize)]) 
                        {
                            match = false;
                            break;
                        }
                    }
                    if (!match) break;
                }

                if (match) 
                {
                    mask[y * width + x] = true;
                }
            }
        }
    }

    return mask;
}

std::vector<int> kMeansClustering(
    const std::vector<std::vector<float>>& data,
    int numOfClusters
) {
    int numOfPoints = data.size();
    int numOfFeatures = data[0].size();

    // 计算两点之间的欧几里得距离
    auto euclideanDistance = [](const std::vector<float>& a, const std::vector<float>& b) -> float {
        float distance = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            distance += (a[i] - b[i]) * (a[i] - b[i]);
        }
        return std::sqrt(distance);
    };

    // 初始化聚类中心（随机选择数据点作为初始中心）
    std::vector<std::vector<float>> centroids(numOfClusters, std::vector<float>(numOfFeatures, 0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, numOfPoints - 1);
    
    // 随机选择初始中心
    for (int i = 0; i < numOfClusters; ++i) {
        centroids[i] = data[dis(gen)];
    }

    std::vector<int> clusterAssignments(numOfPoints, -1);  // 每个数据点的聚类编号
    std::vector<int> pointsInCluster(numOfClusters, 0);     // 每个聚类中点的数量
    std::vector<std::vector<float>> newCentroids(numOfClusters, std::vector<float>(numOfFeatures, 0));

    bool converged = false;
    while (!converged) {
        converged = true;
        
        // Step 1: 为每个点分配聚类
        std::fill(pointsInCluster.begin(), pointsInCluster.end(), 0);
        for (int i = 0; i < numOfPoints; ++i) {
            float minDistance = std::numeric_limits<float>::max();
            int bestCluster = -1;
            for (int j = 0; j < numOfClusters; ++j) {
                float dist = euclideanDistance(data[i], centroids[j]);
                if (dist < minDistance) {
                    minDistance = dist;
                    bestCluster = j;
                }
            }
            
            // Step 2: 更新聚类分配
            if (clusterAssignments[i] != bestCluster) {
                converged = false;
                clusterAssignments[i] = bestCluster;
            }
            pointsInCluster[bestCluster]++;
        }

        // Step 3: 重新计算每个聚类的中心
        std::fill(newCentroids.begin(), newCentroids.end(), std::vector<float>(numOfFeatures, 0));
        for (int i = 0; i < numOfPoints; ++i) {
            int clusterId = clusterAssignments[i];
            for (int j = 0; j < numOfFeatures; ++j) {
                newCentroids[clusterId][j] += data[i][j];
            }
        }

        // 求每个聚类中心的均值
        for (int i = 0; i < numOfClusters; ++i) {
            if (pointsInCluster[i] > 0) {
                for (int j = 0; j < numOfFeatures; ++j) {
                    newCentroids[i][j] /= pointsInCluster[i];
                }
            }
        }

        // 检查是否收敛
        if (!converged) {
            centroids = newCentroids;
        }
    }

    return clusterAssignments;
}
