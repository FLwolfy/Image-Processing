#include <utils.h>

#include <functional>
#include <iostream>
#include <random>
#include <math.h>
#include <stdlib.h>
#include <algorithm>

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

unsigned char Mean(const unsigned char* array, int length)
{
    int sum = 0;

    for (int i = 0; i < length; ++i)
    {
        sum += array[i];
    }

    return sum / length;
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

float Variance(const unsigned char* array, int length) 
{
    float mean = Mean(array, length);
    float sum = 0;

    for (int i = 0; i < length; i++) 
    {
        float diff = array[i] - mean;
        sum += diff * diff;
    }

    return sum / length;
}

float MSE(const unsigned char* data1, const unsigned char* data2, int length)
{
    float sum = 0;

    for (int i = 0; i < length; ++i)
    {
        float diff = data1[i] - data2[i];
        sum += diff * diff;
    }

    return sum / length;
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

std::vector<unsigned char> Convolve( // Actually, this is a correlation
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

std::vector<float> ConvolvePrecise( // Actually, this is a correlation
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

int MaskCount(
    const unsigned char* data,
    int width, int height,
    int bytesPerPixel,
    int channel,
    const Kernel& kernel
) {
    int count = 0;
    int kernelSize = kernel.m_size;

    for (int y = 0; y <= height - kernelSize; y++) 
    {
        for (int x = 0; x <= width - kernelSize; x++) 
        {
            bool match = true;

            for (int i = 0; i < kernelSize; i++) 
            {
                for (int j = 0; j < kernelSize; j++) 
                {
                    int curX = x + j;
                    int curY = y + i;
                    int curIndex = (curY * width + curX) * bytesPerPixel;

                    int pixelValue = data[curIndex + channel] == 0 ? 0 : 1;

                    if (pixelValue != kernel[i * kernelSize + j]) 
                    {
                        match = false;
                        break;
                    }
                }
                if (!match) break;
            }

            if (match) 
            {
                count++;
            }
        }
    }

    return count;
}

///////////////////////// K-Means Clustering /////////////////////////

class Point
{
private:
	int id_point, id_cluster;
	std::vector<double> values;
	int total_values;

public:
	Point(int id_point, std::vector<double>& values)
	{
		this->id_point = id_point;
		total_values = values.size();

		for(int i = 0; i < total_values; i++)
			this->values.push_back(values[i]);

		id_cluster = -1;
	}

	int getID()
	{
		return id_point;
	}

	void setCluster(int id_cluster)
	{
		this->id_cluster = id_cluster;
	}

	int getCluster()
	{
		return id_cluster;
	}

	double getValue(int index)
	{
		return values[index];
	}

	int getTotalValues()
	{
		return total_values;
	}

	void addValue(double value)
	{
		values.push_back(value);
	}
};

class Cluster
{
private:
	int id_cluster;
	std::vector<double> central_values;
	std::vector<Point> points;

public:
	Cluster(int id_cluster, Point point)
	{
		this->id_cluster = id_cluster;

		int total_values = point.getTotalValues();

		for(int i = 0; i < total_values; i++)
			central_values.push_back(point.getValue(i));

		points.push_back(point);
	}

	void addPoint(Point point)
	{
		points.push_back(point);
	}

	bool removePoint(int id_point)
	{
		int total_points = points.size();

		for(int i = 0; i < total_points; i++)
		{
			if(points[i].getID() == id_point)
			{
				points.erase(points.begin() + i);
				return true;
			}
		}
		return false;
	}

	double getCentralValue(int index)
	{
		return central_values[index];
	}

	void setCentralValue(int index, double value)
	{
		central_values[index] = value;
	}

	Point getPoint(int index)
	{
		return points[index];
	}

	int getTotalPoints()
	{
		return points.size();
	}

	int getID()
	{
		return id_cluster;
	}
};

class KMeans
{
private:
	int K; // number of clusters
	int total_values, total_points, max_iterations;
	std::vector<Cluster> clusters;

	// return ID of nearest center (uses euclidean distance)
	int getIDNearestCenter(Point point)
	{
		double sum = 0.0, min_dist;
		int id_cluster_center = 0;

		for(int i = 0; i < total_values; i++)
		{
			sum += pow(clusters[0].getCentralValue(i) -
					   point.getValue(i), 2.0);
		}

		min_dist = sqrt(sum);

		for(int i = 1; i < K; i++)
		{
			double dist;
			sum = 0.0;

			for(int j = 0; j < total_values; j++)
			{
				sum += pow(clusters[i].getCentralValue(j) -
						   point.getValue(j), 2.0);
			}

			dist = sqrt(sum);

			if(dist < min_dist)
			{
				min_dist = dist;
				id_cluster_center = i;
			}
		}

		return id_cluster_center;
	}

public:
	KMeans(int K, int total_points, int total_values, int max_iterations)
	{
		this->K = K;
		this->total_points = total_points;
		this->total_values = total_values;
		this->max_iterations = max_iterations;
	}

	std::vector<int> run(std::vector<Point> & points)
	{
		if(K > total_points)
			throw std::invalid_argument("K is set to a value larger than the number of points.");

		std::vector<int> prohibited_indexes;

		// choose K distinct values for the centers of the clusters
		for(int i = 0; i < K; i++)
		{
			while(true)
			{
				int index_point = rand() % total_points;

				if(find(prohibited_indexes.begin(), prohibited_indexes.end(),
						index_point) == prohibited_indexes.end())
				{
					prohibited_indexes.push_back(index_point);
					points[index_point].setCluster(i);
					Cluster cluster(i, points[index_point]);
					clusters.push_back(cluster);
					break;
				}
			}
		}

		int iter = 1;

		while(true)
		{
			bool done = true;

			// associates each point to the nearest center
			for(int i = 0; i < total_points; i++)
			{
				int id_old_cluster = points[i].getCluster();
				int id_nearest_center = getIDNearestCenter(points[i]);

				if(id_old_cluster != id_nearest_center)
				{
					if(id_old_cluster != -1)
						clusters[id_old_cluster].removePoint(points[i].getID());

					points[i].setCluster(id_nearest_center);
					clusters[id_nearest_center].addPoint(points[i]);
					done = false;
				}
			}

			// recalculating the center of each cluster
			for(int i = 0; i < K; i++)
			{
				for(int j = 0; j < total_values; j++)
				{
					int total_points_cluster = clusters[i].getTotalPoints();
					double sum = 0.0;

					if(total_points_cluster > 0)
					{
						for(int p = 0; p < total_points_cluster; p++)
							sum += clusters[i].getPoint(p).getValue(j);
						clusters[i].setCentralValue(j, sum / total_points_cluster);
					}
				}
			}

			if(done == true || iter >= max_iterations)
			{
				break;
			}

			iter++;
		}

		// return the result of the K-Means algorithm
        std::vector<int> result(total_points);
        for (int i = 0; i < total_points; i++) {
            result[i] = points[i].getCluster();
        }

        return result;
	}
};

std::vector<int> KMEANSClustering(
    const std::vector<std::vector<float>>& featureMatrix,
    int K,
    int max_iterations
) {
    int total_points = featureMatrix.size();
    int total_values = featureMatrix[0].size();

    std::vector<Point> points;

    for(int i = 0; i < total_points; i++)
    {
        std::vector<double> values;

        for(int j = 0; j < total_values; j++)
        {
            values.push_back(featureMatrix[i][j]);
        }

        Point p(i, values);
        points.push_back(p);
    }

    KMeans kmeans(K, total_points, total_values, max_iterations);
    return kmeans.run(points);
}

std::vector<std::pair<int, int>> FindContours(
    const unsigned char* edge,
    int width, int height,
    int bytesPerPixel
) {
    std::vector<std::pair<int, int>> contours;

    for (int y = 1; y < height - 1; y++) 
    {
        for (int x = 1; x < width - 1; x++) 
        {
            int index = (y * width + x) * bytesPerPixel;
            if (edge[index] == 255) 
            {
                bool isEdge = false;
                for (int i = -1; i <= 1; i++) 
                {
                    for (int j = -1; j <= 1; j++) 
                    {
                        if (i == 0 && j == 0) { continue; }
                        if (x + j < 0 || x + j >= width || y + i < 0 || y + i >= height) { continue; }

                        int curX = std::min(std::max(x + j, 0), width - 1);
                        int curY = std::min(std::max(y + i, 0), height - 1);
                        int curIndex = (curY * width + curX) * bytesPerPixel;

                        if (edge[index] > edge[curIndex]) 
                        {
                            isEdge = true;
                            break;
                        }
                    }
                }

                if (isEdge) 
                {
                    contours.push_back(std::make_pair(x, y));
                }
            }
        }
    }

    return contours;
}
