#include <kernel.h>

#include <vector>

///////////////////////// Constructor /////////////////////////

Kernel::Kernel(std::vector<float> initValues) 
    : values(initValues)
{
    size = static_cast<int>(sqrt(values.size()));

    if (size * size != values.size()) 
    {
        throw std::invalid_argument("Kernel size does not match the number of values.");
    }

    if (size < 3 || size % 2 == 0) 
    {
        throw std::invalid_argument("Kernel size must be an odd number greater than or equal to 3.");
    }
}

Kernel::Kernel(std::initializer_list<float> initValues) 
    : values(initValues)
{
    size = static_cast<int>(sqrt(values.size()));

    if (size * size != values.size()) 
    {
        throw std::invalid_argument("Kernel size does not match the number of values.");
    }

    if (size < 3 || size % 2 == 0) 
    {
        throw std::invalid_argument("Kernel size must be an odd number greater than or equal to 3.");
    }
}

Kernel::Kernel(std::vector<float> initValues, int size) 
    : values(initValues), size(size)
{
    if (size * size != values.size()) 
    {
        throw std::invalid_argument("Kernel size does not match the number of values.");
    }

    if (size < 3 || size % 2 == 0) 
    {
        throw std::invalid_argument("Kernel size must be an odd number greater than or equal to 3.");
    }
}

Kernel::Kernel(std::initializer_list<float> initValues, int size) 
    : values(initValues), size(size)
{
    if (size * size != values.size()) 
    {
        throw std::invalid_argument("Kernel size does not match the number of values.");
    }

    if (size < 3 || size % 2 == 0) 
    {
        throw std::invalid_argument("Kernel size must be an odd number greater than or equal to 3.");
    }
}

///////////////////////// Kernel Generating Functions /////////////////////////

Kernel Kernel::Sobel(int size, bool horizontal) 
{
    if (size != 3 && size != 5) 
    {
        throw std::invalid_argument("Kernel size must be 3 or 5.");
    }

    std::vector<float> kernel(size * size, 0);

    if (horizontal) 
    {
        if (size == 3) 
        {
            kernel = {
                1,  0, -1,
                2,  0, -2,
                1,  0, -1
            };
        } 
        else if (size == 5) 
        {
            kernel = {
                1,  1,  0, -1, -1,
                2,  2,  0, -2, -2,
                3,  3,  0, -3, -3,
                2,  2,  0, -2, -2,
                1,  1,  0, -1, -1
            };
        }
    } 
    else 
    {
        if (size == 3) 
        {
            kernel = {
                1,  2,  1,
                0,  0,  0,
                -1, -2, -1
            };
        } 
        else if (size == 5) 
        {
            kernel = {
                1,  2,  3,  2,  1,
                1,  2,  3,  2,  1,
                0,  0,  0,  0,  0,
                -1, -2, -3, -2, -1,
                -1, -2, -3, -2, -1
            };
        }
    }

    return Kernel(kernel, size);
}

Kernel Kernel::Laplacian(int size) 
{
    if (size != 3 && size != 5) 
    {
        throw std::invalid_argument("Kernel size must be 3 or 5.");
    }

    std::vector<float> kernel(size * size, 0);

    if (size == 3) 
    {
        kernel = {
            1,  0, -1,
            2,  0, -2,
            1,  0, -1
        };
    } 
    else if (size == 5) 
    {
        kernel = {
            1,  1,  0, -1, -1,
            2,  2,  0, -2, -2,
            3,  3,  0, -3, -3,
            2,  2,  0, -2, -2,
            1,  1,  0, -1, -1
        };
    }
    
    return Kernel(kernel, size);
}

Kernel Kernel::Mean(int size) 
{
    if (size < 3 || size % 2 == 0) 
    {
        throw std::invalid_argument("Kernel size must be an odd number greater than or equal to 3.");
    }

    std::vector<float> kernel(size * size, 1.0 / (size * size));

    return Kernel(kernel, size);
}

Kernel Kernel::Gaussian(float sigma, int size) 
{
    if (size < 3 || size % 2 == 0) 
    {
        throw std::invalid_argument("Kernel size must be an odd number greater than or equal to 3.");
    }

    int halfSize = size / 2;
    std::vector<float> kernel(size * size, 0);
    float sum = 0;

    for (int i = -halfSize; i <= halfSize; i++) 
    {
        for (int j = -halfSize; j <= halfSize; j++) 
        {
            kernel[(i + halfSize) * size + (j + halfSize)] = exp(-(i * i + j * j) / (2 * sigma * sigma));
            sum += kernel[(i + halfSize) * size + (j + halfSize)];
        }
    }

    // Normalize the kernel
    for (int i = 0; i < size * size; i++) 
    {
        kernel[i] /= sum;
    }

    return Kernel(kernel, size);
}