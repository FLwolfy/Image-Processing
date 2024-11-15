#include <kernel.h>

#include <cmath>
#include <vector>
#include <stdexcept>
#include <functional>

///////////////////////// Constructor /////////////////////////

Kernel::Kernel(std::vector<float> initValues) 
    : values(initValues)
{
    size = (int)(sqrt((double)values.size()));

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
    size = (int)(sqrt((double)values.size()));

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

///////////////////////// Regular Kernel /////////////////////////

Kernel Kernel::Mean(int size) 
{
    if (size < 3 || size % 2 == 0) 
    {
        throw std::invalid_argument("Kernel size must be an odd number greater than or equal to 3.");
    }

    std::vector<float> kernel(size * size, 1.0f / (float)(size * size));

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
            kernel[(i + halfSize) * size + (j + halfSize)] = (float)exp(-(i * i + j * j) / (2 * sigma * sigma));
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

///////////////////////// Edge Detection Kernel /////////////////////////

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

///////////////////////// Morphological Kernel /////////////////////////

struct ConditionalPattern 
{
    static const std::vector<std::vector<float>> S1;
    static const std::vector<std::vector<float>> S2;
    static const std::vector<std::vector<float>> S3;
    static const std::vector<std::vector<float>> TK4;
    static const std::vector<std::vector<float>> STK4;
    static const std::vector<std::vector<float>> ST5;
    static const std::vector<std::vector<float>> ST6;
    static const std::vector<std::vector<float>> STK6;
    static const std::vector<std::vector<float>> STK7;
    static const std::vector<std::vector<float>> STK8;
    static const std::vector<std::vector<float>> STK9;
    static const std::vector<std::vector<float>> STK10;
    static const std::vector<std::vector<float>> K11;
};

struct UnconditionalPattern 
{
    static const std::vector<std::vector<float>> SP1;
    static const std::vector<std::vector<float>> SP2;
    static const std::vector<std::vector<float>> S4C1;
    static const std::vector<std::vector<float>> S4C2;
    static const std::vector<std::vector<float>> LC1;
    static const std::vector<std::vector<float>> LC2;
    static const std::vector<std::vector<float>> C4O;
    static const std::vector<std::vector<float>> SCC;
    static const std::vector<std::vector<float>> CC1;
    static const std::vector<std::vector<float>> CC2;
    static const std::vector<std::vector<float>> TB1;
    static const std::vector<std::vector<float>> TB2;
    static const std::vector<std::vector<float>> VB;
    static const std::vector<std::vector<float>> DB;
};

// 0 -> 0; 1 -> 1.
const std::vector<std::vector<float>> ConditionalPattern::S1 = {
        {0, 0, 1, 
         0, 1, 0, 
         0, 0, 0},

        {1, 0, 0, 
         0, 1, 0, 
         0, 0, 0},

        {0, 0, 0, 
         0, 1, 0, 
         1, 0, 0},

        {0, 0, 0, 
         0, 1, 0, 
         0, 0, 1}
    };
const std::vector<std::vector<float>> ConditionalPattern::S2 = {
        {0, 0, 0, 
         0, 1, 1, 
         0, 0, 0},

        {0, 0, 0, 
         1, 1, 0, 
         0, 0, 0},

        {0, 1, 0, 
         0, 1, 0, 
         0, 0, 0},

        {0, 0, 0, 
         0, 1, 0, 
         0, 1, 0}
    };
const std::vector<std::vector<float>> ConditionalPattern::S3 = {
        {1, 0, 0, 
         1, 1, 0, 
         0, 0, 0},

        {0, 0, 0, 
         0, 1, 1, 
         0, 0, 1},

        {0, 1, 0, 
         1, 1, 0, 
         0, 0, 0},

        {0, 0, 0, 
         1, 1, 0, 
         0, 1, 0},

        {0, 0, 0, 
         1, 1, 0, 
         1, 0, 0},

        {0, 0, 1, 
         0, 1, 1, 
         0, 0, 0},

        {0, 1, 0, 
         0, 1, 1, 
         0, 0, 0},

        {0, 0, 0, 
         0, 1, 1, 
         0, 1, 0}
    };
const std::vector<std::vector<float>> ConditionalPattern::TK4 = {
        {0, 1, 0, 
         0, 1, 1, 
         0, 0, 0},

        {0, 1, 0, 
         1, 1, 0, 
         1, 0, 0},

        {0, 0, 0, 
         1, 1, 0, 
         0, 1, 0},

        {0, 0, 1, 
         0, 1, 1, 
         0, 1, 0}
    };
const std::vector<std::vector<float>> ConditionalPattern::STK4 = {
        {0, 0, 1, 
         0, 1, 1, 
         0, 0, 1},

        {1, 1, 1, 
         0, 1, 0, 
         0, 0, 0},

        {1, 0, 0, 
         1, 1, 0, 
         1, 0, 0},

        {0, 0, 0, 
         0, 1, 0, 
         1, 1, 1}
    };
const std::vector<std::vector<float>> ConditionalPattern::ST5 = {
        {1, 1, 0, 
         0, 1, 1, 
         0, 0, 0},

        {0, 1, 0, 
         0, 1, 1, 
         0, 0, 1},

        {0, 1, 1, 
         1, 1, 0, 
         0, 0, 0},

        {0, 0, 1, 
         0, 1, 1, 
         0, 1, 0},

        {0, 1, 1, 
         0, 1, 1, 
         0, 0, 0},

        {1, 1, 0, 
         1, 1, 0, 
         0, 0, 0},

        {0, 0, 0, 
         1, 1, 0, 
         1, 1, 0},

        {0, 0, 0, 
         0, 1, 1, 
         0, 1, 1}
    };
const std::vector<std::vector<float>> ConditionalPattern::ST6 = {
        {1, 1, 0, 
         0, 1, 1, 
         0, 0, 1},

        {0, 1, 1, 
         1, 1, 0, 
         1, 0, 0}
    };
const std::vector<std::vector<float>> ConditionalPattern::STK6 = {
        {1, 1, 1, 
         0, 1, 1, 
         0, 0, 0},

        {1, 0, 0, 
         1, 1, 0, 
         1, 1, 0},
         
        {0, 0, 0, 
         0, 1, 1, 
         1, 1, 1},

        {0, 0, 1, 
         0, 1, 1, 
         0, 1, 1},

        {1, 1, 1, 
         1, 1, 0, 
         0, 0, 0},

        {1, 1, 0, 
         1, 1, 0, 
         1, 0, 0},

        {0, 0, 0, 
         1, 1, 0, 
         1, 1, 1},

        {0, 1, 1, 
         0, 1, 1, 
         0, 0, 1}
    };
const std::vector<std::vector<float>> ConditionalPattern::STK7 = {
        {1, 1, 1, 
         0, 1, 1, 
         0, 0, 1},

        {1, 1, 1, 
         1, 1, 0, 
         1, 0, 0},

        {0, 0, 1, 
         0, 1, 1, 
         1, 1, 1},

        {1, 0, 0, 
         1, 1, 0, 
         1, 1, 1}
    };
const std::vector<std::vector<float>> ConditionalPattern::STK8 = {
        {0, 1, 1, 
         0, 1, 1, 
         0, 1, 1},

        {1, 1, 1, 
         1, 1, 1, 
         0, 0, 0},

        {1, 1, 0, 
         1, 1, 0, 
         1, 1, 0},

        {0, 0, 0, 
         1, 1, 1, 
         1, 1, 1}
    };
const std::vector<std::vector<float>> ConditionalPattern::STK9 = {
        {1, 1, 1, 
         0, 1, 1, 
         0, 1, 1},

        {1, 1, 1, 
         1, 1, 1, 
         1, 0, 0},

        {1, 1, 1, 
         1, 1, 0, 
         1, 1, 0},

        {1, 0, 0, 
         1, 1, 1, 
         1, 1, 1},

        {0, 1, 1, 
         0, 1, 1, 
         1, 1, 1},

        {1, 1, 1, 
         1, 1, 1, 
         0, 0, 1},

        {1, 1, 0, 
         1, 1, 0, 
         1, 1, 1},

        {0, 0, 1, 
         1, 1, 1, 
         1, 1, 1}
    };
const std::vector<std::vector<float>> ConditionalPattern::STK10 = {
        {1, 1, 1, 
         0, 1, 1, 
         1, 1, 1},

        {1, 1, 1, 
         1, 1, 1, 
         1, 0, 1},

        {1, 1, 1, 
         1, 1, 0, 
         1, 1, 1},

        {1, 0, 1, 
         1, 1, 1, 
         1, 1, 1}
    };
const std::vector<std::vector<float>> ConditionalPattern::K11 = {
        {1, 1, 1, 
         1, 1, 1, 
         0, 1, 1},

        {1, 1, 1, 
         1, 1, 1, 
         1, 1, 0},

        {1, 1, 0, 
         1, 1, 1, 
         1, 1, 1},

        {0, 1, 1, 
         1, 1, 1, 
         1, 1, 1}
    };

// 0 -> 0; 1 -> 1; 2 -> 0 or 1; 3 -> At least one of them are 1.
const std::vector<std::vector<float>> UnconditionalPattern::SP1 = {
        {0, 0, 1, 
         0, 1, 0, 
         0, 0, 0},

        {1, 0, 0, 
         0, 1, 0, 
         0, 0, 0}
    };
const std::vector<std::vector<float>> UnconditionalPattern::SP2 = {
        {0, 0, 0, 
         0, 1, 0, 
         1, 0, 0},

        {0, 0, 0, 
         0, 1, 0, 
         0, 0, 1}
    };
const std::vector<std::vector<float>> UnconditionalPattern::S4C1 = {
        {0, 0, 0, 
         0, 1, 0, 
         0, 1, 0},

        {0, 0, 0, 
         0, 1, 1, 
         0, 0, 0}
    };
const std::vector<std::vector<float>> UnconditionalPattern::S4C2 = {
        {0, 0, 0, 
         1, 1, 0, 
         0, 0, 0},

        {0, 1, 0, 
         0, 1, 0, 
         0, 0, 0}
    };
const std::vector<std::vector<float>> UnconditionalPattern::LC1 = {
        {0, 0, 1, 
         0, 1, 1, 
         0, 0, 0},

        {0, 1, 1, 
         0, 1, 0, 
         0, 0, 0},

        {1, 1, 0, 
         0, 1, 0, 
         0, 0, 0},

        {1, 0, 0, 
         1, 1, 0, 
         0, 0, 0},

        {0, 0, 0, 
         1, 1, 0, 
         1, 0, 0},

        {0, 0, 0, 
         0, 1, 0, 
         1, 1, 0},

        {0, 0, 0, 
         0, 1, 0, 
         0, 1, 1},

        {0, 0, 0, 
         0, 1, 1, 
         0, 0, 1}
    };
const std::vector<std::vector<float>> UnconditionalPattern::LC2 = {
        {0, 0, 0, 
         0, 1, 1, 
         0, 1, 0},

        {0, 1, 0, 
         0, 1, 1, 
         0, 0, 0},

        {0, 1, 0, 
         1, 1, 0, 
         0, 0, 0},

        {0, 0, 0, 
         1, 1, 0, 
         0, 1, 0}
    };
const std::vector<std::vector<float>> UnconditionalPattern::C4O = {
        {0, 1, 1, 
         1, 1, 0, 
         0, 0, 0},
         
        {1, 1, 0, 
         0, 1, 1, 
         0, 0, 0},

        {0, 1, 0, 
         0, 1, 1, 
         0, 0, 1},

        {0, 0, 1, 
         0, 1, 1, 
         0, 1, 0}
    };
const std::vector<std::vector<float>> UnconditionalPattern::SCC = { // Multiples
        {0, 3, 1,
         0, 1, 3,
         1, 0, 0},

        {1, 3, 0,
         3, 1, 0,
         0, 0, 1},

        {0, 0, 1,
         3, 1, 0,
         1, 3, 0},

        {1, 0, 0,
         0, 1, 3,
         0, 3, 1}
    };
const std::vector<std::vector<float>> UnconditionalPattern::CC1 = { // Multiples
        {1, 1, 2,
         1, 1, 2,
         2, 2, 2}
    };
const std::vector<std::vector<float>> UnconditionalPattern::CC2 = { // Multiples
        {2, 2, 2,
         2, 1, 1,
         2, 1, 1}
    };
const std::vector<std::vector<float>> UnconditionalPattern::TB1 = { // Multiples
        {2, 1, 0,
         1, 1, 1,
         2, 0, 0},

        {0, 1, 2,
         1, 1, 1,
         0, 0, 2},

        {0, 0, 2,
         1, 1, 1,
         0, 1, 2},

        {2, 0, 0,
         1, 1, 1,
         2, 1, 0},

        {2, 1, 2,
         1, 1, 0,
         0, 1, 0},

        {0, 1, 0,
         1, 1, 0,
         2, 1, 2},

        {0, 1, 0,
         0, 1, 1,
         2, 1, 2},

        {2, 1, 2,
         0, 1, 1,
         0, 1, 0}
    };
const std::vector<std::vector<float>> UnconditionalPattern::TB2 = { // Multiples
        {2, 1, 2,
         1, 1, 1,
         2, 2, 2},

        {2, 2, 2,
         1, 1, 1,
         2, 1, 2},

        {2, 1, 2,
         1, 1, 2,
         2, 1, 2},

        {2, 1, 2,
         2, 1, 1,
         2, 1, 2}
    };
const std::vector<std::vector<float>> UnconditionalPattern::VB = {  // Multiples
        {1, 2, 1,
         2, 1, 2,
         3, 3, 3},

        {1, 2, 3,
         2, 1, 3,
         1, 2, 3},

        {3, 3, 3,
         2, 1, 2,
         1, 2, 1},

        {3, 2, 1,
         3, 1, 2,
         3, 2, 1}
    };
const std::vector<std::vector<float>> UnconditionalPattern::DB = {  // Multiples
        {2, 1, 0,
         0, 1, 1,
         1, 0, 2},

        {0, 1, 2,
         1, 1, 0,
         2, 0, 1},

        {2, 0, 1,
         1, 1, 0,
         0, 1, 2},

        {1, 0, 2,
         0, 1, 1,
         2, 1, 0}
    };


std::vector<Kernel> Kernel::GenerateAllPatterns(const std::vector<std::vector<float>>& patterns, int size)
{
    std::vector<Kernel> result;

    std::function<void(std::vector<float>, int, bool)> generate = [&](std::vector<float> current, int index, bool hasOne) {
        if (index == current.size()) 
        {
            if (hasOne) 
            {
                result.push_back(Kernel(current, size));
            }
            return;
        }

        float value = current[index];
        if (value == 0.0f || value == 1.0f) 
        {
            // If the value is 0 or 1, keep it as is
            generate(current, index + 1, hasOne);
        } 
        else if (value == 2.0f) 
        {
            // If the value is 2, try both 0 and 1
            current[index] = 0.0f;
            generate(current, index + 1, hasOne);
            current[index] = 1.0f;
            generate(current, index + 1, hasOne);
        } 
        else if (value == 3.0f) 
        {
            // If the value is 3, try both 0 and 1, but need to record if at least one of them is 1
            current[index] = 0.0f;
            generate(current, index + 1, hasOne);
            current[index] = 1.0f;
            generate(current, index + 1, true);
        } 
        else 
        {
            throw std::invalid_argument("Invalid pattern value.");
        }
    };

    for (const std::vector<float>& pattern : patterns) 
    {
        bool hasOne = true;
        for (float value : pattern) 
        {
            if (value == 3.0f) 
            {
                hasOne = false;
                break;
            }
        }

        generate(pattern, 0, hasOne);
    }

    return result;
}

std::vector<Kernel> Kernel::Pattern(const std::string& type, bool conditional) 
{
    std::vector<Kernel> kernels;

    if (conditional)
    {
        if (type == "shrink") 
        {
            kernels.insert(kernels.end(), ConditionalPattern::S1.begin(), ConditionalPattern::S1.end());
            kernels.insert(kernels.end(), ConditionalPattern::S2.begin(), ConditionalPattern::S2.end());
            kernels.insert(kernels.end(), ConditionalPattern::S3.begin(), ConditionalPattern::S3.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK4.begin(), ConditionalPattern::STK4.end());
            kernels.insert(kernels.end(), ConditionalPattern::ST5.begin(), ConditionalPattern::ST5.end());
            kernels.insert(kernels.end(), ConditionalPattern::ST6.begin(), ConditionalPattern::ST6.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK6.begin(), ConditionalPattern::STK6.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK7.begin(), ConditionalPattern::STK7.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK8.begin(), ConditionalPattern::STK8.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK9.begin(), ConditionalPattern::STK9.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK10.begin(), ConditionalPattern::STK10.end());
        } 
        else if (type == "thin") 
        {
            kernels.insert(kernels.end(), ConditionalPattern::TK4.begin(), ConditionalPattern::TK4.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK4.begin(), ConditionalPattern::STK4.end());
            kernels.insert(kernels.end(), ConditionalPattern::ST5.begin(), ConditionalPattern::ST5.end());
            kernels.insert(kernels.end(), ConditionalPattern::ST6.begin(), ConditionalPattern::ST6.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK6.begin(), ConditionalPattern::STK6.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK7.begin(), ConditionalPattern::STK7.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK8.begin(), ConditionalPattern::STK8.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK9.begin(), ConditionalPattern::STK9.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK10.begin(), ConditionalPattern::STK10.end());
        }
        else if (type == "skeletonize")
        {
            kernels.insert(kernels.end(), ConditionalPattern::TK4.begin(), ConditionalPattern::TK4.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK4.begin(), ConditionalPattern::STK4.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK6.begin(), ConditionalPattern::STK6.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK7.begin(), ConditionalPattern::STK7.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK8.begin(), ConditionalPattern::STK8.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK9.begin(), ConditionalPattern::STK9.end());
            kernels.insert(kernels.end(), ConditionalPattern::STK10.begin(), ConditionalPattern::STK10.end());
            kernels.insert(kernels.end(), ConditionalPattern::K11.begin(), ConditionalPattern::K11.end());
        } 
        else 
        {
            throw std::invalid_argument("Invalid conditional pattern type.");
        }
    }
    else
    {
        std::vector<Kernel> SCC = GenerateAllPatterns(UnconditionalPattern::SCC, 3);
        std::vector<Kernel> CC1 = GenerateAllPatterns(UnconditionalPattern::CC1, 3);
        std::vector<Kernel> CC2 = GenerateAllPatterns(UnconditionalPattern::CC2, 3);
        std::vector<Kernel> TB1 = GenerateAllPatterns(UnconditionalPattern::TB1, 3);
        std::vector<Kernel> TB2 = GenerateAllPatterns(UnconditionalPattern::TB2, 3);
        std::vector<Kernel> VB = GenerateAllPatterns(UnconditionalPattern::VB, 3);
        std::vector<Kernel> DB = GenerateAllPatterns(UnconditionalPattern::DB, 3);

        if (type == "shrink") 
        {
            kernels.insert(kernels.end(), UnconditionalPattern::SP1.begin(), UnconditionalPattern::SP1.end());
            kernels.insert(kernels.end(), UnconditionalPattern::S4C1.begin(), UnconditionalPattern::S4C1.end());
            kernels.insert(kernels.end(), UnconditionalPattern::C4O.begin(), UnconditionalPattern::C4O.end());
            kernels.insert(kernels.end(), SCC.begin(), SCC.end());
            kernels.insert(kernels.end(), CC1.begin(), CC1.end());
            kernels.insert(kernels.end(), TB1.begin(), TB1.end());
            kernels.insert(kernels.end(), VB.begin(), VB.end());
            kernels.insert(kernels.end(), DB.begin(), DB.end());
        } 
        else if (type == "thin") 
        {
            kernels.insert(kernels.end(), UnconditionalPattern::SP1.begin(), UnconditionalPattern::SP1.end());
            kernels.insert(kernels.end(), UnconditionalPattern::S4C1.begin(), UnconditionalPattern::S4C1.end());
            kernels.insert(kernels.end(), UnconditionalPattern::LC1.begin(), UnconditionalPattern::LC1.end());
            kernels.insert(kernels.end(), UnconditionalPattern::C4O.begin(), UnconditionalPattern::C4O.end());
            kernels.insert(kernels.end(), SCC.begin(), SCC.end());
            kernels.insert(kernels.end(), CC1.begin(), CC1.end());
            kernels.insert(kernels.end(), TB1.begin(), TB1.end());
            kernels.insert(kernels.end(), VB.begin(), VB.end());
            kernels.insert(kernels.end(), DB.begin(), DB.end());
        } 
        else if (type == "skeletonize") 
        {
            kernels.insert(kernels.end(), UnconditionalPattern::SP1.begin(), UnconditionalPattern::SP1.end());
            kernels.insert(kernels.end(), UnconditionalPattern::SP2.begin(), UnconditionalPattern::SP2.end());
            kernels.insert(kernels.end(), UnconditionalPattern::S4C1.begin(), UnconditionalPattern::S4C1.end());
            kernels.insert(kernels.end(), UnconditionalPattern::S4C2.begin(), UnconditionalPattern::S4C2.end());
            kernels.insert(kernels.end(), UnconditionalPattern::LC2.begin(), UnconditionalPattern::LC2.end());
            kernels.insert(kernels.end(), SCC.begin(), SCC.end());
            kernels.insert(kernels.end(), CC1.begin(), CC1.end());
            kernels.insert(kernels.end(), CC2.begin(), CC2.end());
            kernels.insert(kernels.end(), TB2.begin(), TB2.end());
            kernels.insert(kernels.end(), VB.begin(), VB.end());
            kernels.insert(kernels.end(), DB.begin(), DB.end());
        } 
        else 
        {
            throw std::invalid_argument("Invalid unconditional pattern type.");
        }
    }

    return kernels;
}

