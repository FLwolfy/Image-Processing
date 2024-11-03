# ImageProcessing Library

## Introduction

The ImageProcessing library provides a set of tools for image manipulation and analysis. It includes a C++ backend for efficient image processing and a Python interface for ease of use. This library allows users to perform various image processing tasks, such as loading, saving, and manipulating images, as well as plotting and analyzing histograms.

## Dependencies
- `pybind11`
- `numpy`
- `matplotlib`

## Environment
- Python >= 3.10
- C++ build environment

## Quick Start
Before run the python code, you should compile the C++ scripts into python package. To build the C++ scripts, run the following command:
```bash
python Compile.py
```
Then you have the package of `image_processing` under your `modules` directory. To use the module, use:
```bash
from modules import Image
```
If you want to add other customized python package, put them under the `modules` directory, and add the following in the `__init__.py` file:
```bash
from .YOUR_MODULE import YOUR_IMPORTS
```

## API Documentation

### Image Class (Python, C++ Backend)

The `Image` class is exposed to Python using `pybind11`. Below are the available methods and properties:

#### Constructor
```python
Image(width: int, height: int, bytesPerPixel: int)
```
- `width`: Width of the image.
- `height`: Height of the image.
- `bytesPerPixel`: Number of bytes per pixel.

#### Properties
- `raw_data`: Returns the raw image data as a NumPy array.
- `width`: Width of the image.
- `height`: Height of the image.
- `bytes_per_pixel`: Number of bytes per pixel.

#### Methods
- **[I/O Methods]**
  - `load(input_file_path: str)`: Load an image from a file.
  - `save(output_file_path: str)`: Save the image to a file.

- **[Histogram Methods]**
  - `get_hist()`: Get the histogram of the image.
  - `get_cumulative_hist()`: Get the cumulative histogram of the image.


#### Static Methods

- **[Channel Operations]**
  - `channel_separate(img: Image, channel: int)`: Separate a specific channel from the image.
  - `gray_scale(img: Image)`: Convert the image to grayscale. 
  - `negative(img: Image)`: Convert the image to its negative.

- **[Image Enhancement]**
  - `water_mark(img: Image, watermark: Image, offset_x: int, offset_y: int, filter_white_threshold: int, blend_rate: float)`: Apply a watermark to the image.
  - `linear_scale(img: Image, channel: int, min: int, max: int)`: Apply linear scaling to a specific channel. 
  - `hist_equalize(img: Image, channel: int, bin_size: int)`: Equalize the histogram of a specific channel.

- **[Denoising Methods]**
  - `mean_denoise(img: Image, channel: int, window_size: int)`: Apply mean denoising to a specific channel. 
  - `median_denoise(img: Image, channel: int, window_size: int, pseudo: bool = False)`: Apply median denoising to a specific channel.
  - `gaussian_denoise(img: Image, channel: int, window_size: int, STD: float)`: Apply Gaussian denoising to a specific channel. 
  - `bilateral_denoise(img: Image, channel: int, window_size: int, space_STD: float, color_STD: float)`: Apply bilateral denoising to a specific channel.
 
---

### Plotting Functions (Python)

The `plot.py` module provides functions for visualizing images and histograms.

#### show_image
```python
show_image(image: Image, title: str)
```
Displays a single image with a title.

#### show_images
```python
show_images(images: list[Image], titles: list[str])
```
Displays multiple images in a single figure.

#### plot_histogram
```python
plot_histogram(image: Image, title: str, channel: int = 0, cumulative: bool = False)
```
Plots the histogram of an image.

#### plot_histograms
```python
plot_histograms(images: list[Image], titles: list[str], channels: list[int] = None, cumulative: bool = False)
```
Plots histograms for multiple images.

#### tune_param
```python
tune_param(ref_img: Image, target_img: Image, func_name: str, param_name: str, param_type: type, param_range: np.ndarray, other_param_dict: dict[str, any], channel: int = 0)
```
Plots the Mean Squared Error (MSE) values for a given image method and parameter over a range of values.

#### compare_images
```python
compare_images(img1: Image, img2: Image) -> float
```
Compares two images and calculates the Mean Squared Error (MSE) between their raw data.

## Usage Examples

### Loading and Saving Images
```python
from modules import Image

img = Image(800, 600, 3)
img.load("path/to/image.png")
img.save("path/to/output.png")
```

---

### Displaying Images
```python
from modules import show_image

show_image(img, "Sample Image")
```

### Plotting Histograms
```python
from modules import plot_histogram

plot_histogram(img, "Histogram", channel=0)
```

### Applying Image Processing Functions
```python
from modules import Image

gray_img = Image.gray_scale(img)
```

### Tuning Parameters
```python
from modules import tune_param

tune_param(ref_img, target_img, "mean_denoise", "window_size", int, np.arange(1, 10), {"channel": 0})
```

### Comparing Images
```python
from modules import compare_images

mse = compare_images(img1, img2)
print(f"Mean Squared Error: {mse}")
```

## Notes
This README provides a basic overview of the ImageProcessing library and its functionalities. For more detailed information, please refer to the source code and documentation.
