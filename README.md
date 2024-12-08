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
  - `negative(img: Image, channel: int = 0)`: Convert the image to its negative.

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
- **[Edge Detection Methods]**
  - `sobel_edge(img: Image, channel: int, window_size: int, suppressed_method: str = "none", threshold_method: str = "auto", thresholds: dict[str, float] = {})`: Apply Sobel edge detection to a specific channel.
  - `laplacian_edge(img: Image, channel: int, window_size: int, noise: float)`: Apply Laplacian edge detection to a specific channel.

- **[Morphological Operations]**
  - `shrink(img: Image, channel: int, iterations: int = 8)`: Apply shrinking to a specific channel.
  - `thin(img: Image, channel: int, iterations: int = 8)`: Apply thinning to a specific channel.
  - `skeletonize(img: Image, channel: int, iterations: int = 8)`: Apply skeletonization to a specific channel.
  - `erode(img: Image, channel: int, iterations: int = 8)`: Apply erosion to a specific channel.
  - `dilate(img: Image, channel: int, iterations: int = 8)`: Apply dilation to a specific channel.
  - `open(img: Image, channel: int)`: Apply opening to a specific channel.
  - `close(img: Image, channel: int)`: Apply closing to a specific channel.

- **[Digital Halftoning Methods]**
  - `fixed_dither(img: Image, channel: int, threshold: int)`: Apply fixed threshold dithering to a specific channel.
  - `random_dither(img: Image, channel: int, local_hash: bool = False, seed: int = 0)`: Apply random dithering to a specific channel.
  - `bayer_dither(img: Image, channel: int, window_size: int, num_of_levels: int = 2)`: Apply Bayer matrix dithering to a specific channel.
  - `cluster_dither(img: Image, channel: int, cluster_size: int)`: Apply clustered-dot dithering to a specific channel.
  - `fsed_dither(img: Image, channel: int, method: str, param: float, serpentine: bool = False)`: Apply Floyd-Steinberg error diffusion dithering to a specific channel.

- **[Geometric Modification Methods]**
  - `rotate(img: Image, angle: float, interpolate_method: str = "nearest")`: Rotate the image by a specified angle.
  - `scale(img: Image, scale_x: float, scale_y: float, interpolate_method: str = "nearest")`: Scale the image by specified factors along the x and y axes.
  - `translate(img: Image, offset_x: int, offset_y: int, interpolate_method: str = "nearest")`: Translate the image by specified offsets along the x and y axes.
  - `shear(img: Image, shear_x: float, shear_y: float, interpolate_method: str = "nearest")`: Shear the image by specified factors along the x and y axes.
  - `circle_warp(img: Image, inverse: bool = false)`: Apply a circular warp to the image.

- **[Texture Analysis Methods]**
  - `texture_cluster(imgs: list[Image], filter_size: int, num_of_clusters: int, num_of_iterations: int)`: Perform texture clustering on a list of images.
  - `texture_segment(img: Image, channel: int, filter_size: int, patch_size: int, num_of_clusters: int, num_of_iterations: int)`: Perform texture segmentation on an image.

- **[Feature Extraction Methods]**
  - `segment(img: Image, min_area: int = 9)`: Segment the image based on a minimum area threshold.
  - `get_aspect_ratio()`: Get the aspect ratio of the image.
  - `get_area_rate()`: Get the area rate of the image.
  - `get_perimeter_rate()`: Get the perimeter rate of the image.
  - `get_euler_number(connectivity4: bool = false)`: Get the Euler number of the image.
  - `get_spatial_moment(p: int, q: int)`: Get the spatial moment of the image.
  - `get_centroid()`: Get the centroid of the image.
  - `get_symmetry()`: Get the symmetry of the image.
  - `get_circularity()`: Get the circularity of the image.

### Plotting Functions (Python)

The `plot.py` module provides functions for visualizing images and histograms.

#### show_images
```python
show_image(image: Image, title: str)
show_images(images: list[Image], subtitles: list[str], title: str = None)
```
Displays the image(s) in a single figure.

#### plot_histograms
```python
plot_histogram(image: Image, title: str, channel: int = 0, cumulative: bool = False)
plot_histograms(images: list[Image], subtitles: list[str], title: str = None, channels: list[int] = None, cumulative: bool = False)
```
Plots histogram(s) for image(s).

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
