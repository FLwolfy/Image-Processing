import matplotlib.pyplot as plt
import textwrap
import numpy as np
from modules.image_processing import Image

def show_image(image: Image, title: str):
    if image.bytes_per_pixel == 1:
        plt.imshow(image.raw_data, cmap="gray", vmin=0, vmax=255)
    else:
        plt.imshow(image.raw_data)
    wrapped_title = "\n".join(textwrap.wrap(title, width=20))
    plt.title(wrapped_title)
    plt.axis('off')
    plt.show()
        
def show_images(images: list[Image], titles: list[str]):
    fig, axes = plt.subplots(1, len(images), figsize=(10, 5)) 
    for i, (image, title) in enumerate(zip(images, titles)):
        if image.bytes_per_pixel == 1:
            axes[i].imshow(image.raw_data, cmap="gray", vmin=0, vmax=255)
        else:
            axes[i].imshow(image.raw_data)
        wrapped_title = "\n".join(textwrap.wrap(title, width=20))
        axes[i].set_title(wrapped_title)
        axes[i].axis('off')
    plt.subplots_adjust(wspace=0.4)        
    plt.show()
    
def plot_histogram(image: Image, title: str, channel: int=0, cumulative: bool=False):
    histogram = image.get_cumulative_hist()[channel] if cumulative else image.get_hist()[channel]
    plt.bar(range(256), histogram, width=1, edgecolor='black')
    wrapped_title = "\n".join(textwrap.wrap(title, width=20))
    plt.title(wrapped_title)
    plt.xlabel('Gray Scale Value')
    plt.ylabel('Number of Pixels')
    plt.show()
    
def plot_histograms(images: list[Image], titles: str, channels: int=None, cumulative: bool=False):
    if channels is None:
        channels = [0] * len(images)
    fig, axes = plt.subplots(1, len(images), figsize=(10, 5))
    for i, (image, title, channel) in enumerate(zip(images, titles, channels)):
        histogram = image.get_cumulative_hist()[channel] if cumulative else image.get_hist()[channel]
        axes[i].bar(range(256), histogram, width=1, edgecolor='black')
        wrapped_title = "\n".join(textwrap.wrap(title, width=20))
        axes[i].set_title(wrapped_title)
        axes[i].set_xlabel('Gray Scale Value')
        axes[i].set_ylabel('Number of Pixels')
    plt.subplots_adjust(wspace=0.4)
    plt.show()
    
    
def tune_param(ref_img: 'Image', target_img: 'Image', func_name: str, param_name: str, param_type: type, param_range: np.ndarray, other_param_dict: dict[str, any], channel: int = 0):
    """
    Plot the Mean Squared Error (MSE) values for a given image method and parameter over a range of values.

    Parameters:
        ref_img: The reference image object.
        target_img: The image object to compare against the reference.
        func_name: A string representing the name of the image method (e.g., "gray_scale").
        param_name: A string representing the name of the parameter to vary.
        param_type: The type of the parameter (e.g., int, float).
        param_range: A list or numpy array representing the range of parameter values to iterate over.
        other_param_dict: A dictionary of other parameters to pass to the method.
        channel: An optional integer specifying the image channel to use for comparison (default is 0).
    """
    
    ref_data = np.array(ref_img.raw_data)[:, :, channel]
    mse_values = []

    for param_value in param_range:
        other_param_dict[param_name] = param_type(param_value)
        transformed_img = getattr(Image, func_name)(target_img, **other_param_dict)
        target_data = np.array(transformed_img.raw_data)[:, :, channel]

        mse = np.mean((ref_data - target_data) ** 2)
        mse_values.append(mse)

    # Plot MSE values
    plt.figure(figsize=(10, 6))
    plt.plot(param_range, mse_values, marker='o', label='MSE')

    min_mse = min(mse_values)
    max_mse = max(mse_values)

    min_idx = mse_values.index(min_mse)
    plt.annotate(f'x: {param_range[min_idx]:.2f}', (param_range[min_idx], min_mse), textcoords="offset points", xytext=(0,10), ha='center', color='red')
    max_idx = mse_values.index(max_mse)
    plt.annotate(f'x: {param_range[max_idx]:.2f}', (param_range[max_idx], max_mse), textcoords="offset points", xytext=(0,10), ha='center', color='blue')
    
    plt.xlabel(f'{param_name} ({param_type.__name__})')
    plt.ylabel('Mean Squared Error (MSE)')
    plt.title(f"MSE vs '{param_name}' for '{func_name}' method")
    plt.grid(True)
    plt.legend()
    plt.show()
    
def compare_images(img1: Image, img2: Image) -> float:
    """
    Compares two images and calculates the Mean Squared Error (MSE) between their raw data.
    
    Parameters:
    - img1: The first image to compare.
    - img2: The second image to compare.
    - channel_count: The number of channels to compare (default is 1 for grayscale).
    
    Returns:
    - float: The Mean Squared Error between the two images.
    """
    data1 = img1.raw_data.flatten()
    data2 = img2.raw_data.flatten()

    # Calculate MSE
    mse = np.mean((data1 - data2) ** 2)
    return mse