import os
import json
import sys
import glob
import pybind11
from setuptools import setup, Extension

pybind11_include = pybind11.get_include()
source_files = glob.glob("scripts/*.cpp") + glob.glob("scripts/*.c") + glob.glob("scripts/*.hpp")

# Define C++ extension module
ext_modules = [
    Extension(
        "image_processing",
        source_files,
        include_dirs=[pybind11_include, 'scripts'],
        language="c++",
    ),
]

# Setup
setup(
    name="image_processing",
    ext_modules=ext_modules,
    zip_safe=False,
)