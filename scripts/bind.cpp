#include <image.h>

// Ignore Warning, the include path is automatically added in setup.py"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h> 

namespace py = pybind11;

PYBIND11_MODULE(image_processing, m)
{
    py::class_<Image>(m, "Image")
        // Constructor
        .def(py::init<int, int, int>(), py::arg("width"), py::arg("height"), py::arg("bytesPerPixel"))

        // Properties
        .def_property_readonly("raw_data", [](const Image& img) {
            return py::array_t<unsigned char>(
                {img.m_width, img.m_height, img.m_bytesPerPixel},
                {img.m_height * img.m_bytesPerPixel, img.m_bytesPerPixel, 1},
                img.m_data.data()
            );
        })
        .def_readonly("width", &Image::m_width)
        .def_readonly("height", &Image::m_height)
        .def_readonly("bytes_per_pixel", &Image::m_bytesPerPixel)

        // IO functions
        .def("load", &Image::Load, py::arg("input_file_path"))
        .def("save", &Image::Save, py::arg("output_file_path"))

        // Image processing functions
        .def_static("gray_scale", &Image::GrayScale, py::arg("img"))
        .def_static("water_mark", &Image::WaterMark, py::arg("img"), py::arg("watermark"), py::arg("offset_x"), py::arg("offset_y"), py::arg("filter_white_threshold"), py::arg("blend_rate"))
        .def_static("negative", &Image::Negative, py::arg("img"))
        .def_static("linear_scale", &Image::LinearScale, py::arg("img"), py::arg("min"), py::arg("max"));
}