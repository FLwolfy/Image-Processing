#include <image.h>

// Ignore Warning, the include path is automatically added in setup.py"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

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

        // Histogram functions
        .def("get_hist", [](const Image& img) {
            int bins = 256;
            return py::array_t<unsigned int>(
                {img.m_bytesPerPixel, bins},
                {bins * sizeof(unsigned int), sizeof(unsigned int)},
                img.GetHist().data()
            );
        })
        .def("get_cumulative_hist", [](const Image& img) {
            int bins = 256;
            return py::array_t<unsigned int>(
                {img.m_bytesPerPixel, bins},
                {bins * sizeof(unsigned int), sizeof(unsigned int)},
                img.GetCumulativeHist().data()
            );
        })

        // IO functions
        .def("load", &Image::Load, py::arg("input_file_path"))
        .def("save", &Image::Save, py::arg("output_file_path"))

        // Regular processing functions
        .def_static("channel_separate", &Image::ChannelSeparate, py::arg("img"), py::arg("channel"))
        .def_static("gray_scale", &Image::GrayScale, py::arg("img"))
        .def_static("water_mark", &Image::WaterMark, py::arg("img"), py::arg("watermark"), py::arg("offset_x"), py::arg("offset_y"), py::arg("filter_white_threshold"), py::arg("blend_rate"))
        .def_static("negative", &Image::Negative, py::arg("img"))

        // Enhancement functions
        .def_static("linear_scale", &Image::LinearScale, py::arg("img"), py::arg("channel"), py::arg("min"), py::arg("max"))
        .def_static("hist_equalize", &Image::HistEqualize, py::arg("img"), py::arg("channel"), py::arg("bin_size"))

        // Noise Removal functions
        .def_static("mean_denoise", &Image::MeanDenoise, py::arg("img"), py::arg("channel"), py::arg("window_size"))
        .def_static("median_denoise", &Image::MedianDenoise, py::arg("img"), py::arg("channel"), py::arg("window_size"), py::arg("pseudo") = false)
        .def_static("gaussian_denoise", &Image::GaussianDenoise, py::arg("img"), py::arg("channel"), py::arg("window_size"), py::arg("STD"))
        .def_static("bilateral_denoise", &Image::BilateralDenoise, py::arg("img"), py::arg("channel"), py::arg("window_size"), py::arg("space_STD"), py::arg("color_STD"))

        // Edge Detection functions
        .def_static("sobel_edge", &Image::SobelEdge, py::arg("img"), py::arg("channel"), py::arg("window_size"), py::arg("suppressed_method") = "none", py::arg("threshold_method") = "auto", py::arg("thresholds") = std::unordered_map<std::string, float>())
        .def_static("laplacian_edge", &Image::LaplacianEdge, py::arg("img"), py::arg("channel"), py::arg("window_size"), py::arg("noise"));
}