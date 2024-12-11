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
        .def_static("negative", &Image::Negative, py::arg("img"), py::arg("channel") = 0)

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
        .def_static("laplacian_edge", &Image::LaplacianEdge, py::arg("img"), py::arg("channel"), py::arg("window_size"), py::arg("noise"))
        .def("get_edge_contours", &Image::GetEdgeContours, py::arg("min_length") = 30)

        // Morphological functions
        .def_static("shrink", &Image::Shrink, py::arg("img"), py::arg("channel"), py::arg("iterations") = 8)
        .def_static("thin", &Image::Thin, py::arg("img"), py::arg("channel"), py::arg("iterations") = 8)
        .def_static("skeletonize", &Image::Skeletonize, py::arg("img"), py::arg("channel"), py::arg("iterations") = 8)
        .def_static("erode", &Image::Erode, py::arg("img"), py::arg("channel"), py::arg("iterations") = 8)
        .def_static("dilate", &Image::Dilate, py::arg("img"), py::arg("channel"), py::arg("iterations") = 8)
        .def_static("open", &Image::Open, py::arg("img"), py::arg("channel"))
        .def_static("close", &Image::Close, py::arg("img"), py::arg("channel"))

        // Digital Halftoning functions
        .def_static("fixed_dither", &Image::FixedDither, py::arg("img"), py::arg("channel"), py::arg("threshold"))
        .def_static("random_dither", &Image::RandomDither, py::arg("img"), py::arg("channel"), py::arg("local_hash") = false, py::arg("seed") = 0)
        .def_static("bayer_dither", &Image::BayerDither, py::arg("img"), py::arg("channel"), py::arg("window_size"), py::arg("num_of_levels") = 2)
        .def_static("cluster_dither", &Image::ClusterDither, py::arg("img"), py::arg("channel"), py::arg("cluster_size"))
        .def_static("fsed_dither", &Image::FSEDDither, py::arg("img"), py::arg("channel"), py::arg("method"), py::arg("param"), py::arg("serpentine") = false)

        // Geometric Modification functions
        .def_static("rotate", &Image::Rotate, py::arg("img"), py::arg("angle"), py::arg("interpolate_method") = "nearest")
        .def_static("scale", &Image::Scale, py::arg("img"), py::arg("scale_x"), py::arg("scale_y"), py::arg("interpolate_method") = "nearest")
        .def_static("translate", &Image::Translate, py::arg("img"), py::arg("offset_x"), py::arg("offset_y"), py::arg("interpolate_method") = "nearest")
        .def_static("shear", &Image::Shear, py::arg("img"), py::arg("shear_x"), py::arg("shear_y"), py::arg("interpolate_method") = "nearest")
        .def_static("circle_warp", &Image::CircleWarp, py::arg("img"), py::arg("inverse") = false)
        .def_static("perspective_warp", &Image::PerspectiveWarp, py::arg("img"), py::arg("dst_points"), py::arg("src_points") = std::vector<std::pair<int, int>>())

        // Texture Analysis functions
        .def_static("texture_cluster", &Image::TextureCluster, py::arg("imgs"), py::arg("filter_size"), py::arg("num_of_clusters"), py::arg("num_of_iterations"))
        .def_static("texture_segment", &Image::TextureSegment, py::arg("img"), py::arg("channel"), py::arg("filter_size"), py::arg("patch_size"), py::arg("num_of_clusters"), py::arg("num_of_iterations"))
    
        // Feature Extraction functions
        .def_static("segment", &Image::Segment, py::arg("img"), py::arg("min_area") = 9)
        .def("get_aspect_ratio", &Image::GetAspectRatio)
        .def("get_area_rate", &Image::GetAreaRate)
        .def("get_perimeter_rate", &Image::GetPerimeterRate)
        .def("get_euler_number", &Image::GetEulerNumber, py::arg("connectivity4") = false)
        .def("get_spatial_moment", &Image::GetSpatialMoment, py::arg("p"), py::arg("q"))
        .def("get_centroid", &Image::GetCentroid)
        .def("get_symmetry", &Image::GetSymmetry)
        .def("get_circularity", &Image::GetCircularity)
    ;
}