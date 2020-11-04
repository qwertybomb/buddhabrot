#ifndef PNG_HH
#define PNG_HH
#include <png.h>

namespace png
{
    void write_image(char const *const filename, std::uint32_t const *const image_data, std::uint32_t image_width, std::uint32_t image_height)
    {
        /* create a zeroed out png_image struct */
        png_image output_png;
        std::memset(&output_png, 0, sizeof(output_png));
        output_png.version = PNG_IMAGE_VERSION;
        output_png.format = PNG_FORMAT_RGBA;
        output_png.width = image_width;
        output_png.height = image_height;

        /* write the png file */
        png_image_write_to_file(&output_png, filename, 0, image_data, image_width * sizeof(*image_data), nullptr);

        /* cleanup */
        png_image_free(&output_png);
    }
}
#endif