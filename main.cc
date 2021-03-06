#include <random>
#include <chrono>
#include <array>
#include <vector>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <omp.h>

/* png headers */
#include "png.hh"

namespace
{
    using v8sf = float __attribute__((vector_size(32)));
    using v8si = std::int32_t __attribute__((vector_size(32)));
    using heatmap_t = std::uint32_t;

    constexpr float scale = 1.0f;
    constexpr float min_r = -2.0f * scale;
    constexpr float max_r = 0.828f * scale;
    constexpr float min_i = -1.414f * scale;
    constexpr float max_i = 1.414f * scale;

    constexpr v8sf min_r_v8sf = { min_r, min_r, min_r, min_r, min_r, min_r, min_r, min_r };
    constexpr v8sf max_r_v8sf = { max_r, max_r, max_r, max_r, max_r, max_r, max_r, max_r };
    constexpr v8sf min_i_v8sf = { min_i, min_i, min_i, min_i, min_i, min_i, min_i, min_i };
    constexpr v8sf max_i_v8sf = { max_i, max_i, max_i, max_i, max_i, max_i, max_i, max_i };

    constexpr std::int32_t width = 640 * 2;
    constexpr std::int32_t height = 480 * 2;
    constexpr std::int32_t red_iterations = 800;
    constexpr std::int32_t green_iterations = 200;
    constexpr std::int32_t blue_iterations = 50;
    constexpr std::int32_t max_iterations = std::max({ red_iterations, blue_iterations, green_iterations });
    constexpr std::uint64_t total_samples = std::uint64_t(10000) *std::uint64_t(10000);
    std::uint64_t individual_samples = total_samples / (omp_get_max_threads() * 8);

    void mandelbrot(v8sf cx, v8sf cy, std::array<v8sf[2], max_iterations> &orbit, v8si &iterations)
    {
        auto good = [](v8sf cx, v8sf cy) {
            // test H1 and H2
            v8sf c2 = cx * cx + cy * cy;
            v8sf d[2] = { cx + 1, cx };
            v8sf h1 = 256 * c2 * c2 - 96 * c2 + 32 * cx - 3;
            v8sf h2 = 16 * (d[0] * d[0] + d[1] * d[1]) - 1;
            return h1 >= 0 && h2 >= 0;
        };

        v8sf zx = {};
        v8sf zy = {};
        v8sf zx_temp = {};
        iterations = good(cx, cy) ? iterations : max_iterations;
        std::int32_t orbit_counter = 0;
        v8si active = {1, 1, 1, 1, 1, 1, 1, 1};
        for (;
            orbit_counter < max_iterations &&
            iterations[0] < max_iterations &&
            iterations[1] < max_iterations &&
            iterations[2] < max_iterations &&
            iterations[3] < max_iterations &&
            iterations[4] < max_iterations &&
            iterations[5] < max_iterations &&
            iterations[6] < max_iterations &&
            iterations[7] < max_iterations

            && (active[0] != 0 || active[1] != 0
                || active[2] != 0 || active[3] != 0
                || active[4] != 0 || active[5] != 0
                || active[5] != 0 || active[6] != 0
                || active[7] != 0); orbit_counter++, iterations += -active)
        {
            zx_temp = active != 0 ? zx_temp * zx_temp - zy * zy + cx : zx_temp;
            zy = active != 0 ? zx * zy * 2 + cy : zy;
            zx = zx_temp;
            orbit[orbit_counter][0] = zx;
            orbit[orbit_counter][1] = zy;
            active = active != 0 ? __builtin_convertvector(zx * zx * zy * zy < 4, v8si) : active;
        }

        /* if point escaped don't plot it */
        iterations = iterations >= max_iterations ? 0 : iterations;
    }

    void generate(std::vector<heatmap_t> &image, heatmap_t &red_max_value, heatmap_t &green_max_value, heatmap_t &blue_max_value)
    {
        auto map = [](auto x, auto in_min, auto in_max, auto out_min, auto out_max)
        {
            return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
        };

        std::random_device rd;
        std::mt19937 engine{ rd() };
        std::uniform_real_distribution<float> dist_r(min_r, max_r);
        std::uniform_real_distribution<float> dist_i(min_i, max_i);

        alignas(32) std::array<v8sf[2], max_iterations> orbit;

        for (std::uint64_t sample_index = 0; sample_index < individual_samples; ++sample_index)
        {
            v8sf cx = { dist_r(engine), dist_r(engine), dist_r(engine), dist_r(engine),
                       dist_r(engine), dist_r(engine), dist_r(engine), dist_r(engine) };
            v8sf cy = { dist_i(engine), dist_i(engine), dist_i(engine), dist_i(engine),
                       dist_i(engine), dist_i(engine), dist_i(engine), dist_i(engine) };
            v8si iterations = {};
            mandelbrot(cx, cy, orbit, iterations);
            std::int32_t max_iteration = std::max({ iterations[0], iterations[1], iterations[2], iterations[3],
                                              iterations[4], iterations[5], iterations[6], iterations[7] });
            /* add points in orbit to heatmap */
            for (std::int32_t i = 0; i < max_iteration; ++i, iterations -= 1)
            {
                /* aspect ratio scaling */
                if constexpr (width > height)
                {
                    orbit[i][1] *= width / float(height);
                }
                else if (width < height)
                {
                    orbit[i][0] *= height / float(width);
                }

                /* map the point to pixel */
                v8sf row = map(orbit[i][0], min_r_v8sf, max_r_v8sf, (v8sf) {},
                    (v8sf)
                {
                    width - 1, width - 1, width - 1, width - 1,
                        width - 1, width - 1, width - 1, width - 1
                });
                v8sf col = map(orbit[i][1], min_i_v8sf, max_i_v8sf, (v8sf) {},
                    (v8sf)
                {
                    height - 1, height - 1, height - 1, height - 1,
                        height - 1, height - 1, height - 1, height - 1
                });
                /* check if the point did escape */
                v8si out_of_bounds = row < 0 || row >= width || col < 0 || col >= height || iterations <= 0;
                heatmap_t pixel1 = !out_of_bounds[0] && i < red_iterations ? ++image[std::uint32_t(col[0]) * width + row[0]] : 0;
                heatmap_t pixel2 = !out_of_bounds[1] && i < red_iterations ? ++image[std::uint32_t(col[1]) * width + row[1]] : 0;
                heatmap_t pixel3 = !out_of_bounds[2] && i < red_iterations ? ++image[std::uint32_t(col[2]) * width + row[2]] : 0;
                heatmap_t pixel4 = !out_of_bounds[3] && i < red_iterations ? ++image[std::uint32_t(col[3]) * width + row[3]] : 0;
                heatmap_t pixel5 = !out_of_bounds[4] && i < red_iterations ? ++image[std::uint32_t(col[4]) * width + row[4]] : 0;
                heatmap_t pixel6 = !out_of_bounds[5] && i < red_iterations ? ++image[std::uint32_t(col[5]) * width + row[5]] : 0;
                heatmap_t pixel7 = !out_of_bounds[6] && i < red_iterations ? ++image[std::uint32_t(col[6]) * width + row[6]] : 0;
                heatmap_t pixel8 = !out_of_bounds[7] && i < red_iterations ? ++image[std::uint32_t(col[7]) * width + row[7]] : 0;
                red_max_value = std::max({ red_max_value, pixel1, pixel2, pixel3, pixel4,
                                          pixel5, pixel6, pixel7, pixel8 });
                pixel1 = !out_of_bounds[0] && i < green_iterations ? ++image[std::uint32_t(col[0]) * width + row[0] + width * height] : 0;
                pixel2 = !out_of_bounds[1] && i < green_iterations ? ++image[std::uint32_t(col[1]) * width + row[1] + width * height] : 0;
                pixel3 = !out_of_bounds[2] && i < green_iterations ? ++image[std::uint32_t(col[2]) * width + row[2] + width * height] : 0;
                pixel4 = !out_of_bounds[3] && i < green_iterations ? ++image[std::uint32_t(col[3]) * width + row[3] + width * height] : 0;
                pixel5 = !out_of_bounds[4] && i < green_iterations ? ++image[std::uint32_t(col[4]) * width + row[4] + width * height] : 0;
                pixel6 = !out_of_bounds[5] && i < green_iterations ? ++image[std::uint32_t(col[5]) * width + row[5] + width * height] : 0;
                pixel7 = !out_of_bounds[6] && i < green_iterations ? ++image[std::uint32_t(col[6]) * width + row[6] + width * height] : 0;
                pixel8 = !out_of_bounds[7] && i < green_iterations ? ++image[std::uint32_t(col[7]) * width + row[7] + width * height] : 0;
                green_max_value = std::max({ green_max_value, pixel1, pixel2, pixel3, pixel4,
                                            pixel5, pixel6, pixel7, pixel8 });
                pixel1 = !out_of_bounds[0] && i < blue_iterations ? ++image[std::uint32_t(col[0]) * width + row[0] + width * height * 2] : 0;
                pixel2 = !out_of_bounds[1] && i < blue_iterations ? ++image[std::uint32_t(col[1]) * width + row[1] + width * height * 2] : 0;
                pixel3 = !out_of_bounds[2] && i < blue_iterations ? ++image[std::uint32_t(col[2]) * width + row[2] + width * height * 2] : 0;
                pixel4 = !out_of_bounds[3] && i < blue_iterations ? ++image[std::uint32_t(col[3]) * width + row[3] + width * height * 2] : 0;
                pixel5 = !out_of_bounds[4] && i < blue_iterations ? ++image[std::uint32_t(col[4]) * width + row[4] + width * height * 2] : 0;
                pixel6 = !out_of_bounds[5] && i < blue_iterations ? ++image[std::uint32_t(col[5]) * width + row[5] + width * height * 2] : 0;
                pixel7 = !out_of_bounds[6] && i < blue_iterations ? ++image[std::uint32_t(col[6]) * width + row[6] + width * height * 2] : 0;
                pixel8 = !out_of_bounds[7] && i < blue_iterations ? ++image[std::uint32_t(col[7]) * width + row[7] + width * height * 2] : 0;
                blue_max_value = std::max({ blue_max_value, pixel1, pixel2, pixel3, pixel4,
                                           pixel5, pixel6, pixel7, pixel8 });
            }
        }

    }

    std::uint32_t get_color(heatmap_t const *value, heatmap_t red_max, heatmap_t blue_max, heatmap_t green_max)
    {
        std::uint32_t red = std::min(*value / double(red_max) * 2 * 255.0 + 0.55555 /* round to nearest */, 255.0);
        std::uint32_t green = std::min(value[width * height] / double(blue_max) * 2 * 255.0 + 0.55555 /* round to nearest */, 255.0);
        std::uint32_t blue = std::min(value[width * height * 2] / double(green_max) * 2 * 255.0 + 0.55555 /* round to nearest */, 255.0);
        return red | green << 8 | blue << 16 | 255u << 24;
    }
}

int main()
{
    auto t1 = std::chrono::high_resolution_clock::now();

    std::vector<heatmap_t> image(width * height * 3);
    heatmap_t red_max_value = 0, green_max_value = 0, blue_max_value = 0;

#pragma omp parallel
    generate(image, red_max_value, green_max_value, blue_max_value);

    /* output png image */
    std::vector<std::uint32_t> png_image(width * height);
    for (std::uint32_t i = 0; i < width * height; ++i)
    {
        png_image[i] = get_color(&image[i], red_max_value, green_max_value, blue_max_value);
    }
    png::write_image("result.png", png_image.data(), width, height);


    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    double timeTook = duration.count();
    std::cout << "It took " << timeTook
        << (timeTook == 1.0 ? " second" : "seconds") /* while it is very unlikly I just want to make sure */ << '\n';

    return 0;
}
