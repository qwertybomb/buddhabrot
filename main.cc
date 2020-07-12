#include <iostream>
#include <random>
#include <array>
#include <vector>
#include <cstdint>
#include <omp.h>

/* png headers */
#include <png.hpp>


typedef float v8sf __attribute__ ((vector_size(32)));
typedef int32_t v8si __attribute__ ((vector_size(32)));
typedef uint32_t heatmap_t;

constexpr float min_r = -2.0;
constexpr float max_r = 0.8;
constexpr float min_i = -1.414;
constexpr float max_i = 1.414;

constexpr v8sf min_r_v8sf = {min_r, min_r, min_r, min_r, min_r, min_r, min_r, min_r};
constexpr v8sf max_r_v8sf = {max_r, max_r, max_r, max_r, max_r, max_r, max_r, max_r};
constexpr v8sf min_i_v8sf = {min_i, min_i, min_i, min_i, min_i, min_i, min_i, min_i};
constexpr v8sf max_i_v8sf = {max_i, max_i, max_i, max_i, max_i, max_i, max_i, max_i};

constexpr int32_t width = 750;
constexpr int32_t height = 750;
constexpr int32_t max_iterations = 500;
constexpr uint64_t total_samples = 50000000;
constexpr uint64_t individual_samples = total_samples/64;

static void mandelbrot(v8sf cx,v8sf cy, std::array<v8sf[2], max_iterations>& orbit, v8si& iterations)
{
    auto good = [](v8sf cx, v8sf cy) {
        // test H1 and H2
        v8sf c2 = cx * cx + cy * cy;
        v8sf d[2] = {cx + 1, cx};
        v8sf h1 = 256 * c2 * c2 - 96 * c2 + 32 * cx - 3;
        v8sf h2 = 16 * (d[0]*d[0]+d[1]*d[1]) - 1.0;
        return h1 >= 0 && h2 >= 0;
    };
    v8sf zx = {0};
    v8sf zy = {0};
    v8sf zx_temp = {0};
    iterations = good(cx, cy) ? 0 : max_iterations;
    int32_t orbit_counter = 0;
    v8si active = {1,1,1,1,1,1,1,1};
    for(;iterations[0] < max_iterations &&
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
         || active[7] != 0);orbit_counter++, iterations+=-active)
    {
        zx_temp = active != 0 ? zx_temp*zx_temp-zy*zy+cx : zx_temp;
        zy = active != 0 ? zx * zy * 2 + cy : zy;
        zx = zx_temp;
        orbit[orbit_counter][0] = zx;
        orbit[orbit_counter][1] = zy;
        active = active !=0 ? (zx*zx*zy*zy < 4) : active;
    }
    /* if point escaped don't plot it */
    iterations = iterations >= max_iterations ? iterations^iterations : iterations;
}

static void generate(std::vector<heatmap_t>& image, heatmap_t& max_value)
{

    auto map = [](auto x, auto in_min, auto in_max, auto out_min, auto out_max)
    {
        return (x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
    };

    std::random_device rd;
    std::mt19937 engine {rd()};
    std::uniform_real_distribution<float> dist_r(min_r, max_r);
    std::uniform_real_distribution<float> dist_i(min_i, max_i);

    std::array<v8sf [2], max_iterations> orbit;
    v8si iterations;

    for(size_t i = 0; i < individual_samples; ++i)
    {
        v8sf cx = {dist_r(engine), dist_r(engine), dist_r(engine), dist_r(engine),
                   dist_r(engine), dist_r(engine), dist_r(engine), dist_r(engine)};
        v8sf cy = {dist_i(engine), dist_i(engine), dist_i(engine), dist_i(engine),
                   dist_i(engine), dist_i(engine), dist_i(engine), dist_i(engine)};
        mandelbrot(cx, cy, orbit, iterations);
        int32_t max_iteration = std::max({iterations[0], iterations[1], iterations[2], iterations[3],
                                          iterations[4], iterations[5], iterations[6], iterations[7]});
        /* add points in orbit to heatmap */
        for(int32_t i = 0; i < max_iteration;++i, iterations-=1)
        {
            /* map the point to pixel */
            v8sf row = map(orbit[i][0], min_r_v8sf, max_r_v8sf,
                    (v8sf){0, 0, 0, 0, 0, 0, 0, 0},
                    (v8sf){width -1 , width - 1 , width - 1, width -1,
                           width -1 , width - 1 , width - 1, width -1});
            v8sf col = map(orbit[i][1], min_i_v8sf, max_i_v8sf,
                           (v8sf){0, 0, 0, 0, 0, 0, 0, 0},
                           (v8sf){height -1 , height - 1 , height - 1, height -1,
                                  height -1 , height - 1 , height - 1, height -1});
            /* check if the point did escape */
            v8si out_of_bounds = row < 0 || row >= width || col < 0 || col >= height || iterations <= 0;
            heatmap_t pixel1 = !out_of_bounds[0] ? ++image[uint32_t(col[0])*width+row[0]] : 0;
            heatmap_t pixel2 = !out_of_bounds[1] ? ++image[uint32_t(col[1])*width+row[1]] : 0;
            heatmap_t pixel3 = !out_of_bounds[2] ? ++image[uint32_t(col[2])*width+row[2]] : 0;
            heatmap_t pixel4 = !out_of_bounds[3] ? ++image[uint32_t(col[3])*width+row[3]] : 0;
            heatmap_t pixel5 = !out_of_bounds[4] ? ++image[uint32_t(col[4])*width+row[4]] : 0;
            heatmap_t pixel6 = !out_of_bounds[5] ? ++image[uint32_t(col[5])*width+row[5]] : 0;
            heatmap_t pixel7 = !out_of_bounds[6] ? ++image[uint32_t(col[6])*width+row[6]] : 0;
            heatmap_t pixel8 = !out_of_bounds[7] ? ++image[uint32_t(col[7])*width+row[7]] : 0;
            max_value = std::max({max_value, pixel1, pixel2, pixel3, pixel4,
                                  pixel5, pixel6, pixel7, pixel8});
        }
    }

}

static png::rgb_pixel get_color(heatmap_t value, heatmap_t max_value)
{
    png::byte color =  std::min(value/double(max_value) * 2 * 255.0 + 0.55555 /* round to nearest */, 255.0);
    return png::rgb_pixel(color, color, color);
}

int main()
{
    std::vector<heatmap_t> image(width*height);
    heatmap_t max_value = 0;
    #pragma omp parallel
    generate(image, max_value);

    /* output png image */
    png::image<png::rgb_pixel> png_image(width, height);
    for(uint32_t i = 0; i < width; ++i)
    {
        for (uint32_t j = 0; j < height; ++j)
        {
            png_image.set_pixel(j,i, get_color(image[j * width + i], max_value));
        }
    }
    png_image.write("output.png");
    return 0;
}
