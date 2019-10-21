#include <math.h>
#include "image.h"

float nn_interpolate(image im, float x, float y, int c)
{
    return get_pixel(im, round(x), round(y), c);
}

image nn_resize(image im, int w, int h)
{
    // Create new image
    image new_image = make_image(w, h, im.c);

    // Loop pixels with a new coordinates
    int i, j, k;
    float x, y;
    float value;
    for (i=0; i < w; i++) {
        for (j=0; j < h; j++) {
            for (k=0; k < im.c; k++) {
                // Calculate value
                x = ((i + 0.5) * (im.w) / (float)(w)) - 0.5;
                y = ((j + 0.5) * (im.h) / (float)(h)) - 0.5;
                value = nn_interpolate(im, x, y, k);

                // Set pixel
                set_pixel(new_image, i, j, k, value);
            }
        }
    }
    return new_image;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // Calculate Nearby pixels
    int x_floor = floor(x);
    int x_ceil = ceil(x);
    int y_floor = floor(y);
    int y_ceil = ceil(y);

    // Get values for 4 pixels
    float v1 = get_pixel(im, x_floor, y_floor, c);
    float v2 = get_pixel(im, x_ceil, y_floor, c);
    float v3 = get_pixel(im, x_floor, y_ceil, c);
    float v4 = get_pixel(im, x_ceil, y_ceil, c);

    // Calculate distances
    float d1 = x - x_floor;
    float d2 = x_ceil - x;
    float d3 = y - y_floor;
    float d4 = y_ceil - y;

    // Interpolate
    float q1 = v1 * d2 + v2 * d1;
    float q2 = v3 * d2 + v4 * d1;
    float q = q1 * d4 + q2 * d3;

    return q;
}

image bilinear_resize(image im, int w, int h)
{
    // Create new image
    image new_image = make_image(w, h, im.c);

    // Loop pixels with a new coordinates
    int i, j, k;
    float x, y;
    float value;
    for (i=0; i < w; i++) {
        for (j=0; j < h; j++) {
            for (k=0; k < im.c; k++) {
                // Calculate value
                x = ((i + 0.5) * (im.w) / (float)(w)) - 0.5;
                y = ((j + 0.5) * (im.h) / (float)(h)) - 0.5;
                value = bilinear_interpolate(im, x, y, k);

                // Set pixel
                set_pixel(new_image, i, j, k, value);
            }
        }
    }
    return new_image;
}
