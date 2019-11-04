#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    // Calculate summation
    int i, j, k;
    float total;
    for (k=0; k < im.c; k++) {
        total = 0.0;
        // Get pixel values
        for (i=0; i < im.w; i++) {
            for (j=0; j < im.h; j++) {
                total += get_pixel(im, i, j, k);
            }
        }

        // Set pixel values
        for (i=0; i < im.w; i++) {
            for (j=0; j < im.h; j++) {
                set_pixel(im, i, j, k, get_pixel(im, i, j, k)/total);
            }
        }
    }
}

image make_box_filter(int w)
{
    // Create image
    image new_image = make_image(w, w, 1);

    // Set 1 for every pixel
    int i, j;
    for (i=0; i < new_image.w; i++) {
        for (j=0; j < new_image.h; j++) {
            set_pixel(new_image, i, j, 0, 1.0);
        }
    }

    // Normalize image
    l1_normalize(new_image);

    return new_image;
}

image convolve_image(image im, image filter, int preserve)
{
    // Assumptions:
    // This function assumes filter size is odd in terms of its width and height

    image new_image;
    int output_channels;
    int i, j, k;
    int di, dj;
    float value;
    float pointwise_total;
    float pixel, weight;
    int repeat = 0;

    // Sanity check
    assert(filter.c == 1 || filter.c == im.c);

    // Case 1 - Normal convolution where you apply 3D filter. This will produce a 1 channel image.
    // preserve = 0, repeat = 0

    // Case 2 - we want to apply different 1-channel 2D filter to all channels using separate weight values
    // A Depthwise Convolution. There will be `c` different filters.
    // preserve = 1, repeat = 0

    // Case 3 - we want to apply the same 1-channel 2D filter to all channels.
    // A Depthwise Convolution with a single shared filter that is applied to each of the channels.
    // preserve = 1, repeat = 1  (preserve each filtered channel)
    // preserve = 0, repeat = 1  (sum along the channel)

    // Flag setting
    output_channels = (preserve) ? (im.c) : (1);
    if (filter.c == 1) repeat = 1;

    // Create image
    new_image = make_image(im.w, im.h, output_channels);

    // Convolution
    for (i=0; i < im.w; i++) {
        for (j=0; j < im.h; j++) {
            pointwise_total = 0;
            for (k=0; k < im.c; k++) {
                // Convolve
                value = 0;
                for (di=0; di < filter.w; di++) {
                    for (dj=0; dj < filter.h; dj++) {
                        pixel = get_pixel(im, i+di-filter.w/2, j+dj-filter.h/2, k);
                        if (repeat) {
                            weight = get_pixel(filter, di, dj, 0);
                        } else {
                            weight = get_pixel(filter, di, dj, k);
                        }
                        value += pixel * weight;
                    }
                }
                if (preserve) set_pixel(new_image, i, j, k, value);
                pointwise_total += value;
            }
            if (!preserve) set_pixel(new_image, i, j, 0, pointwise_total);
        }
    }

    return new_image;
}

void helper(image filter, float values[]) {
    int v = 0;
    int i, j;
    for (i=0; i < filter.w; i++) {
        for (j=0; j < filter.h; j++) {
            set_pixel(filter, i, j, 0, values[v]);
            v += 1;
        }
    }
}

image make_highpass_filter()
{
    image filter = make_image(3,3,1);
    float values[] = {0, -1, 0, -1, 4, -1, 0, -1, 0};
    helper(filter, values);

    return filter;
}

image make_sharpen_filter()
{
    image filter = make_image(3,3,1);
    float values[] = {0, -1, 0, -1, 5, -1, 0, -1, 0};
    helper(filter, values);

    return filter;
}

image make_emboss_filter()
{
    image filter = make_image(3,3,1);
    float values[] = {-2, -1, 0, -1, 1, 1, 0, 1, 2};
    helper(filter, values);

    return filter;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
    // Next highest odd integer from 6x sigma
    int n = ceil(6*sigma);
    if (n % 2 == 0) n += 1;

    image filter = make_image(n, n, 1);
    int i, j;
    float x, y;
    float value;
    for (i=0; i < filter.w; i++) {
        for (j=0; j < filter.h; j++) {
            x = i-filter.w/2;
            y = j-filter.h/2;
            value = (1/(TWOPI*sigma*sigma)) * exp(-(x*x + y*y)/(2*sigma*sigma));
            set_pixel(filter, i, j, 0, value);
        }
    }

    l1_normalize(filter);
    return filter;
}

image add_image(image a, image b)
{
    // Sanity check
    assert(a.w == b.w);
    assert(a.h == b.h);
    assert(a.c == b.c);

    // Create new image
    image new_image = make_image(a.w, a.h, a.c);
    int i, j, k;
    float value;
    for (i=0; i < new_image.w; i++) {
        for (j=0; j < new_image.h; j++) {
            for (k=0; k < new_image.c; k++) {
                value = get_pixel(a, i, j, k) + get_pixel(b, i, j, k);
                set_pixel(new_image, i, j, k, value);
            }
        }
    }

    return new_image;
}

image sub_image(image a, image b)
{
    // Sanity check
    assert(a.w == b.w);
    assert(a.h == b.h);
    assert(a.c == b.c);

    // Create new image
    image new_image = make_image(a.w, a.h, a.c);
    int i, j, k;
    float value;
    for (i=0; i < new_image.w; i++) {
        for (j=0; j < new_image.h; j++) {
            for (k=0; k < new_image.c; k++) {
                value = get_pixel(a, i, j, k) - get_pixel(b, i, j, k);
                set_pixel(new_image, i, j, k, value);
            }
        }
    }

    return new_image;
}

image make_gx_filter()
{
    image filter = make_image(3,3,1);
    float values[] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
    helper(filter, values);

    return filter;
}

image make_gy_filter()
{
    image filter = make_image(3,3,1);
    float values[] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
    helper(filter, values);

    return filter;
}

void feature_normalize(image im)
{
    float min = im.data[0];
    float max = im.data[0];
    float range;

    int i, j, k;
    float value;
    for (i=0; i < im.w; i++) {
        for (j=0; j < im.h; j++) {
            for (k=0; k < im.c; k++) {
                // Update min, max
                value = get_pixel(im, i, j, k);
                if (min > value) min = value;
                if (max < value) max = value;
            }
        }
    }

    range = max - min;
    for (i=0; i < im.w; i++) {
        for (j=0; j < im.h; j++) {
            for (k=0; k < im.c; k++) {
                value = get_pixel(im, i, j, k);
                value = (range == 0) ? 0 : (value-min)/range;
                set_pixel(im, i, j, k, value);
            }
        }
    }
}

image *sobel_image(image im)
{
    // Generate filter
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();

    // Convolution
    image gx = convolve_image(im, gx_filter, 0);
    image gy = convolve_image(im, gy_filter, 0);

    // Create two images
    image *two_images = calloc(2, sizeof(image));
    two_images[0] = make_image(im.w, im.h, 1);
    two_images[1] = make_image(im.w, im.h, 1);

    int i, j, k;
    float value;

    // Gradient Magnitude
    float v1, v2;
    for (i=0; i < im.w; i++) {
        for (j=0; j < im.h; j++) {
            for (k=0; k < im.c; k++) {
                v1 = get_pixel(gx, i, j, k);
                v2 = get_pixel(gy, i, j, k);
                value = sqrt(v1*v1 + v2*v2);
                set_pixel(two_images[0], i, j, k, value);
            }
        }
    }

    // Gradient Direction
    for (i=0; i < im.w; i++) {
        for (j=0; j < im.h; j++) {
            v1 = get_pixel(gx, i, j, 0);
            v2 = get_pixel(gy, i, j, 0);
            value = atan2(v2, v1);
            set_pixel(two_images[1], i, j, 0, value);
        }
    }

    // Cleanup
    free_image(gx_filter);
    free_image(gy_filter);
    free_image(gx);
    free_image(gy);

    return two_images;
}

image colorize_sobel(image im)
{
    // Get magnitude and angle for the image
    image *images = sobel_image(im);
    feature_normalize(images[1]);

    // Create image
    image new_image = make_image(im.w, im.h, 3);

    float magnitude, angle;
    int i, j;
    for (i=0; i < im.w; i++) {
        for (j=0; j < im.h; j++) {
            // Get magnitude and angle
            magnitude = get_pixel(images[0], i, j, 0);
            angle = get_pixel(images[1], i, j, 0);

            // Set HSV
            set_pixel(new_image, i, j, 0, angle);
            set_pixel(new_image, i, j, 1, magnitude);
            set_pixel(new_image, i, j, 2, magnitude);
        }
    }

    // Convert hsv to rgb
    hsv_to_rgb(new_image);

    // Cleanup
    free_image(images[0]);
    free_image(images[1]);
    free(images);

    return new_image;
}
