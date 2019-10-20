#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

// Macro for calculating CHW format index
#define chw(C, H, W, x, y, c)   (c * (W*H) + y * (W) + x)

float get_pixel(image im, int x, int y, int c)
{
    // Clamp padding strategy
    if (x < 0) x = 0;
    if (x >= im.w) x = im.w-1;
    if (y < 0) y = 0;
    if (y >= im.h) y = im.h-1;

    // Calculate index
    int idx = chw(im.c, im.h, im.w, x, y, c);
    return im.data[idx];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // Validate input
    if (x < 0 || x >= im.w) return;
    if (y < 0 || y >= im.h) return;
    if (c < 0 || c >= im.c) return;

    // Calculate index
    int idx = chw(im.c, im.h, im.w, x, y, c);
    im.data[idx] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    // Copy attributes
    copy.w = im.w;
    copy.h = im.h;
    copy.c = im.c;
    memcpy(copy.data, im.data, im.w * im.h * im.c * 4);
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    int i, j;
    float r, g, b;
    float grayness;
    for (i=0; i < im.w; i++) {
        for (j=0; j < im.h; j++) {
            r = get_pixel(im, i, j, 0);
            g = get_pixel(im, i, j, 1);
            b = get_pixel(im, i, j, 2);
            grayness = 0.299 * r + 0.587 * g + .114 * b;
            set_pixel(gray, i, j, 0, grayness);
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    int i, j;
    float value;
    // Iterate through width and height
    for (i=0; i < im.w; i++) {
        for (j=0; j < im.h; j++) {
            value = get_pixel(im, i, j, c) + v;
            set_pixel(im, i, j, c, value);
        }
    }
}

void clamp_image(image im)
{
    int i, j, k, idx;
    for (i=0; i < im.w; i++) {
        for (j=0; j < im.h; j++) {
            for (k=0; k < im.c; k++) {
                idx = chw(im.c, im.h, im.w, i, j, k);
                im.data[idx] = MAX(0.0, im.data[idx]);
                im.data[idx] = MIN(1.0, im.data[idx]);
            }
        }
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    int i, j;
    float R, G, B;
    float H, S, V;
    float m, C, h_prime;
    for (i=0; i < im.w; i++) {
        for (j=0; j < im.h; j++) {
            // Get RGB
            R = get_pixel(im, i, j, 0);
            G = get_pixel(im, i, j, 1);
            B = get_pixel(im, i, j, 2);

            // Calculate Value
            V = three_way_max(R, G, B);

            // Calculate Saturation
            m = three_way_min(R, G, B);
            C = V - m;

            // If RGB values are all 0, then Saturation is 0.
            // Otherwise Saturation = C/V
            S = (R == 0 && G == 0 && B == 0) ? (0.0) : (C/V);

            // Calculate h_prime
            if (V == R) {
                h_prime = (G-B)/C;
            } else if (V == G) {
                h_prime = ((B-R)/C) + 2.0;
            } else {
                h_prime = ((R-G)/C) + 4.0;
            }

            // Calculate Hue from h_prime
            H = (h_prime < 0) ? ((h_prime/6) + 1.0) : (h_prime/6);

            // Set Hue to 0 if C = 0
            if (C == 0) H = 0.0;

            // Set HSV
            set_pixel(im, i, j, 0, H);
            set_pixel(im, i, j, 1, S);
            set_pixel(im, i, j, 2, V);
        }
    }
}

// My magic helper function to calculate RGB values from HSV
float magic(int n, float h, float s, float v)
{
    float k = (n + h*6);
    if (k >= 6) k -= 6.0;
    return v - v * s * MAX(three_way_min(k, 4-k, 1), 0);
}

void hsv_to_rgb(image im)
{
    int i, j;
    float R, G, B;
    float H, S, V;
    for (i=0; i < im.w; i++) {
        for (j=0; j < im.h; j++) {
            // Get HSV
            H = get_pixel(im, i, j, 0);
            S = get_pixel(im, i, j, 1);
            V = get_pixel(im, i, j, 2);

            // Do some magic
            R = magic(5, H, S, V);
            G = magic(3, H, S, V);
            B = magic(1, H, S, V);

            // Set RGB
            set_pixel(im, i, j, 0, R);
            set_pixel(im, i, j, 1, G);
            set_pixel(im, i, j, 2, B);
        }
    }
}

void scale_image(image im, int c, float v)
{
    int i, j;
    float value;
    // Iterate through width and height
    for (i=0; i < im.w; i++) {
        for (j=0; j < im.h; j++) {
            value = get_pixel(im, i, j, c) * v;
            set_pixel(im, i, j, c, value);
        }
    }
}
