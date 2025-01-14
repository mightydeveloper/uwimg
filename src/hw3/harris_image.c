#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

// Macros
#define Square(x) ((x)*(x))

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// point p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    // TODO: optional, make separable 1d Gaussian.
    return make_image(1,1,1);
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    if(1){
        image g = make_gaussian_filter(sigma);
        image s = convolve_image(im, g, 1);
        free_image(g);
        return s;
    } else {
        // TODO: optional, use two convolutions with 1d gaussian filter.
        // If you implement, disable the above if check.
        return copy_image(im);
    }
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    // Generate filters
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();

    // Convolution
    image gx = convolve_image(im, gx_filter, 0);
    image gy = convolve_image(im, gy_filter, 0);
    free_image(gx_filter);
    free_image(gy_filter);

    // Calculate Naive Ix^2, Iy^2 and IxIy (before weighted sum)
    image S = make_image(im.w, im.h, 3);
    int i, j;
    float Ix, Iy;
    for (i=0; i < S.w; i++) {
        for (j=0; j < S.h; j++) {
            // Get Ix and Iy
            Ix = get_pixel(gx, i, j, 0);
            Iy = get_pixel(gy, i, j, 0);

            // Set values
            set_pixel(S, i, j, 0, Ix*Ix);
            set_pixel(S, i, j, 1, Iy*Iy);
            set_pixel(S, i, j, 2, Ix*Iy);
        }
    }

    // Calculate weighted sum of these Ix^2, Iy^2, and IxIy using Gaussian blur
    image output = smooth_image(S, sigma);

    // Cleanup
    free_image(gx);
    free_image(gy);
    free_image(S);

    return output;
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    // Hyper-params
    float alpha = 0.06;

    // Make a response map of cornerness calculations for each pixel
    image R = make_image(S.w, S.h, 1);
    int i, j;
    float Ix2, Iy2, IxIy;
    float value;
    for (i=0; i < R.w; i++) {
        for (j=0; j < R.h; j++) {
            // Get structure matrix information
            Ix2 = get_pixel(S, i, j, 0);
            Iy2 = get_pixel(S, i, j, 1);
            IxIy = get_pixel(S, i, j, 2);

            // Calculate our cornerness estimate value
            // "Smallest eigen value estimate" = det(S) - alpha * trace(S)^2, alpha = .06.
            // det(S) = a * d - b * c
            // trace(S) = a + d
            // for [[a b] [c d]] matrix
            value = (Ix2 * Iy2 - Square(IxIy)) - alpha * Square(Ix2 + Iy2);

            // Set the value to the output image map
            set_pixel(R, i, j, 0, value);
        }
    }

    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
    // Constant
    float low_number = -999999.0;

    // Make a new image with initial values having original responses
    image R = copy_image(im);
    // Perform NMS on the response map.
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])
    int i, j, di, dj;
    float self, neighbor;
    bool should_break;
    for (i=0; i < R.w; i++) {
        for (j=0; j < R.h; j++) {
            self = get_pixel(im, i, j, 0);
            should_break = false;
            // Look around neighbors
            for (di=-w; di < w; di++) {
                for (dj=-w; dj < w; dj++) {
                    neighbor = get_pixel(im, i+di, j+dj, 0);
                    // Neighbors are too rich compare to myself. I must be supressed! (Gentrification??? lol)
                    if (neighbor > self) {
                        set_pixel(R, i, j, 0, low_number);
                        // Early exit
                        should_break = true;
                        break;
                    }
                }
                if (should_break) break;
            }
        }
    }
    return R;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);

    // Count number of responses over threshold
    int i, j;
    int idx;
    float response_value;
    int count = 0;
    for (i=0; i < Rnms.w; i++) {
        for (j=0; j < Rnms.h; j++) {
            response_value = get_pixel(Rnms, i, j, 0);
            if (response_value > thresh) {
                count++;
            }
        }
    }

    *n = count; // <- set *n equal to number of corners in image.

    // Fill in array *d with descriptors of corners, use describe_index.
    descriptor *d = calloc(count, sizeof(descriptor));
    int count2 = 0;
    for (i=0; i < Rnms.w; i++) {
        for (j=0; j < Rnms.h; j++) {
            response_value = get_pixel(Rnms, i, j, 0);
            if (response_value > thresh) {
                idx = (im.w * j) + i;
                d[count2++] = describe_index(im, idx);
            }
        }
    }

    // Cleanups
    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}
