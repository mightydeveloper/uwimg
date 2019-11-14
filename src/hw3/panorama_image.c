#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"

// Macros
#define Square(x) ((x)*(x))

// Comparator for matches
// const void *a, *b: pointers to the matches to compare.
// returns: result of comparison, 0 if same, 1 if a > b, -1 if a < b.
int match_compare(const void *a, const void *b)
{
    match *ra = (match *)a;
    match *rb = (match *)b;
    if (ra->distance < rb->distance) return -1;
    else if (ra->distance > rb->distance) return  1;
    else return 0;
}

// Helper function to create 2d points.
// float x, y: coordinates of point.
// returns: the point.
point make_point(float x, float y)
{
    point p;
    p.x = x; p.y = y;
    return p;
}

// Place two images side by side on canvas, for drawing matching pixels.
// image a, b: images to place.
// returns: image with both a and b side-by-side.
image both_images(image a, image b)
{
    image both = make_image(a.w + b.w, a.h > b.h ? a.h : b.h, a.c > b.c ? a.c : b.c);
    int i,j,k;
    for(k = 0; k < a.c; ++k){
        for(j = 0; j < a.h; ++j){
            for(i = 0; i < a.w; ++i){
                set_pixel(both, i, j, k, get_pixel(a, i, j, k));
            }
        }
    }
    for(k = 0; k < b.c; ++k){
        for(j = 0; j < b.h; ++j){
            for(i = 0; i < b.w; ++i){
                set_pixel(both, i+a.w, j, k, get_pixel(b, i, j, k));
            }
        }
    }
    return both;
}

// Draws lines between matching pixels in two images.
// image a, b: two images that have matches.
// match *matches: array of matches between a and b.
// int n: number of matches.
// int inliers: number of inliers at beginning of matches, drawn in green.
// returns: image with matches drawn between a and b on same canvas.
image draw_matches(image a, image b, match *matches, int n, int inliers)
{
    image both = both_images(a, b);
    int i,j;
    for(i = 0; i < n; ++i){
        int bx = matches[i].p.x;
        int ex = matches[i].q.x;
        int by = matches[i].p.y;
        int ey = matches[i].q.y;
        for(j = bx; j < ex + a.w; ++j){
            int r = (float)(j-bx)/(ex+a.w - bx)*(ey - by) + by;
            set_pixel(both, j, r, 0, i<inliers?0:1);
            set_pixel(both, j, r, 1, i<inliers?1:0);
            set_pixel(both, j, r, 2, 0);
        }
    }
    return both;
}

// Draw the matches with inliers in green between two images.
// image a, b: two images to match.
// matches *
image draw_inliers(image a, image b, matrix H, match *m, int n, float thresh)
{
    int inliers = model_inliers(H, m, n, thresh);
    image lines = draw_matches(a, b, m, n, inliers);
    return lines;
}

// Find corners, match them, and draw them between two images.
// image a, b: images to match.
// float sigma: gaussian for harris corner detector. Typical: 2
// float thresh: threshold for corner/no corner. Typical: 1-5
// int nms: window to perform nms on. Typical: 3
image find_and_draw_matches(image a, image b, float sigma, float thresh, int nms)
{
    int an = 0;
    int bn = 0;
    int mn = 0;
    descriptor *ad = harris_corner_detector(a, sigma, thresh, nms, &an);
    descriptor *bd = harris_corner_detector(b, sigma, thresh, nms, &bn);
    match *m = match_descriptors(ad, an, bd, bn, &mn);

    mark_corners(a, ad, an);
    mark_corners(b, bd, bn);
    image lines = draw_matches(a, b, m, mn, 0);

    free_descriptors(ad, an);
    free_descriptors(bd, bn);
    free(m);
    return lines;
}

// Calculates L1 distance between to floating point arrays.
// float *a, *b: arrays to compare.
// int n: number of values in each array.
// returns: l1 distance between arrays (sum of absolute differences).
float l1_distance(float *a, float *b, int n)
{
    float total = 0;
    int i;
    for (i=0; i < n; i++) {
        total += fabsf(a[i] - b[i]);
    }

    return total;
}

// Finds best matches between descriptors of two images.
// descriptor *a, *b: array of descriptors for pixels in two images.
// int an, bn: number of descriptors in arrays a and b.
// int *mn: pointer to number of matches found, to be filled in by function.
// returns: best matches found. each descriptor in a should match with at most
//          one other descriptor in b.
match *match_descriptors(descriptor *a, int an, descriptor *b, int bn, int *mn)
{
    int i,j;
    int bind;
    float lowest_distance;
    float distance;

    // We will have at most an matches.
    *mn = an;
    match *m = calloc(an, sizeof(match));
    for(j = 0; j < an; ++j){
        // For every descriptor in a, find best match in b.
        // record ai as the index in *a and bi as the index in *b.
        lowest_distance = INFINITY;
        bind = 0;
        for (i = 0; i < bn; ++i) {
            distance = l1_distance(a[j].data, b[i].data, a[j].n);
            if (distance < lowest_distance) {
                lowest_distance = distance;
                bind = i;
            }
        }

        m[j].ai = j;
        m[j].bi = bind; // <- should be index in b.
        m[j].p = a[j].p;
        m[j].q = b[bind].p;
        m[j].distance = lowest_distance; // <- should be the smallest L1 distance!
    }

    // We want matches to be injective (one-to-one).
    // Sort matches based on distance using match_compare and qsort.
    qsort(m, an, sizeof(match), match_compare);

    // Then throw out matches to the same element in b. Use seen to keep track.
    int *seen = calloc(bn, sizeof(int));
    int count = 0;
    for (j = 0; j < an; ++j) {
        if (seen[m[j].bi]) {
            continue; // Throw out elements we have already seen
        }
        m[count++] = m[j]; // Shift other elements forward
        seen[m[j].bi] = 1;
    }

    // Set number of matches
    *mn = count;

    // Clean up
    free(seen);

    return m;
}

// Apply a projective transformation to a point.
// matrix H: homography to project point.
// point p: point to project.
// returns: point projected using the homography.
point project_point(matrix H, point p)
{
    matrix c = make_matrix(3, 1);
    c.data[0][0] = p.x;
    c.data[1][0] = p.y;
    c.data[2][0] = 1;
    // Project point p with homography H.
    matrix mat = matrix_mult_matrix(H, c);

    // Remember that homogeneous coordinates are equivalent up to scalar.
    // Have to divide by.... something...
    float w = mat.data[2][0];
    point q = make_point(mat.data[0][0]/w, mat.data[1][0]/w);
    return q;
}

// Calculate L2 distance between two points.
// point p, q: points.
// returns: L2 distance between them.
float point_distance(point p, point q)
{
    return sqrt(Square(p.x - q.x) + Square(p.y - q.y));
}

// Count number of inliers in a set of matches. Should also bring inliers
// to the front of the array.
// matrix H: homography between coordinate systems.
// match *m: matches to compute inlier/outlier.
// int n: number of matches in m.
// float thresh: threshold to be an inlier.
// returns: number of inliers whose projected point falls within thresh of
//          their match in the other image. Should also rearrange matches
//          so that the inliers are first in the array. For drawing.
int model_inliers(matrix H, match *m, int n, float thresh)
{
    int i;
    int count = 0;
    float distance;

    // Count number of matches that are inliers
    // i.e. distance(H*p, q) < thresh
    // Also, sort the matches m so the inliers are the first 'count' elements.
    printf("Calculating inliers\n");
    for (i=0; i < n; i++){
        distance = point_distance(project_point(H, m[i].p), m[i].q);
        if (distance < thresh) {
            m[count++] = m[i];
        }
    }

    return count;
}

// Randomly shuffle matches for RANSAC.
// match *m: matches to shuffle in place.
// int n: number of elements in matches.
void randomize_matches(match *m, int n)
{
    // Fisher-Yates to shuffle the array.
    printf("Shufflin\n");
    int i, j;
    match tmp;
    for (i=n-1; i > 0; i--) {
        j = rand() % (i+1); // Random integer such that 0 <= j <= i
        // Exchange m[j] and m[i]
        tmp = m[j];
        m[j] = m[i];
        m[i] = tmp;
    }
}

// Computes homography between two images given matching pixels.
// match *matches: matching points between images.
// int n: number of matches to use in calculating homography.
// returns: matrix representing homography H that maps image a to image b.
matrix compute_homography(match *matches, int n)
{
    printf("Computing homography\n");
    printf("argument n : %d\n", n);
    matrix M = make_matrix(n*2, 8);
    matrix b = make_matrix(n*2, 1);
    printf("After make matrix\n");

    int i;
    double m_x;
    double m_y;
    double n_x;
    double n_y;

    printf("Filling equations...");
    for(i = 0; i < n; ++i){
        // Get match point values
        m_x = matches[i].p.x;
        m_y = matches[i].p.y;
        n_x = matches[i].q.x;
        n_y = matches[i].q.y;

        // Filling in M
        M.data[2*i][0] = m_x;
        M.data[2*i][1] = m_y;
        M.data[2*i][2] = 1;
        M.data[2*i][3] = 0;
        M.data[2*i][4] = 0;
        M.data[2*i][5] = 0;
        M.data[2*i][6] = -m_x * n_x;
        M.data[2*i][7] = -m_y * n_x;

        M.data[2*i+1][0] = 0;
        M.data[2*i+1][1] = 0;
        M.data[2*i+1][2] = 0;
        M.data[2*i+1][3] = m_x;
        M.data[2*i+1][4] = m_y;
        M.data[2*i+1][5] = 1;
        M.data[2*i+1][6] = -m_x * n_y;
        M.data[2*i+1][7] = -m_y * n_y;

        // Filling in b
        b.data[2*i][0] = n_x;
        b.data[2*i+1][0] = n_y;
    }
    printf("Solving system...\n");
    matrix a = solve_system(M, b);
    //printf("!!!!! Matrix a has rows: %d and cols: %d\n", a.rows, a.cols);
    free_matrix(M); free_matrix(b);

    // If a solution can't be found, return empty matrix;
    matrix none = {0};
    if(!a.data) {
        printf("Solution not found! returning empty matrix!\n");
        return none;
    }

    matrix H = make_matrix(3, 3);
    // Fill in the homography H based on the result in a.
    printf("Filling in the homography based on the result in a\n");
    H.data[0][0] = a.data[0][0];
    H.data[0][1] = a.data[1][0];
    H.data[0][2] = a.data[2][0];

    H.data[1][0] = a.data[3][0];
    H.data[1][1] = a.data[4][0];
    H.data[1][2] = a.data[5][0];

    H.data[2][0] = a.data[6][0];
    H.data[2][1] = a.data[7][0];
    H.data[2][2] = 1;

    free_matrix(a);
    printf("Returning final Homography matrix\n");
    return H;
}

// Perform RANdom SAmple Consensus to calculate homography for noisy matches.
// match *m: set of matches.
// int n: number of matches.
// float thresh: inlier/outlier distance threshold.
// int k: number of iterations to run.
// int cutoff: inlier cutoff to exit early.
// returns: matrix representing most common homography between matches.
matrix RANSAC(match *m, int n, float thresh, int k, int cutoff)
{
    int i;
    int best = 0;
    int number_of_inliers;
    matrix Hb = make_translation_homography(256, 0);
    matrix H;

    // RANSAC algorithm.
    printf("Starting RANSAC algorithm\n");
    printf("We have %d matches\n", n);
    for (i=0; i < k; i++){
        randomize_matches(m, n);
        H = compute_homography(m, 4); // Fit model to sample
        // If solving failed, move on
        if (!H.rows) {
            printf("Empty matrix found!\n");
            continue;
        }
        // Count the number of inliers with modeled homography
        number_of_inliers = model_inliers(H, m, n, thresh);
        printf("We got %d inliers! \n", number_of_inliers);
        if (number_of_inliers > best) {
            H = compute_homography(m, number_of_inliers); // Fit model to all inliers
            // Update best values
            Hb = H;
            best = number_of_inliers;
            if (number_of_inliers > cutoff) {
                printf("Returning H b/c We got %d inliers! \n", number_of_inliers);
                printf("Cutoff was %d\n", cutoff);
                return H;
            }
        }
    }

    // Psudocode
    // for k iterations:
    //     shuffle the matches
    //     compute a homography with a few matches (how many??)
    //     if new homography is better than old (how can you tell?):
    //         compute updated homography using all inliers
    //         remember it and how good it is
    //         if it's better than the cutoff:
    //             return it immediately
    // if we get to the end return the best homography
    return Hb;
}

// Stitches two images together using a projective transformation.
// image a, b: images to stitch.
// matrix H: homography from image a coordinates to image b coordinates.
// returns: combined image stitched together.
image combine_images(image a, image b, matrix H)
{
    matrix Hinv = matrix_invert(H);

    // Project the corners of image b into image a coordinates.
    point c1 = project_point(Hinv, make_point(0,0));
    point c2 = project_point(Hinv, make_point(b.w-1, 0));
    point c3 = project_point(Hinv, make_point(0, b.h-1));
    point c4 = project_point(Hinv, make_point(b.w-1, b.h-1));

    // Find top left and bottom right corners of image b warped into image a.
    point topleft, botright;
    botright.x = MAX(c1.x, MAX(c2.x, MAX(c3.x, c4.x)));
    botright.y = MAX(c1.y, MAX(c2.y, MAX(c3.y, c4.y)));
    topleft.x = MIN(c1.x, MIN(c2.x, MIN(c3.x, c4.x)));
    topleft.y = MIN(c1.y, MIN(c2.y, MIN(c3.y, c4.y)));

    // Find how big our new image should be and the offsets from image a.
    int dx = MIN(0, topleft.x);
    int dy = MIN(0, topleft.y);
    int w = MAX(a.w, botright.x) - dx;
    int h = MAX(a.h, botright.y) - dy;

    // Can disable this if you are making very big panoramas.
    // Usually this means there was an error in calculating H.
    if(w > 7000 || h > 7000){
        fprintf(stderr, "output too big, stopping\n");
        return copy_image(a);
    }

    int i,j,k;
    float value;
    point p;
    image c = make_image(w, h, a.c);

    printf("Pasting image a into c!\n");
    printf("dx and dy are : %d,  %d\n", dx, dy);
    printf("w and h are: %d, %d\n", w, h);
    // Paste image a into the new image(c) offset by dx and dy.
    for(k = 0; k < a.c; ++k){
        for(j = 0; j < a.h; ++j){
            for(i = 0; i < a.w; ++i){
                value = get_pixel(a, i, j, k);
                set_pixel(c, i-dx, j-dy, k, value);
            }
        }
    }

    // Paste in image b to c
    // You should loop over some points in the new image (which? all?)
    // and see if their projection from a coordinates to b coordinates falls
    // inside of the bounds of image b. If so, use bilinear interpolation to
    // estimate the value of b at that projection, then fill in image c.
    printf("Pasting image b into c!\n");
    for (k = 0; k < c.c; ++k) {
        for (j = topleft.y; j < botright.y; ++j) {
            for (i = topleft.x; i < botright.x; ++i) {
                p = project_point(H, make_point(i, j));
                if (0 <= p.x && p.x <= (b.w-1) && 0 <= p.y && p.y <= (b.h-1)) {
                    value = bilinear_interpolate(b, p.x, p.y, k);
                    set_pixel(c, i-dx, j-dy, k, value);
                }
            }
        }
    }

    return c;
}

// Create a panoramam between two images.
// image a, b: images to stitch together.
// float sigma: gaussian for harris corner detector. Typical: 2
// float thresh: threshold for corner/no corner. Typical: 1-5
// int nms: window to perform nms on. Typical: 3
// float inlier_thresh: threshold for RANSAC inliers. Typical: 2-5
// int iters: number of RANSAC iterations. Typical: 1,000-50,000
// int cutoff: RANSAC inlier cutoff. Typical: 10-100
image panorama_image(image a, image b, float sigma, float thresh, int nms, float inlier_thresh, int iters, int cutoff)
{
    srand(10);
    int an = 0;
    int bn = 0;
    int mn = 0;

    // Calculate corners and descriptors
    descriptor *ad = harris_corner_detector(a, sigma, thresh, nms, &an);
    descriptor *bd = harris_corner_detector(b, sigma, thresh, nms, &bn);

    // Find matches
    match *m = match_descriptors(ad, an, bd, bn, &mn);

    // Run RANSAC to find the homography
    matrix H = RANSAC(m, mn, inlier_thresh, iters, cutoff);

    if(1){
        // Mark corners and matches between images
        mark_corners(a, ad, an);
        mark_corners(b, bd, bn);
        image inlier_matches = draw_inliers(a, b, H, m, mn, inlier_thresh);
        save_image(inlier_matches, "inliers");
    }

    free_descriptors(ad, an);
    free_descriptors(bd, bn);
    free(m);

    // Stitch the images together with the homography
    image comb = combine_images(a, b, H);
    return comb;
}

// Project an image onto a cylinder.
// image im: image to project.
// float f: focal length used to take image (in pixels).
// returns: image projected onto cylinder, then flattened.
image cylindrical_project(image im, float f)
{
    //TODO: project image onto a cylinder
    image c = copy_image(im);
    return c;
}
