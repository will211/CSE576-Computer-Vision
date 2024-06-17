#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

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
// ponit p: spot to mark in the image.
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
    // TODO: make separable 1d Gaussian.

    int k = (int)ceilf(6 * sigma);
    if (k % 2 == 0) {
        k += 1;
    }

    image gauss = make_image(k, 1, 1);
    for (int i = 0; i < k; i++) {
        float expo = exp(-((i - k/2) * (i - k/2)) / (2 * sigma * sigma));
        set_pixel(gauss, i, 0, 0, 1 / (sqrt(2 * M_PI) * sigma) * expo);
    }
    l1_normalize(gauss);
    return gauss;
}

image make_1d_gaussian_transposed(float sigma)
{
    // TODO: make separable 1d Gaussian transposed.

    int k = (int)ceilf(6 * sigma);
    if (k % 2 == 0) {
        k += 1;
    }

    image gauss = make_image(1, k, 1);
    for (int i = 0; i < k; i++) {
        float expo = exp(-((i - k/2) * (i - k/2)) / (2 * sigma * sigma));
        set_pixel(gauss, 0, i, 0, 1 / (sqrt(2 * M_PI) * sigma) * expo);
    }
    l1_normalize(gauss);
    return gauss;
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    // TODO: use two convolutions with 1d gaussian filter.
    image gauss_1d = make_1d_gaussian(sigma);
    image horizontal = convolve_image(im, gauss_1d, 1);
    
    image gauss_1d_transposed =  make_1d_gaussian_transposed(sigma);

    image vertical = convolve_image(horizontal, gauss_1d_transposed, 1);
    free_image(gauss_1d);
    free_image(horizontal);
    free_image(gauss_1d_transposed);

    return vertical;
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    image S = make_image(im.w, im.h, 3);
    // TODO: calculate structure matrix for im.

    image gx = make_gx_filter();
    image gy = make_gy_filter();
    image Ix = convolve_image(im, gx, 0);
    image Iy = convolve_image(im, gy, 0);

    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            set_pixel(S, i, j, 0, pow(get_pixel(Ix, i, j, 0), 2));
            set_pixel(S, i, j, 1, pow(get_pixel(Iy, i, j, 0), 2));
            set_pixel(S, i, j, 2, get_pixel(Ix, i, j, 0) * get_pixel(Iy, i, j, 0));
        }
    }
    free_image(gx);
    free_image(gy);
    free_image(Ix);
    free_image(Iy);

    image gaussian_filter = make_gaussian_filter(sigma);
    return convolve_image(S, gaussian_filter, 1);
    // return smooth_image(S, sigma);
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.
    float det = 0, trace = 0;
    float alpha = 0.06;

    for (int i = 0; i < S.w; i++) {
        for (int j = 0; j < S.h; j++) {
            float Ix = get_pixel(S, i, j, 0);
            float Iy = get_pixel(S, i, j, 1);
            float IxIy = get_pixel(S, i, j, 2);

            det = Ix * Iy - powf(IxIy, 2);
            trace = Ix + Iy;
            set_pixel(R, i, j, 0, det - alpha * trace * trace);
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
    image r = copy_image(im);
    // TODO: perform NMS on the response map.
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])
    float LOW = -999999;

    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j< im.h; j++) {
            int found = 0;
            float curr = get_pixel(im, i, j, 0);

            for (int x = -w; x <= w; x++) {
                for (int y = -w; y <= w; y++) {
                    int x1 = i + x;
                    int y1 = j + y;

                    if (x1 >= 0 && y1 >=0 && x1 < im.w && y1 < im.h) {
                        if (get_pixel(im, x1, y1, 0) > curr) {
                            set_pixel(r, i, j, 0, LOW);
                            found = 1;
                            break;
                        }
                    }
                }
                if (found) break;
            }
        }
    }

    return r;
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


    //TODO: count number of responses over threshold
    int count = 0; // change this
    for (int i = 0; i < Rnms.w; i++) {
        for (int j = 0; j < Rnms.h; j++) {
            if (get_pixel(Rnms, i, j, 0) > thresh) {
                count++;
            }
        }
    }

    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));
    //TODO: fill in array *d with descriptors of corners, use describe_index.
    int i = 0;
    for (int j = 0; j < Rnms.w * Rnms.h; j++) {
        if (Rnms.data[j] > thresh) {
            d[i] = describe_index(im, j);
            i++;
        }
    }

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
