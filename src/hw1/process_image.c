#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

// I IMPLEMENTED THE EXTRA CREDIT PART

static int clamp(int value, int min, int max) {
    if (value < min){
        return min;
    } else if (value > max) {
        return max;
    }
    return value;
}

static float clamp_float(float value, float min, float max) {
    if (value < min){
        return min;
    } else if (value > max) {
        return max;
    }
    return value;
}


float get_pixel(image im, int x, int y, int c)
{
    // TODO Fill this in
    x = clamp(x, 0, im.w -1);
    y = clamp(y, 0, im.h -1);
    c = clamp(c, 0, im.c -1);
    return im.data[x + y * im.w + c * im.h * im.w];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // TODO Fill this in
    im.data[x + y * im.w + c * im.h * im.w] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    // TODO Fill this in
    memcpy(copy.data, im.data, im.w * im.h * im.c * sizeof(float));
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            float r = get_pixel(im, i, j, 0);
            float g = get_pixel(im, i, j, 1);
            float b = get_pixel(im, i, j, 2);
            set_pixel(gray, i, j, 0, (0.299 * r + 0.587 * g + 0.114 * b));
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            float val = get_pixel(im, i, j, c);
            set_pixel(im, i, j, c, val + v);
        }
    }
}

void clamp_image(image im)
{
    // TODO Fill this in
    for (int i = 0; i < im.w * im.h * im.c; i++) {
        im.data[i] = clamp_float(im.data[i], 0.0, 1.0);
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
    // TODO Fill this in
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            float r = get_pixel(im, i, j, 0);
            float g = get_pixel(im, i, j, 1);
            float b = get_pixel(im, i, j, 2);

            float v = three_way_max(r, g, b);
            float m = three_way_min(r, g, b);
            float c = v - m;
            float s = (v == 0) ? 0 : c / v;

            set_pixel(im, i, j, 1, s);
            set_pixel(im, i, j, 2, v);

            float h;
            if (c == 0) {
                h = 0;
            } else if (v == r) {
                h = (g - b) / c;
            } else if (v == g) {
                h = (b - r) / c + 2;
            } else {
                h = (r - g) / c + 4;
            }

            if (h < 0) {
                h = h / 6 + 1;
            } else {
                h = h / 6;
            }
            set_pixel(im, i, j, 0, h);

        }
    }
}

void hsv_to_rgb(image im)
{
    // TODO Fill this in
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            float h = get_pixel(im, i, j, 0);
            float s = get_pixel(im, i, j, 1);
            float v = get_pixel(im, i, j, 2);

            h = 6 * h;
            float Hi = floor(h);
            float f = h - Hi;
            float p = v * (1 - s);
            float q = v * (1 - f * s);
            float t = v * (1 - (1 - f) * s);


            if (Hi == 0) {
                set_pixel(im, i, j, 0, v);
                set_pixel(im, i, j, 1, t);
                set_pixel(im, i, j, 2, p);
            } else if (Hi == 1) {
                set_pixel(im, i, j, 0, q);
                set_pixel(im, i, j, 1, v);
                set_pixel(im, i, j, 2, p);
            } else if (Hi == 2) {
                set_pixel(im, i, j, 0, p);
                set_pixel(im, i, j, 1, v);
                set_pixel(im, i, j, 2, t);
            } else if (Hi == 3) {
                set_pixel(im, i, j, 0, p);
                set_pixel(im, i, j, 1, q);
                set_pixel(im, i, j, 2, v);
            } else if (Hi == 4) {
                set_pixel(im, i, j, 0, t);
                set_pixel(im, i, j, 1, p);
                set_pixel(im, i, j, 2, v);
            } else if (Hi == 5) {
                set_pixel(im, i, j, 0, v);
                set_pixel(im, i, j, 1, p);
                set_pixel(im, i, j, 2, q);
            }
            
        }
    }
}

void scale_image(image im, int c, float v) {
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            float val = get_pixel(im, i, j, c);
            set_pixel(im, i, j, c, val * v);
        }
    }
}

void rgb_to_lch(image im) {
    // srgb to linear rgb
    for (int i = 0; i < im.w * im.h * im.c; i++) {
        if (im.data[i] <= 0.4045) {
            im.data[i] = im.data[i] / 12.92;
        } else {
            im.data[i] = powf((im.data[i] + 0.055)/1.055, 2.4);
        }
    }

    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            float r = get_pixel(im, i, j, 0);
            float g = get_pixel(im, i, j, 1);
            float b = get_pixel(im, i, j, 2);

            // rgb to xyz 
            float X = 0.4124 * r + 0.3576 * g + 0.1805 * b;
            float Y = 0.2126 * r + 0.7152 * g + 0.0722 * b;
            float Z = 0.0193 * r + 0.1192 * g + 0.9505 * b;

            // xyz to CIELUV
            float L, u , v;
            if (Y <= powf((6.0/29), 3)) {
                L = powf((29.0/3.0), 3) * Y;
            } else {
                L = 116 * powf(Y, 1.0/3.0) -16;
            }

            u = 13 * L * (((4 * X) / (X + 15 * Y + 3 * Z)) - 0.2009);
            v = 13 * L * (((9 * Y) / (X + 15 * Y + 3 * Z)) - 0.4610);


            set_pixel(im, i, j, 0, L);
            set_pixel(im, i, j, 1, sqrtf(powf(u, 2.0) + powf(v, 2.0)));
            set_pixel(im, i, j, 2, atan2f(v, u));
        }
    }
}