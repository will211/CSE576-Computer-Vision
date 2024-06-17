#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

int cmpfunc(const void *a, const void *b);

/******************************** Resizing *****************************
  To resize we'll need some interpolation methods and a function to create
  a new image and fill it in with our interpolation methods.
************************************************************************/

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs nearest-neighbor interpolation on image "im"
      given a floating column value "x", row value "y" and integer channel "c",
      and returns the interpolated value.
    ************************************************************************/
    int int_x = roundf(x);
    int int_y = roundf(y);

    int_x = MAX(0, MIN(int_x, im.w - 1));
    int_y = MAX(0, MIN(int_y, im.h - 1));

    return get_pixel(im, int_x, int_y, c);
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix the return line)
    /***********************************************************************
      This function uses nearest-neighbor interpolation on image "im" to a new
      image of size "w x h"
    ************************************************************************/
    
    image resize = make_image(w, h, im.c);

    float ratio_w = (float) im.w / w;
    float ratio_h = (float) im.h / h;

    for (int i = 0; i < w; i++) {
      for (int j = 0; j < h; j++) {
        for (int c = 0; c < im.c; c++) {
          float ori_x = (i + 0.5) * ratio_w - 0.5;
          float ori_y = (j + 0.5) * ratio_h - 0.5;

          set_pixel(resize, i, j, c, nn_interpolate(im, ori_x, ori_y, c));
        }
      }
    }

    return resize;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs bilinear interpolation on image "im" given
      a floating column value "x", row value "y" and integer channel "c".
      It interpolates and returns the interpolated value.
    ************************************************************************/

    int x1 = floorf(x);
    int y1 = floorf(y);
    int x2 = ceilf(x);
    int y2 = ceilf(y);

    float t1 = get_pixel(im, x1, y1, c);
    float t2 = get_pixel(im, x2, y1, c);
    float b1 = get_pixel(im, x1, y2, c);
    float b2 = get_pixel(im, x2, y2, c);

    float dx1 = x - x1;
    float dy1 = y - y1;
    float dx2 = 1 - dx1;
    float dy2 = 1 - dy1;

    return t1 * dx2 * dy2 + t2 * dx1 * dy2
          + b1 * dx2 * dy1 + b2 * dx1 * dy1;

}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    /***********************************************************************
      This function uses bilinear interpolation on image "im" to a new image
      of size "w x h". Algorithm is same as nearest-neighbor interpolation.
    ************************************************************************/
    image resize = make_image(w, h, im.c);

    float ratio_w = (float) im.w / w;
    float ratio_h = (float) im.h / h;

    for (int i = 0; i < w; i++) {
      for (int j = 0; j < h; j++) {
        for (int c = 0; c < im.c; c++) {
          float ori_x = (i + 0.5) * ratio_w - 0.5;
          float ori_y = (j + 0.5) * ratio_h - 0.5;

          set_pixel(resize, i, j, c, bilinear_interpolate(im, ori_x, ori_y, c));
        }
      }
    }

    return resize;
}


/********************** Filtering: Box filter ***************************
  We want to create a box filter. We will only use square box filters.
************************************************************************/

void l1_normalize(image im)
{
    // TODO
    /***********************************************************************
      This function divides each value in image "im" by the sum of all the
      values in the image and modifies the image in place.
    ************************************************************************/
  float sum = 0;
  
  for (int i = 0; i < im.w; i++) {
    for (int j = 0; j < im.h; j++) {
      for (int c = 0; c < im.c; c++) {
        sum += get_pixel(im, i, j, c);
      }
    }
  }

  for (int i = 0; i < im.w; i++) {
    for (int j = 0; j < im.h; j++) {
      for (int c = 0; c < im.c; c++) {
        set_pixel(im, i, j, c, get_pixel(im, i, j ,c) / sum);
      }
    }
  }

}

image make_box_filter(int w)
{
    // TODO
    /***********************************************************************
      This function makes a square filter of size "w x w". Make an image of
      width = height = w and number of channels = 1, with all entries equal
      to 1. Then use "l1_normalize" to normalize your filter.
    ************************************************************************/
    image im = make_image(w, w, 1);
    for (int i = 0; i < w * w; i++) {
      im.data[i] = 1.0 / (w * w);
    }

    return im;
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    /***********************************************************************
      This function convolves the image "im" with the "filter". The value
      of preserve is 1 if the number of input image channels need to be 
      preserved. Check the detailed algorithm given in the README.  
    ************************************************************************/
    assert(filter.c == 1 || filter.c == im.c);
    image result;

    if (preserve) {
      result = make_image(im.w, im.h, im.c);
    } else {
      result = make_image(im.w, im.h, 1);
    }

    for (int i = 0; i < im.w; i++) {
      for (int j = 0; j < im.h; j++) {
        float convolve_val_1 = 0;
        for (int c = 0; c < im.c; c++) {
          float convolve_val = 0;

          for (int fw = 0; fw < filter.w; fw++) {
            for (int fh = 0; fh < filter.h; fh++) {
              int x = i + fw - (filter.w / 2);
              int y = j + fh - (filter.h / 2);
              convolve_val += get_pixel(im, x, y, c) * get_pixel(filter, fw, fh, 0);
            }
          }
          if (preserve) {
            set_pixel(result, i, j, c, convolve_val);
          } else {
            convolve_val_1 += convolve_val;
          }
        }

        if (preserve == 0) {
          set_pixel(result, i, j ,0, convolve_val_1);
        }
      }
    }

  return result;
}

image make_highpass_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with highpass filter values using image.data[]
    ************************************************************************/
    image im = make_image(3, 3, 1);
    im.data[0] = 0;
    im.data[1] = -1;
    im.data[2] = 0;
    im.data[3] = -1;
    im.data[4] = 4;
    im.data[5] = -1;
    im.data[6] = 0;
    im.data[7] = -1;
    im.data[8] = 0;
    
    return im;
}

image make_sharpen_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with sharpen filter values using image.data[]
    ************************************************************************/
    image im = make_image(3, 3, 1);
    im.data[0] = 0;
    im.data[1] = -1;
    im.data[2] = 0;
    im.data[3] = -1;
    im.data[4] = 5;
    im.data[5] = -1;
    im.data[6] = 0;
    im.data[7] = -1;
    im.data[8] = 0;
    
    return im;
}

image make_emboss_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with emboss filter values using image.data[]
    ************************************************************************/
    image im = make_image(3, 3, 1);
    im.data[0] = -2;
    im.data[1] = -1;
    im.data[2] = 0;
    im.data[3] = -1;
    im.data[4] = 1;
    im.data[5] = 1;
    im.data[6] = 0;
    im.data[7] = 1;
    im.data[8] = 2;
    
    return im;
}

// Question 2.3.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO
// We should use the sharpen and highpass filters since the sharpen can make the image appear more detailed and sharper.
// While the highpass is good for detecting edges.
// Do not use the emboss filter, since this filter is good for detecting texture.

// Question 2.3.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO
// We should do the post-processing on the highpass filter, since this filter has gray background with edges highlighted, 
// we will adjust the brightness and the contract of the image to make it more visible.

image make_gaussian_filter(float sigma)
{
    // TODO
    /***********************************************************************
      sigma: a float number for the Gaussian.
      Create a Gaussian filter with the given sigma. Note that the kernel size 
      is the next highest odd integer from 6 x sigma. Return the Gaussian filter.
    ************************************************************************/
    int k = (int)ceilf(6 * sigma);
    if (k % 2 == 0) {
      k += 1;
    }

    image im = make_image(k, k, 1);
    for (int i = 0; i < k; i++) {
      for (int j = 0; j < k; j++) {
        float expo = exp(-((i - k/2) * (i - k/2) + (j - k/2) * (j - k/2)) / (2 * sigma * sigma));
        set_pixel(im, i, j, 0, 1 / (2 * M_PI * sigma * sigma) * expo);
      }
    }
    l1_normalize(im);
    return im;
}

image add_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input images a and image b have the same height, width, and channels.
      Sum the given two images and return the result, which should also have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    image im = make_image(a.w, a.h, a.c);

    for (int i = 0; i < a.w; i++) {
      for (int j = 0; j < a.h; j++) {
        for (int c = 0; c < a.c; c++) {
          set_pixel(im, i, j, c, get_pixel(a, i, j, c) + get_pixel(b, i, j, c));
        }
      }
    }

    return im;
}

image sub_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input image a and image b have the same height, width, and channels.
      Subtract the given two images and return the result, which should have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    image im = make_image(a.w, a.h, a.c);

    for (int i = 0; i < a.w; i++) {
      for (int j = 0; j < a.h; j++) {
        for (int c = 0; c < a.c; c++) {
          set_pixel(im, i, j, c, get_pixel(a, i, j, c) - get_pixel(b, i, j, c));
        }
      }
    }

    return im;
}

image make_gx_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gx filter and return it
    ************************************************************************/
    image im = make_image(3, 3, 1);
    im.data[0] = -1;
    im.data[1] = 0;
    im.data[2] = 1;
    im.data[3] = -2;
    im.data[4] = 0;
    im.data[5] = 2;
    im.data[6] = -1;
    im.data[7] = 0;
    im.data[8] = 1;
    
    return im;
}

image make_gy_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gy filter and return it
    ************************************************************************/
    image im = make_image(3, 3, 1);
    im.data[0] = -1;
    im.data[1] = -2;
    im.data[2] = -1;
    im.data[3] = 0;
    im.data[4] = 0;
    im.data[5] = 0;
    im.data[6] = 1;
    im.data[7] = 2;
    im.data[8] = 1;
    
    return im;
}

void feature_normalize(image im)
{
    // TODO
    /***********************************************************************
      Calculate minimum and maximum pixel values. Normalize the image by
      subtracting the minimum and dividing by the max-min difference.
    ************************************************************************/

    float min_val = MAXFLOAT;
    float max_val = -MAXFLOAT;

    for (int i = 0; i < im.w * im.h * im.c; i++) {
      if (im.data[i] < min_val) {
        min_val = im.data[i];
      }

      if (im.data[i] > max_val) {
        max_val = im.data[i];
      }
    }

    if (max_val == min_val) {
      for (int i = 0; i < im.w * im.h * im.c; i++) {
        im.data[i] = 0;
      }
    } else {
      float different = max_val - min_val;
      for (int i = 0; i < im.w * im.h * im.c; i++) {
        im.data[i] = (im.data[i] - min_val) / different;
      } 
    }
}

image *sobel_image(image im)
{
    // TODO
    /***********************************************************************
      Apply Sobel filter to the input image "im", get the magnitude as sobelimg[0]
      and gradient as sobelimg[1], and return the result.
    ************************************************************************/
    image *sobelimg = calloc(2, sizeof(image));

    image gx = convolve_image(im, make_gx_filter(), 0);
    image gy = convolve_image(im, make_gy_filter(), 0);

    sobelimg[0] = make_image(im.w, im.h, 1);
    sobelimg[1] = make_image(im.w, im.h, 1);

    for (int i = 0; i < im.w; i++) {
      for (int j = 0; j < im.h; j++) {
        float gx_val = get_pixel(gx, i, j, 0);
        float gy_val = get_pixel(gy, i, j, 0);
        set_pixel(sobelimg[0], i, j, 0, sqrt(gx_val * gx_val + gy_val * gy_val));
        set_pixel(sobelimg[1], i, j, 0, atan2(gy_val, gx_val));
      }
    }

    free_image(gx);
    free_image(gy);

    return sobelimg;
}

image colorize_sobel(image im)
{
  // TODO
  /***********************************************************************
    Create a colorized version of the edges in image "im" using the 
    algorithm described in the README.
  ************************************************************************/
  image *sobelimg = sobel_image(im);
  image mag = sobelimg[0];
  image dir = sobelimg[1];
  
  feature_normalize(mag);
  feature_normalize(dir);

  image result = make_image(im.w, im.h, 3);
  for (int i = 0; i < im.w; i++) {
    for (int j = 0; j < im.h; j++) {
      set_pixel(result, i, j, 0, get_pixel(dir, i, j, 0));
      set_pixel(result, i, j, 1, get_pixel(mag, i, j, 0));
      set_pixel(result, i, j, 2, get_pixel(mag, i, j, 0));
    }
  }
  hsv_to_rgb(result);
  // free(sobel_image);
  // free_image(mag);
  // free_image(dir);

  return result;
}

// EXTRA CREDIT: Median filter

image apply_median_filter(image im, int kernel_size)
{
  image result = make_image(im.w, im.h, im.c);
  int edge = kernel_size / 2;

  for (int i = 0; i < im.w; i++) {
    for (int j = 0; j < im.h; j++) {
      for (int c = 0; c < im.c; c++) {
        int index = 0;
        float window[kernel_size * kernel_size];

        // Get the pixel value of those filters
        for (int m = -edge; m <= edge; m++) {
          for (int n = -edge; n <= edge; n++) {
            int x = i + m;
            int y = j + n;

            // Check the image boundary
            if (x >= 0 && y >= 0 && x < im.w && y < im.h) {
              window[index] = get_pixel(im, x, y, c);
              index += 1;
            }
          }
        }

        qsort(window, index, sizeof(float), cmpfunc);
        set_pixel(result, i, j, c, window[index / 2]);
      }
    }
  }
  return result;
}


// SUPER EXTRA CREDIT: Bilateral filter


image apply_bilateral_filter(image im, float sigma1, float sigma2)
{
  // image result = make_image(im.w, im.h, im.c);
  // image gauss = make_gaussian_filter(sigma1);

  // int k = (int)ceilf(6 * sigma1);
  // if (k % 2 == 0) {
  //   k += 1;
  // }

  // for (int i = 0; i < im.w; i++) {
  //   for (int j = 0; j < im.h; j++) {
  //     for (int c = 0; c < im.c; c++) {
  //       float val = 0;

  //       for (int m = 0; m < k; m++) {
  //         for (int n = 0; n < k; n++) {
  //           float gauss_val = get_pixel(gauss, m, n, 0);
  //           float f_x_y = get_pixel(im, i, j, c);
  //           float f_x_y_1 = get_pixel(im, m + i - k/2, n + j - k/2, c);
  //           float guass_val_1 = (M_2_PI * powf(sigma2, 2)) * exp(-(powf(f_x_y - f_x_y_1, 2) / 2 * powf(sigma2, 2)));
  //           val += gauss_val * guass_val_1;
  //         }
  //       }
  // 
  //     }
  //   }
  // }

  // return im;
  return make_image(1, 1, 1);
}

int cmpfunc(const void *a, const void *b) {
    float fa = *(const float*)a;
    float fb = *(const float*)b;
    
    if (fa < fb) {
      return -1;
    }
    if(fa > fb) {
      return 1;
    }
    return 0;
}