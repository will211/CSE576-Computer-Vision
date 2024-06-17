from uwimg import *

im = load_image("data/dogsmall.jpg")
a = nn_resize(im, im.w*4, im.h*4)
save_image(a, "dog4x-nn")

im = load_image("data/dogsmall.jpg")
a = bilinear_resize(im, im.w*4, im.h*4)
save_image(a, "dog4x-bl")

im = load_image("data/dog.jpg")
a = nn_resize(im, im.w//7, im.h//7)
save_image(a, "dog7th-bl")

im = load_image("data/dog.jpg")
f = make_box_filter(7)
f1 = make_highpass_filter()
f2 = make_sharpen_filter()
f3 = make_emboss_filter()
blur = convolve_image(im, f, 1)
save_image(blur, "dog-box7")

im = load_image("data/dog.jpg")
f = make_box_filter(7)
blur = convolve_image(im, f, 1)
thumb = nn_resize(blur, blur.w//7, blur.h//7)
save_image(thumb, "dogthumb")

im = load_image("data/dog.jpg")
f = make_gaussian_filter(2)
blur = convolve_image(im, f, 1)
save_image(blur, "dog-gauss2")

im = load_image("data/dog.jpg")
f = make_gaussian_filter(2)
lfreq = convolve_image(im, f, 1)
hfreq = im - lfreq
reconstruct = lfreq + hfreq
save_image(lfreq, "low-frequency")
save_image(hfreq, "high-frequency")
save_image(reconstruct, "reconstruct")

im = load_image("data/dog.jpg")
res = sobel_image(im)
mag = res[0]
feature_normalize(mag)
save_image(mag, "magnitude")

im = load_image("data/dog.jpg")
res = colorize_sobel(im)
mag = res
feature_normalize(mag)
save_image(mag, "magnitude_color")

im = load_image("data/landscape.jpg")
res = apply_median_filter(im, 3)
save_image(res, "median.jpg")

im = load_image("data/bilateral_raw.png")
res = apply_bilateral_filter(im, 2, 0.2)
save_image(res, "bilateral.jpg")

im = load_image("data/dog.jpg")
color = colorize_sobel(im)
save_image(color, "dog_color")