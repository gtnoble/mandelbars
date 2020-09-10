# Mandelbars

## Introduction
Mandelbars is a fractal rendering program that focuses on features that improve image quality. Examples include: high dynamic range images and anti-aliasing.

The currently implemented anti-aliasing algorithm is an adaptive Monte Carlo algorithm, with multiple variance reduction techniques available to increase speed.

The program is targeted at Linux / Unix / POSIX systems and currently has a command line interface.

## Dependencies
The only dependencies are the C and Unix standard libraries as well as the GNU Scientific Library. GSL is used for random number generation and running statistics functionality used for the anti-aliasing algorithm. The C compiler must support VLAs and complex numbers.

## Build instructions
Type `make` and copy the resulting mandelbars executable to the desired location. Currently only supports GCC, but I imagine the code will compile on any C99 compliant compiler.

### Usage
`mandelbars [-ach] [-f parameter_filename | -L real_location,imag_location,
scene_zoom] [-d exterior_stop_distance] [-e stop_std_err_mean] [-E std_err_mean_output_filenam
e] [-k kernel_scale] [-l periodicity_check_length] [-n maximum_iterations] [-p num_initial_poi
nts] [-P npoints_output_filename] [-s scattering_distribution] [-v visualization_algorithm] [-x image_x_dimension] [-y image_y_dimension] image_output_filename` 

## Command line options
`-a` - Use antithetic sampling. Improves speed in images with large numbers of monotonic sections.

`-c` - Use control variate technique. Uses information from previously calculated pixels to improve speed. Recommended when the fractal visualization is approximately continuous.

`-h` - Print the usage statement.

`-f` - Specifies Fractint parameter filename. Fractint can be used to explore the fractal to find interesting regions. The Fractint location can be exported as a parameter file and this parameter file can be loaded into mandelbars to generate a high-quality image.

`-L` - Specifies location and zoom for the fractal scene.

`-d` - Specifies the stop distance for exterior distance visualization method. You likely won't need to change this. Default is `2^64`.

`-e` - Specifies the maximum standard error for the adaptive anti-aliasing algorithm. Lower numbers produce a less noisy image. Default is `0.05`.

`-E` - Output a csv file with the standard error values for the image pixels. Used for diagnostic purposes to evaluate the anti-aliasing algorithm.

`-k` - Specifies the scale of the anti-aliasing kernel. Value is in pixels. Larger values reduce aliasing more at the expense of image sharpness. Default is `0.3`.

`-l` - Specifies the number of iterations to look back to check for periodicities in escape time visualization algorithms. Default is `20`.

`-n` - Specifies the number of iterations to perform for a given visualization algorithm. Larger numbers give better image quality at the expense of computation time. Default is `1024`.

`-p` - Specifies the number of initial points calculated by the anti-aliasing algorithm to determine the initial variance of the pixel. Default is `2`.

`-P` - Output a csv file with the number of points calculated for the individual image pixels. Used for diagnostic purposes to evaluate the anti-aliasing algorithm.

`-s` - Specifies the kernel type for the anti-aliasing algorithm. Options are `uniform`, `gaussian`, and `null`. Uniform is a rectangular "boxcar" kernel. Scale specifies the width of the kernel. Gaussian specifies a gaussian kernel. Scale specifies the standard deviation. "null" turns off anti-aliasing. Default is `gaussian`.

`-v` - Specifies the visualization algorithm. Options are `escape_time` and `exterior_distance`. Default is `escape_time`.

`-x` - Specifies the x dimension for the image, in pixels. Default is `128`.

`-y` - Specifies the y dimension for the image, in pixels. Default is `128`.

`image_output_filename` - Specifies the output filename. The output filetype is 16-bit PGM.


