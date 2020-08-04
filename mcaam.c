/*
Mandelbars - A command line mandelbrot set renderer for Unix-like systems
By: Garret Noble
Updated: June 25 2020

I am writing this program primarily to explore various methods of anti-aliasing graphics. Fractal images suffer greatly from aliasing, yet, they are simple to generate. This makes fractal images ideal for my purposes. I am also using this project to learn C and good programming practices. Hence, the extensive commenting.
*/


#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_randist.h>
#include <stdbool.h>
#include <regex.h>
#include <string.h>
#include "mcaam.h"

//maximum escape time algorithm iterations
//#define N_MAX 100
#define N_MAX 1024
//output image X dimension (pixels)
#define X_DIM 640
//output image Y dimension (pixels)
#define Y_DIM 480
//length of the imaginary interval covered by the y axis of the screen.
//the imaginary and complex distance between pixels is calculated from this.
//#define ZOOM 2.0
#define ZOOM (1 / 105.0408)
//where the center of the screen points to on the complex plane
//#define CENTER CMPLX(-0.5,0.0)
#define CENTER CMPLX(-1.25716839231694433, 0.38008325342230764)
//for the adaptive anti-aliasing algortithm
//how many how many samples are taken for the initial variance estimates
#define NUM_INITIAL_PTS 3
//what the standard error of the mean must be to stop the adaptive algorithm for a pixel
#define STOP_SD_MEAN 0.05
// Scale of antialiasing kernel. Length for boxcar, SD for gaussian
#define KERNEL_SCALE 0.3 
// Are we using planar control variate technique?
#define USE_CONTROL_VARIATE false
// Are we using antithetic variate technique?
#define USE_ANTITHETIC true
// What is the max number of iterations that we look back to check for periodicity
#define PERIODICITY_LENGTH 20
// For exterior distance method, how far away can the iterated point be before stopping
#define STOP_DISTANCE pow(2,64) 
#define SCATTERER scatter_gaussian
// How do we visualize the fractal?
#define VISUALIZER visualize_escape_time
// Output image filename
#define OUTPUT_FILENAME "image.pgm"
// Defines the maximum pixel value for a 16-bit PGM file
#define MAX_PGM_PIXEL_VALUE ((1 << 16) - 1)
// Pattern used to detect the fractint parameter file row with center and zoom data
#define FRACTINT_COORD_ZOOM_FIELD "center-mag="
// Path to fractint parameter file for testing
#define FRACTINT_SCENE_FILE "./test_param_files/fractint_params.par"

int main() {

  struct scene_params scene;
  //scene = generate_scene(X_DIM, Y_DIM, ZOOM, CENTER);
  scene = read_fractint_param_file(FRACTINT_SCENE_FILE, X_DIM, Y_DIM);
  struct render_params render;
  render.iter_max = N_MAX;
  render.num_initial_pts = NUM_INITIAL_PTS;
  render.stop_std_err_mean = STOP_SD_MEAN; 
  render.kernel_scale = KERNEL_SCALE;
  render.use_control_variate = USE_CONTROL_VARIATE;
  render.use_antithetic = USE_ANTITHETIC;
  render.periodicity_check_length = PERIODICITY_LENGTH;
  render.exterior_stop_distance = STOP_DISTANCE;
  render.visualizer = VISUALIZER;
  
  double image[scene.x_dim][scene.y_dim];
  unsigned short uiimage[scene.x_dim][scene.y_dim];
  render_image(scene, render, image);
  convert_image_to_unit(scene, image, uiimage);
  write_pgm(OUTPUT_FILENAME, scene, uiimage);
  }

struct scene_params read_fractint_param_file(const char *filename, int xdim, int ydim) {

  struct scene_params scene;

  FILE * file_pointer;
  char *line = NULL;
  size_t length = 0;
  double real_coord;
  double imag_coord;
  double fractint_zoom;
  bool is_coord_zoom_field_detected = false;
  int scan_ret;

  file_pointer = fopen(filename, "r");
  if(file_pointer == NULL) {
    printf("failed to open parameter file\n");
    exit(EXIT_FAILURE);
  }

  // Find the zoom and center data in the fractint parameter file and read it
  while(getline(&line, &length, file_pointer) != -1) {
   if(strstr(line, FRACTINT_COORD_ZOOM_FIELD) != NULL) {
    scan_ret = sscanf(line, " center-mag=%lf/%lf/%lf", 
                      &real_coord, &imag_coord, &fractint_zoom);
    is_coord_zoom_field_detected = true;
   }
  }

  if(is_coord_zoom_field_detected && scan_ret != EOF) {
     double zoom = 1 / fractint_zoom;
     double complex center = CMPLX(real_coord, imag_coord);
     scene = generate_scene(xdim, ydim, zoom, center);
  }
  else {
     printf("failed to parse parameter file!\n");
     exit(EXIT_FAILURE);
  }

  return(scene);
}

int convert_image_to_unit(struct scene_params scene, 
                          const double fimage[scene.x_dim][scene.y_dim], 
                          unsigned short uiimage[scene.x_dim][scene.y_dim]) {

  // Minimum pixel value
  double min = fimage[0][0];
  // Maximum pixel value
  double max = fimage[0][0];
  
  // Find maximim and minimum pixel values for entire image
  for(int x = 0; x < scene.x_dim; x++) {
    for(int y = 0; y < scene.y_dim; y++) {
      if(fimage[x][y] > max) max = fimage[x][y];
      if(fimage[x][y] < min) min = fimage[x][y];
      //printf("%f %f\n", max, min);
    }
  }
  //printf("%f %f\n", max, min); 

  // Scale image and convert to uint short. PGM pixels only take unit values up to 2^16
  for(int x = 0; x < scene.x_dim; x++) {
    for(int y = 0; y < scene.y_dim; y++) {
      double normalized = (fimage[x][y] - min) / (max - min); 
      uiimage[x][y] = (unsigned short) (normalized * MAX_PGM_PIXEL_VALUE);
     //printf("%hu ", uiimage[x][y]);
    }
    //printf("\n");
  }
}


int write_pgm(const char *filename, 
              struct scene_params scene, 
              unsigned short image[scene.x_dim][scene.y_dim]) {
  FILE * file_pointer;
  file_pointer = fopen(filename, "w");
  
  // Specifies PGM file type
  fprintf(file_pointer, "P2\n");
  // Specifies image dimensions
  fprintf(file_pointer, "%i %i\n", scene.x_dim, scene.y_dim); 
  // Specifies max pixel value
  fprintf(file_pointer, "%i\n", MAX_PGM_PIXEL_VALUE);
  // Write pixels
  for(int y = 0; y < scene.y_dim; y++) {
    for(int x = 0; x < scene.x_dim; x++) {
      fprintf(file_pointer, "%hu ", image[x][y]);
    }
    fprintf(file_pointer, "\n");
  }

  fclose(file_pointer);
}

void render_image(struct scene_params scene, 
                  struct render_params render,
                  double image[scene.x_dim][scene.y_dim]) {

  // Using a null plane for the control variate estimation disables control variate method
  struct plane_params null_plane;
  null_plane.dzdx = 0;
  null_plane.dzdy = 0;
  null_plane.z0 = 0;

  double complex point_coords;
  gsl_rng * rng;
  rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rstat_workspace *sample_stat_accumulator = gsl_rstat_alloc();
  struct plane_params plane;
  
    for(int x = 0; x < scene.x_dim; x++) {
      for(int y = 0; y < scene.y_dim; y++) {
        point_coords = screen_to_scene_coords(x, y, scene);
        double point;
        // Control variate plane cannot be calculated on the edge 
        bool is_on_edge = (x == 0 || y == 0);
        // Control variate technique disabled on edge and when instructed by params
        bool is_control_variate_inactive = is_on_edge || (! render.use_control_variate);
        //bool is_on_edge = true;

       int i = 0;
        while((i < render.num_initial_pts) || 
               (gsl_rstat_sd_mean(sample_stat_accumulator) > render.stop_std_err_mean)) {
          
          // Render point without control variate, if disabled
          if(is_control_variate_inactive) {
            point = sample_naive(SCATTERER, 
                                 null_plane,
                                 point_coords, 
                                 render, 
                                 scene, 
                                 rng);

                                           }
          else {
            plane = estimate_plane_params(x, y, scene, false, image);
            point = sample_naive(SCATTERER, 
                                 plane,
                                 point_coords, 
                                 render, 
                                 scene, 
                                 rng);
                                 
          }

          gsl_rstat_add(point, sample_stat_accumulator);
          i++;
        }

        if(is_control_variate_inactive) image[x][y] = gsl_rstat_mean(sample_stat_accumulator);
        // Add the MC CV integrated portion to the analytic portion (planar)
        else image[x][y] = gsl_rstat_mean(sample_stat_accumulator) + plane.z0;
        //printf("%f ", image[x][y]);
        gsl_rstat_reset(sample_stat_accumulator);
      }
      //printf("\n");
    }
  gsl_rng_free(rng);
  gsl_rstat_free(sample_stat_accumulator);
}


double visualize_escape_time(double complex complex_coordinate, 
                             struct render_params params) {
  
  double complex z = complex_coordinate;
  //cardioid and bulb check
  double x = creal(z);
  double y = cimag(z);
  double q = pow(x - .25, 2) + pow(y, 2);
  bool is_in_cardioid = (q * (q + (x - .25))) <= (.25 * pow(y, 2));
  bool is_in_bulb = pow(x + 1.0, 2) + pow(y, 2) <= 1.0/16.0;
  if(is_in_cardioid || is_in_bulb) return(1.0);

  double complex historical_z = CMPLX(0, 0);
  int count = 0;
  while((count < params.iter_max) && ((pow(creal(z), 2) + pow(cimag(z), 2)) < 4)) { 
    z = cpow(z,2) + complex_coordinate;
    
    //periodicity check
    if(z == historical_z) return(1.0);
    if((count % params.periodicity_check_length) == 0) historical_z = z;
    
    count++;
  }
  double scaled_count = (double)count / params.iter_max; 
  return(scaled_count);
}

double visualize_exterior_distance(double complex complex_coordinate,
                                   struct render_params params) {

  double complex dPdc = CMPLX(0.0, 0.0);
  double complex Pc = complex_coordinate;
  double distance;
  int stop = params.iter_max;
  bool has_bailed_out = false;

  for(int n = 0; n < stop; n++) {
   // Can't let |Pc| get too big or there will be numerical issues with log calculation
   if(pow(creal(Pc), 2) + pow(cimag(Pc), 2) >= params.exterior_stop_distance) {
   has_bailed_out = true;
   break;
   }

   dPdc = CMPLX(2.0, 0) * Pc * dPdc + CMPLX(1.0, 0);
   Pc = cpow(Pc,2) + complex_coordinate;
  }

  if(! has_bailed_out) distance = 0.0;
  else distance = cabs(Pc) / cabs(dPdc) * log(cabs(Pc)); 
  return(distance);

}
//generates a uniformly distrubuted x,y coordinate over the kernel area
//centered at center of kernel, interval is -.5 -> .5 pixel
struct scattered_point scatter_naive(gsl_rng * rng) {
  struct scattered_point scatter; 
  scatter.x = (gsl_rng_uniform(rng) - .5);
  scatter.y = (gsl_rng_uniform(rng) - .5);
  scatter.weight = 1;
  return(scatter); 
}

//generates a gaussian distributed x,y coordinate over the kernel area
//centered at center of kernel, std. dev is 1 pixel
struct scattered_point scatter_gaussian(gsl_rng * rng) {
  struct scattered_point scatter;
  scatter.x = gsl_ran_gaussian(rng, 1.0);
  scatter.y = gsl_ran_gaussian(rng, 1.0);
  scatter.weight = 1;
  return(scatter);
}

//generates a non-scattered x,y coordinate
//effectively turns off anti-aliasing
struct scattered_point scatter_null(gsl_rng * rng) {
  struct scattered_point scatter;
  scatter.x = 0;
  scatter.y = 0;
  scatter.weight = 1;
  return(scatter);
}

double sample_naive(struct scattered_point (*scatterer)(gsl_rng *),
                    struct plane_params control_plane,
                    double complex complex_coordinates, 
                    struct render_params render,
                    struct scene_params scene,
                    gsl_rng * rng) {

                      
  double point;
  struct scattered_point scatter;
  scatter = scatterer(rng);

  double scaled_scatter_x = scatter.x * render.kernel_scale;
  double scaled_scatter_y = scatter.y * render.kernel_scale;

  double complex complex_scene_scatter = CMPLX(scaled_scatter_x * scene.dreal_dx, 
                                               scaled_scatter_y * scene.dimag_dy);

  point = (render.visualizer(complex_coordinates + complex_scene_scatter, render) -
                      evaluate_plane(scaled_scatter_x, scaled_scatter_y, control_plane)) / 
                      scatter.weight;

  if(render.use_antithetic) {
    double point_antithetic;
    point_antithetic = (render.visualizer(complex_coordinates - complex_scene_scatter, 
                                          render) -
                        evaluate_plane(-scaled_scatter_x,
                                       -scaled_scatter_y, 
                                       control_plane)) / 
                        scatter.weight;

    point = (point + point_antithetic) / 2;
  }

return(point);
}


double complex screen_to_scene_coords(int x, int y, struct scene_params scene) {
  double complex coords = CMPLX((x + scene.shift_x) * scene.dreal_dx, 
                                (y + scene.shift_y) * scene.dimag_dy) + scene.center;
  return(coords);
}

struct scene_params generate_scene(int x_dim, 
                                   int y_dim, 
                                   double zoom, 
                                   double complex center) {
  struct scene_params scene;
  scene.x_dim = x_dim;
  scene.y_dim = y_dim;
  scene.zoom = zoom;
  scene.center = center;
  scene.dimag_dy = (zoom / y_dim);
  scene.dreal_dx = (zoom / y_dim);
  scene.shift_x = (x_dim / 2 - x_dim); 
  scene.shift_y = (y_dim / 2 - y_dim);
  return(scene);
}

//calculate parameters for planar function z = dzdx * x + dzdy * y + z0
//plane centered at current pixel coordinate
struct plane_params estimate_plane_params(int screen_x,
                                          int screen_y,
                                          struct scene_params scene,
                                          bool is_normalized,
                                          double image[scene.x_dim][scene.y_dim]) {

  struct plane_params plane;

  double z_corner = image[screen_x - 1][screen_y - 1];
  plane.dzdx = image[screen_x - 1][screen_y] - z_corner;
  plane.dzdy = image[screen_x][screen_y - 1] - z_corner;
  plane.z0 = z_corner - plane.dzdx - plane.dzdy;

  if(is_normalized) { 
    //make sure the extreme point in the AA kernel bounds is positive
    //if not, bias the plane so that it is
    double z_far_point = plane.dzdx * .5 + plane.dzdy * .5 + plane.z0;
    if(z_far_point < 0) plane.z0 = z_corner - z_far_point;

    //normalize the area under the plane within the kernel bounds
    //normalization is required to treat the function as a distribution
    double normalization_constant = plane.z0;
    plane.dzdx = plane.dzdx / normalization_constant;
    plane.dzdy = plane.dzdy / normalization_constant;
    plane.z0 = plane.z0 / normalization_constant;
    //calculate the greatest plane point within the AA kernel bounds
    //needed for rejection sampling of planar distribution
    plane.max_val = 0;
    for(int x = -0.5; x < 0.5; x++) {
      for(int y = -0.5; y < 0.5; y++) {
        double z = plane.dzdx * x + plane.dzdy * y + plane.z0;
        if(z > plane.max_val) plane.max_val = z;
      }
    }
} 
    return(plane);
}

double evaluate_plane(int x, int y, struct plane_params plane) {
  double value = x * plane.dzdx + y * plane.dzdy + plane.z0;
  return(value);
}

