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
#include <unistd.h>
#include <inttypes.h>
#include <ctype.h>
#include "mcaam.h"

//Default Options ------------------------------------------------------------

// Default x dimension of output image, in pixels
#define IMAGE_X_DIMENSION 128

// Default y dimension of output image, in pixels
#define IMAGE_Y_DIMENSION 128

//maximum escape time algorithm iterations
#define N_MAX 1024

//for the adaptive anti-aliasing algortithm
//how many how many samples are taken for the initial variance estimates
#define NUM_INITIAL_PTS 2

//what the standard error of the mean must be to stop the adaptive algorithm for a pixel
#define STOP_SD_MEAN 0.05

// Scale of antialiasing kernel. Length for boxcar, SD for gaussian
#define KERNEL_SCALE 0.3 

// What is the max number of iterations that we look back to check for periodicity
#define PERIODICITY_LENGTH 20

// For exterior distance method, how far away can the iterated point be before stopping
#define STOP_DISTANCE pow(2,64) 

// How should the scattered points be distributed?
#define SCATTERER scatter_gaussian

// How do we visualize the fractal?
#define VISUALIZER visualize_escape_time

// Magic Numbers / Values ---------------------------------------------------

// Defines the maximum pixel value for a 16-bit PGM file
#define MAX_PGM_PIXEL_VALUE ((1 << 16) - 1)
// Pattern used to detect the fractint parameter file row with center and zoom data
#define FRACTINT_COORD_ZOOM_FIELD "center-mag="
// How data is arranged in the fractint parameter file location/zoom row
#define FRACTINT_COORD_ZOOM_FORMAT " center-mag=%lf/%lf/%lf"
// Command line option flags and flag types
#define CLI_OPTION_FORMAT "hn:p:e:k:cal:d:v:x:y:f:s:L:P:E:"

int main(int argc, char *argv[]) {

  struct scene_params scene;
  struct render_params render;
  char *image_filename = NULL;
  char *npoints_filename = NULL;
  char *std_err_mean_filename = NULL;
  parse_cli(argc, argv, &render, &scene, 
            &image_filename, &npoints_filename, &std_err_mean_filename);

  double (*image)[scene.x_dim] = malloc(sizeof(*image) * scene.y_dim);
  if(image == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for framebuffer.");
    exit(EXIT_FAILURE);
  }

  unsigned short (*uiimage)[scene.x_dim] = malloc(sizeof(*uiimage) * scene.y_dim);
  if(uiimage == NULL) {
    fprintf(stderr, "ERROR: Could not allocate memory for output image.");
    exit(EXIT_FAILURE);
  }

  // If requested, allocate memory for the npoints data matrix
  double (*npoints)[scene.x_dim] = NULL;
  if(npoints_filename != NULL) {
    if((npoints = malloc(sizeof(*npoints) * scene.y_dim)) == NULL) {
      fprintf(stderr, "ERROR: Could not allocate memory for n_points array");
      exit(EXIT_FAILURE);
    }
  }
  else npoints = NULL;

  // If requested, allocate memory for the standard error of the mean data matrix
  double (*std_err_mean)[scene.x_dim] = NULL;
  if(std_err_mean_filename != NULL) {
    if((std_err_mean = malloc(sizeof(*std_err_mean) * scene.y_dim)) == NULL) {
      fprintf(stderr, "ERROR: Could not allocate memory for std. err. mean array");
      exit(EXIT_FAILURE);
    }
  }
  else std_err_mean = NULL;

  render_image(scene, render, image, npoints, std_err_mean);
  convert_image_to_unit(scene, image, uiimage);
  write_pgm(image_filename, scene, uiimage);

  // Write statistical data files, if the memory has been allocated for them.
  if(npoints != NULL) write_csv(npoints_filename, scene, npoints);
  if(std_err_mean != NULL) write_csv(std_err_mean_filename, scene, std_err_mean);

  free(image);
  free(uiimage);
  free(npoints);
  free(std_err_mean);
}

void parse_cli(int argc,
               char *argv[], 
               struct render_params *render, 
               struct scene_params *scene,
               char *image_output_filename[],
               char *npoints_output_filename[],
               char *std_err_mean_output_filename[]) {

                             
  // Default render options
  render->iter_max = N_MAX;
  render->num_initial_pts = NUM_INITIAL_PTS;
  render->stop_std_err_mean = STOP_SD_MEAN; 
  render->kernel_scale = KERNEL_SCALE;
  render->use_control_variate = false;
  render->use_antithetic = false;
  render->periodicity_check_length = PERIODICITY_LENGTH;
  render->exterior_stop_distance = STOP_DISTANCE;
  render->visualizer = VISUALIZER;
  render->scatterer = SCATTERER;

  char *parameter_filename;

  // Scene setup varaibles
  // Used when the scene is specified via the CLI
  int x_screen_dimension = IMAGE_Y_DIMENSION;
  int y_screen_dimension = IMAGE_X_DIMENSION;
  double real_location;
  double imag_location;
  double scene_zoom;
  
  // We need to turn off anti-aliasing for the null scatterer, so we keep track of it
  bool is_scatterer_null = false;

  // The scene parameters can be specified directly in the CLI or in a parameter file
  // We need to keep track of this because it doesn't make sense to have both options
  // selected
  bool does_cli_specify_parameter_file = false;
  bool does_cli_specify_scene = false;

  *image_output_filename = NULL;
  *npoints_output_filename = NULL;
  *std_err_mean_output_filename = NULL;
  parameter_filename = NULL;
  
  int command_flag;

  // Helper function to parse command line options that yield a double value
  void double_parse(double *render_parameter) {
    if(! sscanf(optarg, "%lf", render_parameter)) {
      fprintf(stderr, "Option %c must have a real number argument\n", command_flag);
      exit(EXIT_FAILURE);
    }
  }
  
  //Helper function to parse command line options that yield a positive int value
  void int_parse(int *render_parameter) {
    if(! sscanf(optarg, "%d", render_parameter) || *render_parameter <= 0) {
      fprintf(stderr, "Option %c must have a positive integer argument\n", command_flag);
      exit(EXIT_FAILURE);
    }
  }
  
  void print_usage() {

    char usage_parameters[] = "[-ach] [-f parameter_filename | -L real_location,imag_location,scene_zoom] [-d exterior_stop_distance] [-e stop_std_err_mean] [-E std_err_mean_output_filename] [-k kernel_scale] [-l periodicity_check_length] [-n maximum_iterations] [-p num_initial_points] [-s scattering_distribution] [-v visualization_algorithm] [-x image_x_dimension] [-y image_y_dimension] image_output_filename";

    printf("Usage: %s %s\n", argv[0], usage_parameters);
  }
  
  void suggest_help() {
    fprintf(stderr, "Try '%s -h' for usage information.\n", argv[0]);
  }

  while((command_flag = getopt(argc, argv, CLI_OPTION_FORMAT)) != -1) {
    switch(command_flag) {
      // Print usage statement
      case 'h':
        print_usage();
        exit(EXIT_SUCCESS);
        break;
      // Select max iterations
      case 'n':
        int_parse(&render->iter_max);
        break;
      // Select number of points used to calculate initial pixel error
      case 'p':
        int_parse(&render->num_initial_pts);
        if(render->num_initial_pts < 2) {
          fprintf(stderr, "ERROR: The number of initial points must be 2 or greater\n");
          suggest_help();
          exit(EXIT_FAILURE);
        }
        break;
      // Output filename for number of points calculated for image pixels
      case 'P':
        *npoints_output_filename = optarg;
        break;
      // Set the maximum standard error of the mean for the AA algorithm
      case 'e':
        double_parse(&render->stop_std_err_mean);
        break;
      // Output filename for standard error of mean for image pixels
      case 'E':
        *std_err_mean_output_filename = optarg;
        break;
      // Set the AA kernel scale
      case 'k':
        double_parse(&render->kernel_scale);
        break;
      // Activate control variate technique
      case 'c':
        render->use_control_variate = true;
        break;
      // Activate antithetic variate technique
      case 'a':
        render->use_antithetic = true;
        break;
      // Set the periodicity check length
      case 'l':
        int_parse(&render->periodicity_check_length);
        break;
      // Set the stop distance for the exterior distance visualization method
      case 'd':
        double_parse(&render->exterior_stop_distance);
        break;
      // Select the scattering distribution / kernel shape
      case 's':
        if(! strcmp(optarg, "uniform"))
          render->scatterer = scatter_naive;
        else if(! strcmp(optarg, "gaussian"))
          render->scatterer = scatter_gaussian;
        else if(! strcmp(optarg, "null")) {
          render->scatterer = scatter_null;
          is_scatterer_null = true;
        }
        else {
          fprintf(stderr, "%s is not a valid scattering distribution\n", optarg);
          suggest_help(); 
          exit(EXIT_FAILURE);
        }
        break;
      // Select the fractal visualization algorithm
      case 'v':
        if(! strcmp(optarg, "escape_time"))
          render->visualizer = visualize_escape_time;
        else if(! strcmp(optarg, "exterior_distance"))
          render->visualizer = visualize_exterior_distance;
        else {
          fprintf(stderr, "%s is not a valid visualization method\n", optarg);
          suggest_help();
          exit(EXIT_FAILURE);
        }
        break;
      // Set the output image x dimension
      case 'x':
        int_parse(&x_screen_dimension);
        break;
      // Set the output image y dimension
      case 'y':
        int_parse(&y_screen_dimension);
        break;
      // Select fractint parameter filename (if desired)
      case 'f':
        parameter_filename = optarg;
        does_cli_specify_parameter_file = true;
        break;
      // Select fractal location and scale (if not using fractint parameter file)
      // Parses location parameters from CLI using the CSV format "x,y,zoom"
      case 'L':
        if(! sscanf(optarg, "%lf,%lf,%lf", &real_location, &imag_location, &scene_zoom)) {
          fprintf(stderr, "Failed to parse location parameters\n");
          suggest_help();
          exit(EXIT_FAILURE);
        }
        else does_cli_specify_scene = true;
        break;

      case '?':
        if(optopt == command_flag) {
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
          suggest_help();
          exit(EXIT_FAILURE);
        }
        else if(isprint(optopt)) {
          fprintf(stderr, "Unknown option '-%c'.\n", optopt);
          suggest_help();
          exit(EXIT_FAILURE);
        }
        else{
          fprintf(stderr, "Unknown option character '\\x%x'.\n", optopt);
          suggest_help();
          exit(EXIT_FAILURE);
        }
        default:
          exit(EXIT_FAILURE);
          suggest_help();
    }
  }

  if(optind >= argc) {
    fprintf(stderr, "No output filename specified.\n");
    suggest_help();
    exit(EXIT_FAILURE);
  }

  *image_output_filename = argv[optind];

  // Turn off anti-aliasing, because it is meaningless if you aren't scattering samples
  if(is_scatterer_null) {
  render->use_antithetic = false;
  render->use_control_variate = false;
  }

  // Parses location parameters from either CLI or from a parameter file
  if(does_cli_specify_parameter_file && does_cli_specify_scene) {
    fprintf(stderr, "Must specify either scene parameters or a parameter file, not both\n");
    suggest_help();
    exit(EXIT_FAILURE);
  }
  if(! (does_cli_specify_parameter_file || does_cli_specify_scene)) {
    fprintf(stderr, "Must specify the scene parameters directly (-z) or specify a parameter file (-f)\n");
    suggest_help();
    exit(EXIT_FAILURE);
  }
  if(does_cli_specify_parameter_file) {
    *scene = read_fractint_param_file(parameter_filename, 
                                     x_screen_dimension, 
                                     y_screen_dimension);
  }
  if(does_cli_specify_scene) {
    *scene = generate_scene(x_screen_dimension,
                           y_screen_dimension,
                           scene_zoom,
                           CMPLX(real_location, imag_location));
  }

}

struct scene_params read_fractint_param_file(const char *filename, int xdim, int ydim) {

  struct scene_params scene;

  FILE *file_pointer;
  char *line = NULL;
  size_t length = 0;

  double real_coord;
  double imag_coord;
  double fractint_zoom;
  int scan_ret = 0;

  file_pointer = fopen(filename, "r");
  if(file_pointer == NULL) {
    fprintf(stderr, "failed to open parameter file\n");
    exit(EXIT_FAILURE);
  }

  // Attempt to find the zoom and center data in the fractint parameter file and read it
  while(getline(&line, &length, file_pointer) != -1) {
   if(strstr(line, FRACTINT_COORD_ZOOM_FIELD) != NULL) {
    scan_ret = sscanf(line, FRACTINT_COORD_ZOOM_FORMAT, 
                      &real_coord, &imag_coord, &fractint_zoom);

    if(scan_ret != 3) {
      fprintf(stderr, "Could not parse location field in parameter file.\n");
      exit(EXIT_FAILURE);
    }

    break;
   }
  }

  // If the parameter file was successfully read, convert fractint parameters to mandelbars
  // parameters and generate a scene
  if(scan_ret) {
     double zoom = 1 / fractint_zoom;
     double complex center = CMPLX(real_coord, imag_coord);
     scene = generate_scene(xdim, ydim, zoom, center);
  }
  else {
     fprintf(stderr, "Failed to detect location parameters in parameter file!\n");
     exit(EXIT_FAILURE);
  }

  return(scene);
}

void convert_image_to_unit(struct scene_params scene, 
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


void write_pgm(const char *filename, 
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

// Write a csv file containing the matrix of statistical data.
void write_csv(const char *filename, struct scene_params scene,
               double data[scene.x_dim][scene.y_dim]) {

  FILE * file_pointer;
  file_pointer = fopen(filename, "w");

  for(int y = 0; y < scene.y_dim; y++) {
    for(int x = 0; x < scene.x_dim; x++) {
      // If the data point is the last on the row, omit the trailing comma
      if(x == (scene.x_dim - 1)) fprintf(file_pointer, "%f", data[x][y]);
      else fprintf(file_pointer, "%f,", data[x][y]);
    }
    fprintf(file_pointer, "\n");
  }

  fclose(file_pointer);
}


//TODO Multithread this function, perhaps on a per-row basis
void render_image(struct scene_params scene, 
                  struct render_params render,
                  double image[scene.x_dim][scene.y_dim],
                  double npoints[scene.x_dim][scene.y_dim],
                  double std_err_mean[scene.x_dim][scene.y_dim]) {

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
            point = sample_naive(null_plane,
                                 point_coords, 
                                 render, 
                                 scene, 
                                 rng);
          }
          // Otherwise use planar control variate technique
          else {
            // Estimate the parameters of a plane that approximates the fractal visualization
            //function.
            plane = estimate_plane_params(x, y, scene, image);
            // Subtract the planar approximation from the fractal visualization and evaluate
            // a scattered point (control variate technique).
            // This is useful because subtracting the approximation can reduce the variance
            // of the resulting function. This makes the monte-carlo integration of that 
            // portion more precise for a given number of points. The MC integrated function 
            // is then added to the analytically integrated planar approximation and the sum 
            // is the full convolution integral.
            point = sample_naive(plane,
                                 point_coords, 
                                 render, 
                                 scene, 
                                 rng);
          }

          gsl_rstat_add(point, sample_stat_accumulator);
          i++;
        }

        if(is_control_variate_inactive) image[x][y] = gsl_rstat_mean(sample_stat_accumulator);
        // If control variate technique used, 
        //add the MC CV integrated portion to the analytic portion (planar)
        else image[x][y] = gsl_rstat_mean(sample_stat_accumulator) + plane.z0;
        
        if(npoints != NULL) npoints[x][y] = gsl_rstat_n(sample_stat_accumulator);
        if(std_err_mean != NULL) 
          std_err_mean[x][y] = gsl_rstat_sd_mean(sample_stat_accumulator);

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

  // This never escapes, so its a safe initialization value
  double complex historical_z = CMPLX(0, 0);
  
  // Typical mandelbrot set escape time algorithm
  int count = 0;
  while((count < params.iter_max) && ((pow(creal(z), 2) + pow(cimag(z), 2)) < 4)) { 
    z = cpow(z,2) + complex_coordinate;
    
    //periodicity check
    if(z == historical_z) return(1.0);
    if((count % params.periodicity_check_length) == 0) historical_z = z;
    
    count++;
  }

  // Scaling the iteration count relative to the max iteration count makes specifying the
  //stopping std. err. mean more intuitive. It can be speficied as a fraction of the 
  //range of values.
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
  return(scatter); 
}

//generates a gaussian distributed x,y coordinate over the kernel area
//centered at center of kernel, std. dev is 1 pixel
struct scattered_point scatter_gaussian(gsl_rng * rng) {
  struct scattered_point scatter;
  scatter.x = gsl_ran_gaussian(rng, 1.0);
  scatter.y = gsl_ran_gaussian(rng, 1.0);
  return(scatter);
}

//generates a non-scattered x,y coordinate
//effectively turns off anti-aliasing
struct scattered_point scatter_null(gsl_rng * rng) {
  struct scattered_point scatter;
  scatter.x = 0;
  scatter.y = 0;
  return(scatter);
}

double sample_naive(struct plane_params control_plane,
                    double complex complex_coordinates, 
                    struct render_params render,
                    struct scene_params scene,
                    gsl_rng * rng) {

                    
                      
  double point_value;
  struct scattered_point scatter;

  // Scatter a point according to kernel distribution and specified scaling
  scatter = render.scatterer(rng);
  double scaled_scatter_x = scatter.x * render.kernel_scale;
  double scaled_scatter_y = scatter.y * render.kernel_scale;
  double complex complex_scene_scatter = CMPLX(scaled_scatter_x * scene.dreal_dx, 
                                               scaled_scatter_y * scene.dimag_dy);

  // Evaluate point
  point_value = (render.visualizer(complex_coordinates + complex_scene_scatter, render) -
                      evaluate_plane(scaled_scatter_x, scaled_scatter_y, control_plane)); 

  // Evaluate a point antithetic to the original point and average, if specified.
  // Antithetic sampling offers improved performance if the fractal visualization function
  // is convex (or approximately so) over much of the image.
  if(render.use_antithetic) {
    double point_antithetic_value;
    point_antithetic_value = (render.visualizer(complex_coordinates - complex_scene_scatter, 
                                          render) -
                        evaluate_plane(-scaled_scatter_x,
                                       -scaled_scatter_y, 
                                       control_plane)); 

    point_value = (point_value + point_antithetic_value) / 2;
  }

return(point_value);
}

// Converts screen coordinates to scene coordinates, with the center of the screen
// corresponding to the scene center location parameter
double complex screen_to_scene_coords(int x, int y, struct scene_params scene) {
  double complex coords = CMPLX((x + scene.screen_to_cartesian_x) * scene.dreal_dx, 
                                (y + scene.screen_to_cartesian_y) * scene.dimag_dy) + scene.center;
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
  scene.screen_to_cartesian_x = (x_dim / 2 - x_dim); 
  scene.screen_to_cartesian_y = (y_dim / 2 - y_dim);
  return(scene);
}

// For a given pixel location, calculate parameters for a planar approximation function in
// screen coordinates: z = dzdx * x + dzdy * y + z0 using values from the three previously 
// calculated adjacent pixels.
// plane centered at given pixel location
struct plane_params estimate_plane_params(int screen_x,
                                          int screen_y,
                                          struct scene_params scene,
                                          double image[scene.x_dim][scene.y_dim]) {

  struct plane_params plane;

  double z_corner = image[screen_x - 1][screen_y - 1];
  plane.dzdx = image[screen_x - 1][screen_y] - z_corner;
  plane.dzdy = image[screen_x][screen_y - 1] - z_corner;
  plane.z0 = z_corner - plane.dzdx - plane.dzdy;

  return(plane);
}

double evaluate_plane(int x, int y, struct plane_params plane) {
  double value = x * plane.dzdx + y * plane.dzdy + plane.z0;
  return(value);
}

