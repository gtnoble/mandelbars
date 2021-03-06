struct scene_params {
  //output image Y dimension (pixels)
  int y_dim;
  //output image X dimension (pixels)
  int x_dim;
  
  //length of the imaginary interval covered by the y axis of the screen.
  //the imaginary and complex distance between pixels is calculated from this.
  double zoom;
  //where the center of the screen points to on the complex plane
  double complex center;

  double dreal_dx; //increment of real component of complex plane per pixel in screen x
  double dimag_dy; //same as above, but for the imaginary component in screen y 

  // Add these values to screen coordinates to convert to cartesian
  double screen_to_cartesian_x;
  double screen_to_cartesian_y;
};

struct render_params {
  //maximum escape time algorithm iterations
  int iter_max;
  // How many fractal iterations should we "look back on" to check for periodicity
  int periodicity_check_length;
  // Bailout distance for exterior distance visualization technique
  double exterior_stop_distance;
  double (*visualizer)(double complex, struct render_params);

  //how many how many samples are taken for the initial variance estimates
  int num_initial_pts;
  //what the standard error of the mean must be to stop the adaptive algorithm for a pixel
  double stop_std_err_mean;
  // What is the width of the AA kernel? Measure of scale varies by kernel type
  double kernel_scale;
  bool use_control_variate;
  bool use_antithetic;
  // What AA kernel should we use?
  struct scattered_point (*scatterer)(gsl_rng *);
};

struct plane_params {
  double dzdx;
  double dzdy;
  double z0;
  double max_val;
};

struct scattered_point {
  double x;
  double y;
};

struct cli_options {
  struct render_params render;
  int screen_x_dim;
  int screen_y_dim;
};

struct scene_params read_fractint_param_file(const char *filename, int xdim, int ydim); 

void parse_cli(int argc, char *argv[], struct render_params *render,
               struct scene_params *scene, char *image_output_filename[],
               char *npoints_output_filename[], char *std_err_mean_output_filename[]);

void convert_image_to_unit(struct scene_params scene, 
                          const double fimage[scene.x_dim][scene.y_dim], 
                          unsigned short uiimage[scene.x_dim][scene.y_dim]);

void write_pgm(const char *filename, 
              struct scene_params scene, 
              unsigned short image[scene.x_dim][scene.y_dim]); 

void write_csv(const char *filename, struct scene_params scene,
               double data[scene.x_dim][scene.y_dim]);

struct scene_params generate_scene(int x_dim, 
                                   int y_dim, 
                                   double zoom, 
                                   double complex center);

void render_image(struct scene_params scene, 
                  struct render_params render, 
                  double image[scene.x_dim][scene.y_dim],
                  double npoints[scene.x_dim][scene.y_dim],
                  double std_err_mean[scene.x_dim][scene.y_dim]); 

                  
double complex screen_to_scene_coords(int x, int y, struct scene_params scene);


struct scattered_point scatter_naive(gsl_rng * rng);
struct scattered_point scatter_null(gsl_rng * rng);
struct scattered_point scatter_gaussian(gsl_rng * rng); 

double visualize_escape_time(double complex complex_cooridnates, 
                             struct render_params render);

double visualize_exterior_distance(double complex complex_coordinates,
                                   struct render_params render);

double sample_naive(struct plane_params control_plane,
                    double complex complex_coordinates, 
                    struct render_params render,
                    struct scene_params scene,
                    gsl_rng * rng);

                    
struct plane_params estimate_plane_params(int nsample__x,
                                          int nsample__y,
                                          struct scene_params scene,
                                          double image[scene.x_dim][scene.y_dim]);

double evaluate_plane(int x, int y, struct plane_params plane);
