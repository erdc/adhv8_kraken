#ifndef H_STOOLS_
#define H_STOOLS_


int doesFileExist(const char *fname);

// time-stepping tools
int tc_end(SMODEL_DESIGN *);
void tc_init(SMODEL_DESIGN *);
void tc_timeunits(FILE *, int);
double tc_eval_series(SSERIES, int, double, int);
#ifdef _ADH_GROUNDWATER
double tc_eval_series_slope(SSERIES, int, double, int);
#endif
double tc_conversion_factor(int, int);
double tc_trap_area(double, double, double, double);
void tc_scale(double *
#ifdef _MESSG
    , MPI_Comm
#endif
    );


// io stuff
char * build_filename(char *, const int, const char *, const char *);
char * build_filename2(char *, const int, const char *, const char *, const int, const char *, const int);
char * build_filename3(char *, const int, const char *, const char *, const int, const char *, const int, const char *, const int);
FILE * io_fopen(const char *filename, const char *mode, const int abort);
int read_int_field(SIO io, char **data);
int read_query_int_field(SIO io, char **data, int *query);
int read_int_field_custom(SIO io, char **data, int *query, char *field_name,int field_idx, int abort);

int io_read_error(SIO io, const char *msg, const int abort);
int io_save_line(SIO *io, FILE *fp, const char *filename, const char *line);
int strip_comments(char *fileline);
int strip_white(char **string);
int get_index_as_text(const int index, char *text);
char * convert_to_uppercase(char *text);
void pdata(char *proj_name, char *run_name, FILE *fp_out, int o_flag, char *title, char *initials, char *begtype, char *meshdim, int lnodes, int lnelems);
void pdata_index(SIO *io, FILE *fp_out, int o_flag, char *title, char *initials, char *begtype, char *meshdim, int lnodes, int lnelems, int index);
void phot(FILE *fp_out, char *title, char *begtype, char *meshdim, int lnodes, int lelems);
void phot_index(FILE *fp_out, char *title, char *begtype, char *meshdim, int lnodes, int lelems, int index);
void print_header_sw2d(SMODEL_SUPER *mod, FILE * fp_out,  int ps_flag);
void print_header_sw3d(SMODEL_SUPER *mod, FILE * fp_out,  int ps_flag);
void print_header_sediment(SMODEL_SUPER *mod, FILE *fp_out,  int ps_flag, int index);
#ifdef _ADH_GROUNDWATER
void print_header_gw(SMODEL_SUPER *mod, FILE * fp_out, int ps_flag);
#endif
void print_header(SMODEL_SUPER *mod, FILE * fp_out, int ps_flag);
SGRID create_rectangular_grid(double xmin, double xmax, double ymin, double ymax, int npx, int npy,
 double theta, int dz, double a0, double ax, double ax2, double ay, double ay2, double axy,
 double ax2y, double axy2, double ax2y2, int flag3d );

//random routines for arrays and integers
int solv_isnan(double value);
int solv_isinf(double value);
double l2_error(double *v1, double *v2, int n);
double linf_error(double *v1, double *v2, int n);
bool is_near(double x, double y);
double l2_norm(double *v, int size);
int binary_search_part(int *arr, int start, int end, int target);
double max_dbl(double a, double b);
int compare_ints(const void *a, const void *b);


void printScreen_dble_array(char * descript, double *array, int size,  int linenumber, char *filename);
void printScreen_int_array(char * descript, int *array, int size, int linenumber, char *filename);

void printScreen_debug2_dbl(char *descript, double *f, int n, int *global_nd_ids);
void printScreen_debug_int(char *descript, int *f, int n);
void printScreen_debug_dbl(char *descript, double *f, int n);
void printScreen_debug2_dbl(char *descript, double *f, int n, int *global_nd_ids);
void printScreen_rhs_3dof(char *string, int nnodes, int ie, int *nodes, double *elem_rhs);
void Is_Double_Inf_or_NaN(double X,char *filename,int linenumber);
void Is_DoubleArray_Inf_or_NaN(double *X,int arraybounds,char *filename,int linenumber);

SVECT2D tl_bendway_correction(
  SVECT2D * grad_shp,
  SVECT2D * elem_vel,
  double *elem_head,
  double *elem_c,
  double prop1,
  double prop2,
  double drying_lower_limit,
  double roughness,
  double fluid_density,
  int routine_flag,     /* flag indicator indicating what routine is call  bwc  */
  int ie
);
double get_coriolis_angular_speed(double coriolis_factor);
int is_double_small( double value ) ;
double fe_get_supg_tau_sw(int nnodes, SVECT *nodes, double g, double elem_avg_depth, double elem_avg_u, double elem_avg_v,
                                 double elem_avg_w,double *grad_shp_x, double *grad_shp_y, double *grad_shp_z, double djac,
                                 double alpha, int ndim, int tau_method_flag, int le_method_flag);
double get_element_length(int nnodes, SVECT *node, double elem_avg_u, double elem_avg_v, double elem_avg_w,
                          double *grad_shp_x, double *grad_shp_y, double *grad_shp_z, double volume_or_area,
                          int ndim, int method_flag);

double tl_find_grid_mass_elem2d(double density, STR_VALUE *str, SSERIES *series_head, double *depth, SGRID *grid, SFLAGS flags);
double tl_find_grid_mass_error_elem2d(double density, double *depth, SVECT2D *vel, SGRID *grid, SFLAGS flags, double initial_grid_mass, SSERIES *series_head, STR_VALUE *str, double dt, double *total_time_mass_flux_T);
double tl_find_grid_mass_elem3d(double density, SGRID *grid, double *displacement);
double tl_find_3d_grid_mass_error(STR_VALUE *str, SSERIES *series_head, double initial_grid_mass, double density, SGRID *grid, SVECT *vel, double *displacement, double *old_displacement, double *older_displacement, double dt, double *new_grid_mass, double *total_time_mass_flux_T);
double fe_3m2d_dry_wet_factornew(SVECT2D *x, double s, double *h, double djac, int *num_dry);

#endif
