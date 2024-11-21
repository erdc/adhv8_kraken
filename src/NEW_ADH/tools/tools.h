#ifndef H_STOOLS_
#define H_STOOLS_


int doesFileExist(const char *fname);

// time-stepping tools
int tc_end(SMODEL_SUPER *);
void tc_init(SMODEL_SUPER *);
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

#endif
