#ifndef H_SVECT2D_
#define H_SVECT2D_

typedef struct{
    double x, y;          /* the coordinates of the vector */
} SVECT2D;    

/*********************************************************/
/* struct methods -------------------------------------- */

void svect2d_init(SVECT2D *);
void svect2d_init_array(SVECT2D *, int);
void svect2d_print(FILE *, SVECT2D);
void svect2d_print_array(FILE *, SVECT2D *, int);
void svect2d_copy_array(SVECT2D *, SVECT2D *, int);
void svect2d_printScreen(SVECT2D);
void svect2d_printScreen_array(char *, SVECT2D *, int, int, char *);
double svect2d_mag(SVECT2D);
double svect2d_mag_safe(SVECT2D v);
SVECT2D svect2d_subtract(SVECT2D, SVECT2D);
SVECT2D svect2d_add(SVECT2D, SVECT2D);
void svect2d_add_array(SVECT2D *v, SVECT2D *v1, SVECT2D *v2, int size);
SVECT2D svect2d_scale(SVECT2D, double);
void svect2d_scale_array(SVECT2D *v, int size, double scale);
void svect2d_scale_replace_array(SVECT2D *v, double scale, int size);
void svect2d_nscale_array(SVECT2D *v, int size, double *scale);
double svect2d_dotp(SVECT2D, SVECT2D);
void svect2d_integrity_check(SVECT2D, int, char *);
void svect2d_integrity_check_array(SVECT2D *, int, int, char *);
SVECT2D svect2d_avg(SVECT2D, SVECT2D);
void svect2d_init_array_value_range(SVECT2D *, double, double, int, int);
SVECT2D svect2d_average_array(SVECT2D *vect, int size);
void dumpVector2D(SVECT2D *vel, int nnodes, double *u, double *v);
void printScreen_debug_svec2d(char *descript, SVECT2D *v, int n, int *global_nd_ids);

#endif
