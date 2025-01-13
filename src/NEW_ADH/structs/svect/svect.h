#ifndef H_SVECT3D_
#define H_SVECT3D_

typedef struct {
      double x, y, z;       /* the coordinates of the vector */
} SVECT;          

/*********************************************************/
/* struct methods -------------------------------------- */

void svect_init(SVECT *);
void svect_init_value(SVECT *v, double x, double y, double z);
void svect_init_array(SVECT *, int);
SVECT svect_avg(SVECT v1, SVECT v2);
void svect_print(FILE *, SVECT);
void svect_print_array(FILE *, SVECT *, int);
void svect_copy_array(SVECT *, SVECT *, int);
void svect_copy(SVECT *vdest, SVECT vorig);
void svect_printScreen(SVECT, char *);
void svect_printScreen_array(char *, SVECT *, char *, int, int, char *);
double svect_mag(SVECT);
double svect_mag_safe(SVECT);
SVECT svect_subtract(SVECT, SVECT);
SVECT svect_add(SVECT, SVECT);
SVECT svect_scale(SVECT, double);
void svect_scale_replace(SVECT *v, double scale);
void svect_scale_replace_array(SVECT *v, double scale, int size);
void svect_add_array2(SVECT *v, SVECT *v1, SVECT *v2, int size);
SVECT svect_subtract_array(SVECT *vect, int size);
void svect_subtract_array2(SVECT *v, SVECT *v1, SVECT *v2, int size);
SVECT svect_average_array(SVECT *vect, int size);
void svect_scale_replace(SVECT *, double);
double svect_dotp(SVECT, SVECT);
void svect_dotp_array(SVECT *v1, SVECT *v2, int size, double *result);
SVECT svect_sum_array(SVECT *, int);
void svect_integrity_check(SVECT, int, char *);
void svect_integrity_check_array(SVECT *, int, int, char *);
SVECT svect_cross(SVECT v1, SVECT v2);
void printScreen_debug_svect(char *descript, SVECT *v, int n, int *global_nd_ids);
void printScreen_debug_vec(char *descript, SVECT *f, int n);

#endif
