#ifndef H_SARRAY_
#define H_SARRAY_

/*********************************************************/
/* struct methods -------------------------------------- */

void sarray_init_int(int *, int);
void sarray_init_dbl(double *, int);

void sarray_init_value_int(int *, int, int);
void sarray_init_value_dbl(double *, int, double);

void sarray_init_value_range_int(int *, int, int, int);
void sarray_init_value_range_dbl(double *, double, int, int);

void sarray_init2one_int(int *, int);
void sarray_init2one_dbl(double *, int);

void sarray_copy_int(int *to, int *from, int nsize);
void sarray_copy_dbl(double *to, double *from, int nsize);

int sarray_sum_int(int *array, int nsize);
double sarray_sum_dbl(double *array, int nsize);

void sarray_add_replace_int(int *array1, int *array2, int nsize);
void sarray_add_replace_dbl(double *array1, double *array2, int nsize);

void sarray_add_int(int *result, int *array1, int *array2, int nsize);
void sarray_add_dbl(double *result, double *array1, double *array2, int nsize);

void sarray_subtract_replace_int(int *array1, int *array2, int nsize);
void sarray_subtract_replace_dbl(double *array1, double *array2, int nsize);

void sarray_subtract_int(int *result, int *array1, int *array2, int nsize);
void sarray_subtract_dbl(double *result, double *array1, double *array2, int nsize);

void sarray_add_scalar_int(int *array, int scalar, int nsize);
void sarray_add_scalar_dbl(double *array, double scalar, int nsize);

void sarray_subtract_scalar_int(int *array, int scalar, int nsize);
void sarray_subtract_scalar_dbl(double *array, double scalar, int nsize);

double sarray_avg_int(int *array, int nsize);
double sarray_avg_dbl(double *array, int nsize);

void sarray_scale_int(int *scaled_array, int *array, int scale, int nsize);
void sarray_scale_dbl(double *scaled_array, double *array, double scale, int nsize);

void sarray_scale_replace_int(int *array, int scale, int nsize);
void sarray_scale_replace_dbl(double *array, double scale, int nsize);

void sarray_integrity_check_int(int *array, int nsize, int linenumber, char *filename);
void sarray_integrity_check_dbl(double *array, int nsize, int linenumber, char *filename);

void sarray_printScreen_int(int *array, int nsize, char *name);
void sarray_printScreen_dbl(double *array, int nsize, char *name);

void sarray_print_int(SGRID *, FILE *, FILE *, int *, int **, int, int *, int);
void sarray_print_dbl(SGRID *, FILE *, FILE *, double *, int **, int, int *, int);

void sarray_2column_print_int(FILE *fout, int *array1, int *array2, int size);
void sarray_2column_print_dbl(FILE *fout, double *array1, double *array2, int size);

void sarray_3column_print_int(FILE *fout, int *array1, int *array2, int *array3, int size);
void sarray_3column_print_dbl(FILE *fout, double *array1, double *array2, double *array3, int size);

void sarray_vector_matrix_multiply(double *v, double **matrix, double *vout, int ndim);
void sarray_matrix_matrix_multiply(double **matrix1, double **matrix2, double **matrix_out, int ndim1, int ndim2);

int sarray_max_int(int a[], int num_elements);
double sarray_max_dbl(double a[], int num_elements);

void sarray_init_int_2d(int **to, int nrows, int ncols);
void sarray_init_double_2d(double **to, int nrows, int ncols);

void sarray_y_plus_ax_dbl(double *y, double alpha, double *x,  int n);
double sarray_dot_dbl(double *x, double *y,int n);
double sarray_l_infty_norm(double *v1, int n);
int sarray_unique_int(int *arr, int size);
int sarray_is_in_int(int *arr, int size, int val);
int sarray_argmin_int(int *arr, int n);
int sarray_argsort_int(int *my_int_arr, int *my_index_arr, int size);
int sarray_shuffle_int(int *my_int_arr, int *my_index_arr, int size);
int sarray_reverse_argsort_int(int *my_int_arr, int *my_index_arr, int size);

#endif
