#ifndef H_GRID_INTERPOLATOR_
#define H_GRID_INTERPOLATOR_

/***********************************************************/
/***********************************************************/
/***********************************************************/

void read_adh_grid(char *adh_root, SGRID **adh_grid);

FILE* utility_open_input_file(const char* filename, int is_bin);

FILE* utility_open_output_file(const char* filename, int is_bin);

void interpolate_grid(double *out_grid, SGRID *adh_grid, float *vals,
        int nnodes, double x1, double x2, double y1, double y2, int ni, int nj);

void utility_write_winds(FILE *fp, double *x_wind, double *y_wind, int nnodes, int time, int is_bin);

///***********************************************************/
///***********************************************************/
///***********************************************************/

#endif
