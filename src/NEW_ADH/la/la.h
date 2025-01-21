#ifndef H_LA_
#define H_LA_

//matrix_utilities
//void allocate_adh_system(SMODEL_SUPER *sm);
//void fe_allocate_initialize_linear_system(SMODEL_SUPER *sm);
//void create_sparsity_split_CSR(SMODEL_SUPER *sm, SGRID *grid);
//void create_sparsity_split_CSR(SLIN_SYS *lin_sys, SMODEL_SUPER *sm, SGRID *grid);
//int unique(int *arr, int size);
//int compare_ints(const void *a, const void *b);
void Screen_print_CSR(int *indptr, int *cols, double *vals, int nrow);
void apply_Dirichlet_BC(SMODEL_SUPER *sm);
//void allocate_petsc_objects(SMODEL_SUPER *sm);
//void allocate_petsc_objects(SLIN_SYS *lin_sys);
void split_CSR_mat_vec_mult(double *Ax, int *indptr_diag, int *cols_diag, double *vals_diag, 
  int *indptr_off_diag, int *cols_off_diag, double *vals_off_diag,
  double *x, int nrows, int *ghosts, int nghost);
void check_diag(SMODEL_SUPER *sm);


// bcgstab solver routines
int solve_linear_sys_bcgstab(double *x, int *indptr_diag, int *cols_diag, double *vals_diag, 
  int *indptr_off_diag, int *cols_off_diag,double *vals_off_diag, double *b,
   double *scale_vect, int local_size, int size, int rank,int *ghosts, int nghost);
void free_bcgstab(void);
int prep_umfpack(int *indptr_diag, int *cols_diag, double *vals_diag, int nrow);
int solve_umfpack(double *x, int *indptr_diag, int *cols_diag, double *vals_diag, double *b, int nrow);
void umfpack_clear(void);
//void copy_array_to_from(double *vec_to, double *vec_from, int n);
//void zero_dbl_array(double *v, int n);
//void y_plus_ax(double *y, double alpha, double *x,int n );
//double get_global_max(double x);
//double l_infty_norm(int n, double *v1);
//double dot_dbl_array(int n, double *x,double *y);
//double messg_dsum(double x);
//void scale_dbl_array(double *v,  double alpha, int n);
//double max_dbl(double a, double b);

//linear system scaler
void scale_linear_system(int *indptr_diag, int *cols_diag, double *vals_diag, 
  int *indptr_off_diag, int *cols_off_diag,double *vals_off_diag, double *b, 
  double *x, double *scale_vect, int local_size, int size, int rank,int *ghosts, int nghost);
void unscale_linear_system(double *x, double *x0, double *scale_vect, int local_size);


#endif
