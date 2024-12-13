// An AdH System of linear equations
#ifndef H_SLIN_SYS_
#define H_SLIN_SYS_

typedef struct {
    //Mark making changes, the structures should all be the same now
#ifdef _PETSC
    Mat A;
    KSP ksp;
    Vec B; //Petsc residual
    Vec X; //Persc solution
#endif
    //Split CSR format
    //diagonal is locally owned dofs to process
    int *indptr_diag; //pointers to number of nnz each row
    int *cols_diag; //column addresses of all nonzero entries (local to process)
    double *vals_diag; //actual matrix values

    //this actually will be empty if serial run
    int *indptr_off_diag; //pointers to number of nnz each row
    int *cols_off_diag; //column addresses of all nonzero entries (global)
    double *vals_off_diag; //actual matrix values

    //Also needs nnz for allocation purposes
    int nnz_diag, nnz_off_diag, nnz_diag_old, nnz_off_diag_old;

    //vectors
    double *residual;
    double *sol;
    double *scale_vect;

    //mark, proposes using unified solution variable
    double *sol_old;
    double *sol_older;
    //actual solution of linear system is an increment within Newton iteration
    double *dsol;

    //things that are important to transforming local to global
    int *ghosts;
    int nghost;
    int local_size; // same as my_ndofs
    int global_size; //same as macro_ndofs
    int size; //same as ndofs (SHOULD THESE THEN ALL BE POINTERS? NEED TO DECIDE SOON)
    int local_range[2];
    int local_range_old[2];

} SLIN_SYS;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void slin_sys_alloc_init_array(SLIN_SYS **lin_sys, int nlin_sys);
void lin_sys_free(SLIN_SYS *lin_sys, int nlin_sys);
void lin_sys_CSR_printScreen(SLIN_SYS *lin_sys);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif
