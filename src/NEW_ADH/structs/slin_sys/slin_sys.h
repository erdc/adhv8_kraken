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
    double *scale_vect;
    //actual solution of linear system is an increment within Newton iteration
    double *dsol;

    //things that are important to transforming local to global
    int *ghosts;
    int nghost;
    int *local_size; // same as my_ndofs, pointer back to design model
    int *global_size; //same as macro_ndofs
    int *size; //same as ndofs
    int local_range[2];
    int local_range_old[2];

    //old ones for refinement maybe, should also be pointer??
    int *local_size_old;
    int *size_old; 
    int *global_size_old;

} SLIN_SYS;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Methods
void slin_sys_alloc_array(SLIN_SYS **lin_sys, int nlin_sys);
void slin_sys_init_ghosts(SLIN_SYS *lin_sys, SGRID *grid, int *fmap);
void slin_sys_init_ptrs(SLIN_SYS *lin_sys, int *my_ndof_ptr, int *ndof_ptr, int *macro_ndof_ptr,
    int *my_ndof_ptr_old, int *ndof_ptr_old, int *macro_ndof_ptr_old,
    int start, int end, int nghost);
//void slin_sys_init_sparsity_mono(SLIN_SYS *lin_sys, int *elem3d_physics_mat_id, 
//    int *elem2d_physics_mat_id, int *elem1d_physics_mat_id,
//    SMAT_PHYSICS *elem3d_physics_mat, SMAT_PHYSICS *elem2d_physics_mat,
//    SMAT_PHYSICS *elem1d_physics_mat, int *node_physics_mat_id,
//    SMAT_PHYSICS *node_physics_mat, SGRID *grid, int *fmap);
void slin_sys_init_sparsity_mono(SLIN_SYS *lin_sys, int *elem3d_physics_mat_id, 
    int *elem2d_physics_mat_id, int *elem1d_physics_mat_id,
    SMAT_PHYSICS *elem3d_physics_mat, SMAT_PHYSICS *elem2d_physics_mat,
    SMAT_PHYSICS *elem1d_physics_mat,
    SMAT_PHYSICS **node_physics_mat, SGRID *grid, int *fmap);
void slin_sys_allocate_petsc_objects(SLIN_SYS *lin_sys);
void slin_sys_free_array(SLIN_SYS *lin_sys, int nlin_sys);
void slin_sys_free(SLIN_SYS *lin_sys);
void lin_sys_free(SLIN_SYS *lin_sys, int nlin_sys);
void lin_sys_CSR_printScreen(SLIN_SYS *lin_sys);
void slin_sys_init_array(SLIN_SYS *lin_sys, int nlin_sys);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif
