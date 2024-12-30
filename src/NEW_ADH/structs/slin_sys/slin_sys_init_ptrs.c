#include "adh.h"
void slin_sys_init_ptrs(SLIN_SYS *lin_sys, int *my_ndof_ptr, int *ndof_ptr, int *macro_ndof_ptr,
    int *my_ndof_ptr_old, int *ndof_ptr_old, int *macro_ndof_ptr_old,
    int start, int end, int nghost){

    //sets up lin sys that has a non-trivial smat physics set up
    int i = *ndof_ptr;
    int j = *my_ndof_ptr;
    // allocate and initialize some things
    //also set the lin sys variables
    lin_sys->local_range[1] = end;
    lin_sys->local_range[0] = start;
    lin_sys->local_range_old[1] = 0;
    lin_sys->local_range_old[0] = 0;
    lin_sys->nghost = nghost;

    //should be pointers back to dm->ndofs[]
    lin_sys->size = ndof_ptr;
    lin_sys->local_size = my_ndof_ptr;
    lin_sys->global_size = macro_ndof_ptr;
    lin_sys->size_old = ndof_ptr_old;
    lin_sys->local_size_old = my_ndof_ptr_old;
    lin_sys->global_size_old = macro_ndof_ptr_old;

    //when in init need to change this to check if we have more than one processor or not
    lin_sys->indptr_off_diag=NULL;
    lin_sys->cols_off_diag=NULL;
    lin_sys->vals_off_diag=NULL;
    lin_sys->ghosts=NULL;
    //should this go somewhere else? Like still in super model since we have pointer now?
    //but this should only happen over each unique model
    lin_sys->residual = (double*) tl_alloc(sizeof(double), i);
    lin_sys->indptr_diag = (int*) tl_alloc(sizeof(int), j+1);
    //actual solution of linear system is an increment within Newton iteration
    lin_sys->dsol = (double*) tl_alloc(sizeof(double), i);
    lin_sys->scale_vect = (double*) tl_alloc(sizeof(double), i);
    
    if(nghost > 0){
        lin_sys->ghosts= (int*) tl_alloc(sizeof(int), nghost);
        lin_sys->indptr_off_diag = (int*) tl_alloc(sizeof(int), j+1);
        sarray_init_int(lin_sys->indptr_off_diag,j+1);
        sarray_init_int(lin_sys->ghosts,nghost);
    }

    //init to 0's here
    sarray_init_dbl(lin_sys->residual,i);
    sarray_init_dbl(lin_sys->dsol,i);
    sarray_init_dbl(lin_sys->scale_vect,i);
}
