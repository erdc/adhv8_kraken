/* this routine gets a new 3d element */
#include "global_header.h"

int elem3d_new(SGRID *grid, int nalloc_inc, int nnodes_on_elem) {
    
    int the_elem;   /* the element */
    int ie;			/* loop counter */
    int i;          /* loop counter */
    
    /* checks to see if it is time to allocate more elements */
    if(grid->nelems3d == grid->max_nelems3d) {
        
        /* gets more elements */
        grid->max_nelems3d += nalloc_inc;
        
        grid->elem3d = (SELEM_3D *) tl_realloc(sizeof(SELEM_3D), grid->max_nelems3d, grid->nelems3d, grid->elem3d);
        grid->elem_error = (double *) tl_realloc(sizeof(double), grid->max_nelems3d, grid->nelems3d, grid->elem_error);
        grid->hyd_eddy = (double *) tl_realloc(sizeof(double), grid->max_nelems3d, grid->nelems3d, grid->hyd_eddy);
        grid->trn_diff = (double *) tl_realloc(sizeof(double), grid->max_nelems3d, grid->nelems3d, grid->trn_diff); 
        /* commented out seh, seems that we only need to init for the individual element */
        for(ie = grid->nelems3d; ie < grid->max_nelems3d; ie++) {

            /* Gajanan gkc adding
             * The old version did not need to allocate because it used NDONTET in the struct in selem3d.h 
             * The new version uses pointers in selem3d.h so we have to allocate them before calling init.
             */
            selem3d_alloc(&(grid->elem3d[ie]), nnodes_on_elem);
            /* --------- */

            selem3d_init(&(grid->elem3d[ie]));
            grid->elem_error[ie]=0.;
            grid->hyd_eddy[ie]=0.;
            grid->trn_diff[ie]=0.;
        }
    }
    
    
    /* gets the element */
    the_elem = grid->nelems3d;
    grid->nelems3d++;
    
    /* initialize the element data */
    selem3d_init(&(grid->elem3d[the_elem]));
    
    /* returns the element */
    return (the_elem);
}
