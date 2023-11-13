/* this routine gets a new 2d element */
#include "global_header.h"

int elem2d_new(SGRID *grid, int nalloc_inc, int nnodes_on_elem) {
    int i;            /*loop counter*/
    int the_elem;		/* the element */
    int ie;			/* loop counter */
    
    /* checks to see if it is time to allocate more elements */
    if(grid->nelems2d == grid->max_nelems2d) {
        
        /* gets more elements */
        grid->max_nelems2d += nalloc_inc;
        grid->elem2d = (SELEM_2D *) tl_realloc(sizeof(SELEM_2D), grid->max_nelems2d, grid->nelems2d, grid->elem2d);
        for(ie = grid->nelems2d; ie < grid->max_nelems2d; ie++) {

            /* Gajanan gkc adding
             * The old version did not need to allocate because it used NDONTET in the struct in selem3d.h 
             * The new version uses pointers in selem3d.h so we have to allocate them before calling init.
             */
            selem2d_alloc(&(grid->elem2d[ie]), nnodes_on_elem);
            /* --------- */

            selem2d_init(&(grid->elem2d[ie]));
        }

        if(grid->wd_flag != NULL) grid->wd_flag = (int *)tl_realloc(sizeof(int), grid->max_nelems2d, grid->nelems2d, grid->wd_flag);
        if(grid->ndim == 2)grid->elem_error = (double *)tl_realloc(sizeof(double), grid->max_nelems2d, grid->nelems2d, grid->elem_error);
    }
    
    /* gets the element */
    the_elem = grid->nelems2d;
    grid->nelems2d++;
    
    /* initialize the element data */
    selem2d_init(&(grid->elem2d[the_elem]));
    
    /* returns the element */
    return (the_elem);
}
