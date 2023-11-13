/* this routine gets a new 1d element */
#include "global_header.h"

int elem1d_new(SGRID *grid, int nalloc_inc) {
    
    int the_elem;		/* the element */
    int ie;			/* loop counter */
        
    /* checks to see if it is time to allocate more elements */
    if(grid->nelems1d == grid->max_nelems1d) {
        
        /* gets more elements */
        grid->max_nelems1d += nalloc_inc;
        selem1d_realloc_array(&(grid->elem1d), grid->nelems1d, grid->max_nelems1d, NDONSEG);
        for(ie = grid->nelems1d; ie < grid->max_nelems1d; ie++) {
            selem1d_init(&(grid->elem1d[ie]));
        }
    }
    
    /* gets the element */
    the_elem = grid->nelems1d;
    grid->nelems1d++;
    
    /* initialize the element data */
    selem1d_init(&(grid->elem1d[the_elem]));
     
    /* returns the element */
    return (the_elem);
}
