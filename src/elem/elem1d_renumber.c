/* this routine renumbers the 1d elements */
#include "global_header.h"

static int DEBUG = OFF;

void elem1d_renumber(SGRID *grid) {
    
    int ie;			/* loop counter over the elements */
    int next;		/* the next available node number */
    
    if(DEBUG) {
        printf("Renumbering the 1D elements.\n");
    }
    
    /* find connection to 2d elements and outward normal */
    elem1d_find_elem2d_init(grid);
    
    for(ie = 0, next = 0; ie < grid->max_nelems1d; ie++) {
        if(grid->elem1d[ie].string != UNSET_INT) {
            selem1d_copy(&(grid->elem1d[next++]), grid->elem1d[ie]);
        }
    }
    grid->nelems1d = next;
    for(ie = next; ie < grid->max_nelems1d; ie++) {
        selem1d_init(&(grid->elem1d[ie]));
    }

    for(ie = 0; ie < grid->max_nelems1d; ie++) {
        grid->elem1d[ie].id = ie;
    }
}
