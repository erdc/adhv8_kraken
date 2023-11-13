/* this routine renumbers the 2d elements */
#include "global_header.h"

static int DEBUG = OFF;

void elem2d_renumber(SGRID *grid) {
    int ie;			/* loop counter over the elements */
    int next;			/* the next available node number */
    int i, j, k, target = UNSET_INT;
    
    if(DEBUG) {
        printf("Renumbering the 2D elements.\n");
    }
   
    for (ie = 0, next = 0; ie < grid->max_nelems2d; ie++) {
        if (grid->elem2d[ie].string != UNSET_INT) {
   
            target = next++;
            selem2d_copy(&(grid->elem2d[target]), grid->elem2d[ie]);
        
            // element errors (not loaded in elem2d struct)
            if ((grid->elem_error) != NULL && (grid->ndim==2)) grid->elem_error[target] = grid->elem_error[ie];
        
            // wetting and drying factors (not loaded in elem2d struct)
            if (grid->ndim == 2) grid->wd_flag[target] = grid->wd_flag[ie];
        }
    }
    grid->nelems2d = next;
    for(ie = next; ie < grid->max_nelems2d; ie++) {
        selem2d_init(&(grid->elem2d[ie]));
    }

    for(ie = 0; ie < grid->max_nelems2d; ie++) {
        grid->elem2d[ie].id = ie;
    }

    return;
}
