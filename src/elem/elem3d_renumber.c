/* this routine renumbers the 3d elements */
#include "global_header.h"

static int DEBUG = OFF;

void elem3d_renumber(SGRID *grid) {
    int ie;			/* loop counter over the elements */
    int next;			/* the next available node number */
    int target =  UNSET_INT;
    
    if(DEBUG) {
        printf("Renumbering the 3D elements.\n");
    }
    
    for(ie = 0, next = 0; ie < grid->max_nelems3d; ie++) {
        if(grid->elem3d[ie].mat != UNSET_INT) {
            target = next++;
            selem3d_copy(&(grid->elem3d[target]), grid->elem3d[ie]);
            if (grid->elem_error != NULL) grid->elem_error[target] = grid->elem_error[ie];
        }
        
    }
    grid->nelems3d = next;
    for(ie = next; ie < grid->max_nelems3d; ie++) {
        selem3d_init(&(grid->elem3d[ie]));
    }
    for(ie = 0; ie < grid->max_nelems3d; ie++) {
        grid->elem3d[ie].id = ie;
    }

    return;
}
