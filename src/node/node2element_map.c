#include "global_header.h"

void node2elem2d_map(SGRID *grid) {
    int ie=0, inode=0, gnode=0;
    
    int DEBUG = OFF;
    
    // reset
    for (inode=0; inode<grid->nnodes; inode++) {
	if (grid->node[inode].elemID != NULL) grid->node[inode].elemID  = (int *) tl_free(sizeof(int), grid->node[inode].nelems_connected, grid->node[inode].elemID);
	grid->node[inode].nelems_connected = 0;
    }

     // count the total number of elements connected to each node
    for (ie=0; ie<grid->nelems2d; ie++) {
        for (inode=0; inode<NDONTRI; inode++) {
            gnode = grid->elem2d[ie].nodes[inode];
            (grid->node[gnode].nelems_connected)++;
        }
    }
    
    // allocate elemIDs
    for (inode=0; inode<grid->nnodes; inode++) {
	grid->node[inode].elemID = (int *)tl_alloc(sizeof(int), grid->node[inode].nelems_connected);
    }
    
    // store element IDs
    int *elem_count = (int *)tl_alloc(sizeof(int), grid->nnodes);
    sarray_init_int(elem_count, grid->nnodes);
    for (ie=0; ie<grid->nelems2d; ie++) {
        
        for (inode=0; inode<NDONTRI; inode++) {
            gnode = grid->elem2d[ie].nodes[inode];
            grid->node[gnode].elemID[elem_count[gnode]] = ie;
            elem_count[gnode]++;
        }
    }
    
    elem_count = (int *) tl_free(sizeof(int), grid->nnodes, elem_count);

#ifdef _DEBUG
    if (DEBUG) {
        for (inode=0; inode<grid->nnodes; inode++) {
            printf("node: %d eIDs ::",inode);
            for (ie=0; ie<grid->node[inode].nelems_connected; ie++) {
                printf(" %d \t",grid->node[inode].elemID[ie]);
            }
            printf("\n");
        }
    }
#endif
    
}

