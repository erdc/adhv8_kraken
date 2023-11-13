/* unrefines the mesh */

#include "global_header.h"

/* Gajanan gkc - shortening these functions while adding support for prisms. */

void adpt_set_flags(
                    SGRID *grid,                /* the model grid */
                    int *node_unref_flags		/* flags nodes to be eliminated */
)
{
    int i,j;			/* loop counters over the nodes */
    int ie;			/* loop counter over the elements */
    int elem_level;		/* the level of the element */
#ifdef _MESSG
    int ibrdr;			/* counter for number of owned nodes staying */
#endif
    
    /* initialize flags to YES */
    for(i = 0; i < grid->nnodes; i++) {
        node_unref_flags[i] = YES;
    }
    
    
    if (grid->ndim == 3) {
        /* loop over the 3D elements */
        for(ie = 0; ie < grid->nelems3d; ie++) {
            if (grid->elem3d[ie].mat != UNSET_INT) {
                /* for elements with large errors or nodes that can't be removed for geometric reasons, set node_unref_flags to NO */
                elem_level = elem3d_level(grid->elem3d[ie]);
                if(grid->elem_error[ie] > UNREF_TOL || elem_level == 0) {
                    for (i=0; i<grid->elem3d[ie].nnodes; i++)
                        node_unref_flags[grid->elem3d[ie].nodes[i]] = NO;
                } else {
                    for (i=0; i<grid->elem3d[ie].nnodes; i++){
                        if(grid->elem3d[ie].levels[i] < elem_level)
                            node_unref_flags[grid->elem3d[ie].nodes[i]] = NO;
                    }
                }
            }
//            if (grid->elem3d[ie].string != UNSET_INT && grid->elem_error[ie] < UNREF_TOL && elem3d_level(grid->elem3d[ie]) > 0)
//                printf("4D ie: %d mat: %d \t error: %20.10e level: %d \t UNREF_TOL: %20.10e \t node levls: %d %d %d %d \t node_unref_flags: %d %d %d %d\n",
//                   ie,grid->elem3d[ie].string,grid->elem_error[ie],elem_level,UNREF_TOL,
//                       grid->elem3d[ie].levels[0],
//                       grid->elem3d[ie].levels[1],
//                       grid->elem3d[ie].levels[2],
//                       grid->elem3d[ie].levels[3],
//                       node_unref_flags[grid->elem3d[ie].nodes[0]],
//                       node_unref_flags[grid->elem3d[ie].nodes[1]],
//                       node_unref_flags[grid->elem3d[ie].nodes[2]],
//                       node_unref_flags[grid->elem3d[ie].nodes[3]]);
        }
        //for (i=0; i<grid->nnodes; i++) printf("node_unref_flags[%d]: %d\n",i,node_unref_flags[i]);
        
    } else if (grid->ndim == 2) {
        /* loop over the 2D elements */
        for(ie = 0; ie < grid->nelems2d; ie++) {
            if(grid->elem2d[ie].mat != UNSET_INT) {
                /* for elements with large errors or nodes that can't be removed for geometric reasons, set node_unref_flags to NO */
                elem_level = elem2d_level(grid->elem2d[ie]);
                if(grid->elem_error[ie] > UNREF_TOL || elem_level == 0) {
                    for (i=0; i<grid->elem2d[ie].nnodes; i++)
                        node_unref_flags[grid->elem2d[ie].nodes[i]] = NO;
                } else {
                    for (i=0; i<grid->elem2d[ie].nnodes; i++){
                        if(grid->elem2d[ie].levels[i] < elem_level)
                            node_unref_flags[grid->elem2d[ie].nodes[i]] = NO;
                    }
                }
            }
        }
        
    }  else if (grid->ndim == 1) {
        /* loop over the 1D elements */
        for(ie = 0; ie < grid->nelems1d; ie++) {
            if(grid->elem1d[ie].string != UNSET_INT) {
                /* for elements with large errors or nodes that can't be removed for geometric reasons, set node_unref_flags to NO */
                elem_level = elem1d_level(grid->elem1d[ie]);
                if(grid->elem_error[ie] > UNREF_TOL || elem_level == 0) {
                    for (i=0; i<grid->elem1d[ie].nnodes; i++)
                        node_unref_flags[grid->elem1d[ie].nodes[i]] = NO;
                } else {
                    for (i=0; i<grid->elem1d[ie].nnodes; i++){
                        if(grid->elem1d[ie].levels[i] < elem_level)
                            node_unref_flags[grid->elem1d[ie].nodes[i]] = NO;
                    }
                }
            }
        }
    }
    
    
    
#ifdef _MESSG
    ibrdr = 0;
    for(i = 0; i < grid->my_nnodes; i++)
        if(node_unref_flags[i] == NO)
            ibrdr++;
    if(ibrdr == 0)		/* if no owned nodes are staying, pick a border node to stay */
        for(i = 0; i < grid->my_nnodes; i++)
            if(grid->node[i].parent_res_pe[0] != grid->smpi->myid || grid->node[i].parent_res_pe[1] != grid->smpi->myid)
            {
                node_unref_flags[i] = NO;
                goto one_brdr_staying;
            }
    if(ibrdr == 0)
        tl_error("All owned nodes are set to be removed from processor in adpt_set_flags.");
one_brdr_staying:;
    
    /* update the node_unref_flags */
    comm_update_int(node_unref_flags, 1, grid->smpi);
#endif
}
