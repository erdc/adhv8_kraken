/* this routine computes the new order for the nodes
 - should add bandwidth renumbering here
 - metis has fill reducing ordering */

#include "global_header.h"
typedef int (*compfn)(const void*, const void*);

void node_order(SGRID *grid, int nblock, int *new_node_number) {
    
    int i;                        /* loop counter */
    int iblk;                     /* loop counter over the blocks */
    int blk_inc;                  /* the number of nodes in each block */
    int ifirst, ilast;            /* bounds on the loops to load the blocks */
    int inode;                    /* loop counter through the nodes */
    int my_block;                 /* the block number being assigned */
    int iblk_shift;               /* the block shift for parallel calculations */
    
    /* allocate and initialize the sort array */
    SNODE *tmp_node = (SNODE *) tl_alloc(sizeof(SNODE),grid->max_nnodes);
    
    for (i = 0; i < grid->max_nnodes; i++) {
        grid->node[i].myid = grid->smpi->myid;
        snode_copy(&(tmp_node[i]), grid->node[i]);
    }


    /* sorts the nodes NOTE:  WILL NOT CHANGE IDS */
    qsort((void *) tmp_node, grid->max_nnodes, sizeof(SNODE), (compfn)node_cmp);
  
    /* sets the new node numbers */
    for (i = 0; i < grid->max_nnodes; i++) {
        if((tmp_node[i].id >= 0) && (tmp_node[i].id < grid->max_nnodes)){
          
          new_node_number[tmp_node[i].id] = i;
        }
        else{
          new_node_number[i]=i;
        }
         
    }

    /* fix the node numbers */
    for (i = 0, grid->my_nnodes = 0, grid->nnodes = 0; i < grid->max_nnodes; i++) {
        
        if (grid->node[i].string != UNSET_INT) {
            grid->nnodes++;
#ifdef _MESSG
            if (grid->node[i].resident_pe == grid->smpi->myid) {
                grid->my_nnodes++;
            }
#else
            grid->my_nnodes++;
#endif
        }
    }
   /* fix 2D surafce and bed nodes, for 3D this is in the columns routine*/
    if(grid->ndim==2){
      grid->my_nnodes_sur=grid->my_nnodes;
      grid->my_nnodes_bed=grid->my_nnodes;
      grid->nnodes_sur=grid->nnodes;
      grid->nnodes_bed=grid->nnodes;
    }
    /* calculate the block increment */
    blk_inc = grid->my_nnodes / nblock;
    if (blk_inc < 1) {
        tl_error("Too many blocks for the number of equations in node_order.");
    }
    
    /* calculate the block shift for parallel calculations -
     NOTE:  this is tied into the shift in solve_coarse_load */
#ifdef _MESSG
    iblk_shift = grid->smpi->myid * nblock;
   
#else
    iblk_shift = 0;
#endif
    
    /* split the nodes into blocks - this should be done with Metis - weights for better performance */
    for (iblk = 0, ilast = 0, inode = 0; iblk < nblock; iblk++) {
        /* sets the first and last nodes of the subdomain - this assignment should be done with metis or some such */
        ifirst = ilast;
        if (iblk == nblock - 1)
            ilast = grid->my_nnodes;
        else
            ilast += blk_inc;
        
        /* sets the block */
        my_block = iblk_shift + iblk;
        
        /* sets the block number */
        for (i = ifirst; i < ilast; inode++) {
            if (grid->node[inode].string != UNSET_INT) {
#ifdef _MESSG
                if (grid->node[inode].resident_pe == grid->smpi->myid) {
                    i++;
                    grid->node[inode].block = my_block;
                }
#else
                i++;
                grid->node[inode].block = my_block;
#endif
            }
        }
    }
    
    /* free the memory */
    tmp_node = (SNODE *) tl_free(sizeof(SNODE), grid->max_nnodes, tmp_node);
}
