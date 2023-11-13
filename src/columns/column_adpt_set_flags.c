/* ADH Version 2.0.0 6-04 */

#include "global_header.h"

// flags nodes to be eliminated
void column_adpt_set_flags(SGRID *grid, int *surf_unref_flags, int *node_unref_flags) {
    
    int i;                        /* loop counters over the nodes */
    int ie;                       /* loop counter over the elements */
    ID_LIST_ITEM *ptr;
#ifdef _MESSG
    int ibrdr;                    /* counter for number of owned nodes staying */
#endif
    
    /* initialize flags to NO */
    sarray_init_value_int(node_unref_flags, grid->nnodes, NO);
    
    /* loop over the surface nodes */
    for (i = 0; i < grid->nnodes_sur; i++) {
        if (surf_unref_flags[i] == YES) {
            ptr = grid->vertical_list[i];
            while (ptr->next != NULL) {
                node_unref_flags[ptr->id] = YES;
                ptr = ptr->next;
            }
        }
    }
    
#ifdef _MESSG // cjt :: I'm not certain this is working like it should ... not sure who added it
    ibrdr = 0;
    for (i = 0; i < grid->my_nnodes; i++)
        if (node_unref_flags[i] == NO)
            ibrdr++;
    if (ibrdr == 0)   /* if no owned nodes are staying, pick a border node to stay */
        for (i = 0; i < grid->my_nnodes; i++)
            if (grid->node[i].parent_res_pe[0] != grid->smpi->myid || grid->node[i].parent_res_pe[1] != grid->smpi->myid) {
                node_unref_flags[i] = NO;
                goto one_brdr_staying;
            }
    if (ibrdr == 0)
        tl_error("All owned nodes are set to be removed from processor in adpt_set_flags.");
one_brdr_staying:;
    
    /* update the node_unref_flags */
    comm_update_int(node_unref_flags, 1, grid->smpi);
#endif
}
