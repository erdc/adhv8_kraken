/* ADH Version 2.0.0 6-04 */
/* Unrefines the columnar mesh never coarsening past the initial grid. */
/* The key here is to use the 2D surface mesh as a guide to unrefinement */

#include "global_header.h"

#ifdef _MESSG
static ELEM_REF_LIST_ITEM **send_3d_elem = NULL;    /* linked list of 3D elements to send */
static ELEM_REF_LIST_ITEM **send_2d_elem = NULL;    /* linked list of 2D elements to send */
static ELEM_REF_LIST_ITEM **send_1d_elem = NULL;    /* linked list of 1D elements to send */
#endif

void column_adpt_unref(SMODEL *mod) {
    
    int check;                    /* flag for reseting matrix */
    int i, jj;                    /* loop counter */
    int ie;
    int *node_unref_flags;        /* flags nodes to be eliminated */
    int ielem;
    int inode;
    int current_level, hi_level;
    int new_node;
    ID_LIST_ITEM *ptr;
    int *surf_unref_flags;        /* flag columns to be eliminated */
    int elem_level;
    int new_local;
    int column_id[3];             /* the columns that the 3 nodes of an element are in */
    
    // save the original nnodes and max_nnodes for physics reallocation
    assert(mod->grid->nnodes_sur == mod->grid->nnodes_bed); // sanity
    int old_nnodes = mod->grid->nnodes;
    int old_max_nnodes = mod->grid->max_nnodes;
    int old_nnodes2d = mod->grid->nnodes_sur;
    int old_max_nnodes2d = mod->grid->max_nnodes_sur;
    int old_nnodes_sur = mod->grid->nnodes_sur;
    ELEM1D_LIST_ITEM *elem1d_hashtab[HASHSIZE];   /* 1D element hash table */
    ELEM2D_LIST_ITEM *elem2d_hashtab[HASHSIZE];   /* 2D element hash table */
    ELEM3D_LIST_ITEM *elem3d_hashtab[HASHSIZE];   /* 3D element hash table */
    
#ifdef _MESSG
    int isd;                      /* loop counter for the processor number */
    int new_size_node_unref_flags;    /* new size for realloc of node_unref_flags */
    NODE_LIST_ITEM *node_hashtab[HASHSIZE];   /* node hash table */
#else
    NODE_LIST_ITEM *node_hashtab = NULL;
#endif
    
    /* set check = nnode to determine whether nodes were removed */
    check = mod->grid->nnodes;
    
#ifdef _MESSG
    if (send_3d_elem == NULL) send_3d_elem = (ELEM_REF_LIST_ITEM **) tl_alloc(sizeof(ELEM_REF_LIST_ITEM *), mod->grid->smpi->npes);
    if (send_2d_elem == NULL) send_2d_elem = (ELEM_REF_LIST_ITEM **) tl_alloc(sizeof(ELEM_REF_LIST_ITEM *), mod->grid->smpi->npes);
    if (send_1d_elem == NULL) send_1d_elem = (ELEM_REF_LIST_ITEM **) tl_alloc(sizeof(ELEM_REF_LIST_ITEM *), mod->grid->smpi->npes);
    for (isd = 0; isd < mod->grid->smpi->npes; isd++) {
        send_3d_elem[isd] = NULL;
        send_2d_elem[isd] = NULL;
        send_1d_elem[isd] = NULL;
    }
#endif
    
    /* initialize hash table */
    for (i = 0; i < HASHSIZE; i++) {
#ifdef _MESSG
        node_hashtab[i] = NULL;
#endif
        elem3d_hashtab[i] = NULL;
        elem2d_hashtab[i] = NULL;
        elem1d_hashtab[i] = NULL;
    }
    
    /* fill hash tables */
#ifdef _MESSG
    for (i = 0; i < mod->grid->nnodes; i++)
        node_hash_add_entry(mod->grid->node[i], node_hashtab, i, mod->grid->smpi->npes);
#endif
    for (ie = 0; ie < mod->grid->nelems3d; ie++) {
        elem3d_hash_add_entry(mod->grid->elem3d[ie].nodes[0],
                              mod->grid->elem3d[ie].nodes[1],
                              mod->grid->elem3d[ie].nodes[2],
                              mod->grid->elem3d[ie].nodes[3], elem3d_hashtab, ie);
    }
    for (ie = 0; ie < mod->grid->nelems2d; ie++) {
        elem2d_hash_add_entry(mod->grid->elem2d[ie].nodes[0],
                              mod->grid->elem2d[ie].nodes[1],
                              mod->grid->elem2d[ie].nodes[2], elem2d_hashtab, ie);
    }
    for (ie = 0; ie < mod->grid->nelems1d; ie++) {
        elem1d_hash_add_entry(mod->grid->elem1d[ie].nodes[0],
                              mod->grid->elem1d[ie].nodes[1], elem1d_hashtab, ie);
    }
    
    /* allocate unref flags */
    node_unref_flags = (int *) tl_alloc(sizeof(int), mod->grid->nnodes);
    
    /* surf_unref_flags indicates which surface nodes are to be removed.  */
    surf_unref_flags = (int *) tl_alloc(sizeof(int), mod->grid->nnodes_sur);
    int nsurf_unref_flags = mod->grid->nnodes_sur;
    
    /* set flags which inidicate which surface nodes can be removed.
     * We key on the (2D) surface elements and we don't need to consult the levels on
     * the 3D elements at all
     */
    
    for (i = 0; i < mod->grid->nnodes_sur; i++) {
        surf_unref_flags[i] = YES;
    }
    
#ifdef _MESSG
    /***************************************************/
    /***************************************************/
    /* cjt :: 2_4_2016 :: This was done to to avoid unrefining
     to zero residential nodes.  It says that if the domain
     has less than 20 surface residential nodes, do not unrefine.
     This  is an ad hoc number though, and some cases could occur
     where is should be increased */
    int min_node = 20;
    int icol, nd;
    int resid_surf_count = 0;
    for (icol = 0; icol < mod->grid->nnodes_sur; icol++) {
        ptr = mod->grid->vertical_list[icol];
        nd = ptr->id;
        if (mod->grid->node[nd].resident_pe == mod->grid->smpi->myid) resid_surf_count++;
    }
    if (resid_surf_count < min_node) {
        for (i = 0; i < mod->grid->nnodes_sur; i++) {
            surf_unref_flags[i] = NO;
        }
    }
    /***************************************************/
    /***************************************************/
#endif
    
    /* Determine which surface nodes can be removed */
    /* Loops over columns */
    for (i = 0; i < mod->grid->ncolumns; i++) {
        ie = mod->grid->elem2d_sur[i];
        elem_level = elem2d_level(mod->grid->elem2d[ie]);
        for (jj = 0; jj < 3; jj++) {
            column_id[jj] = find_vertical_segment(mod->grid, mod->grid->elem2d[ie].nodes[jj], mod->grid->vertical_hash);
        }
        
        /* If original element, do not unrefine */
        if (elem_level == 0) {
            for (jj = 0; jj < 3; jj++) {
                surf_unref_flags[column_id[jj]] = NO;
            }
            continue;
        }
        else {
            for (jj=0; jj<mod->grid->elem2d[ie].nnodes; jj++){
                if (mod->grid->elem2d[ie].levels[jj] < elem_level)
                    surf_unref_flags[column_id[jj]] = NO;
            }
        }
        
        /* check error tolerance in entire column
         * If the error in any of the elements is too large, then
         * we do not unrefine anything in the column */
        ptr = mod->grid->column_list[i];
        while (ptr->next != NULL) {
            ielem = ptr->id;
            /* Gajanan gkc question / doubt - shouldn't this be UNREF_TOL and not 0.9? */
            if (mod->grid->elem_error[ielem] > UNREF_TOL) {
                //if (mod->grid->elem_error[ielem] > 0.9) {
                for (jj=0; jj<mod->grid->elem2d[ie].nnodes; jj++){
                    surf_unref_flags[column_id[jj]] = NO;
                }
                break;
            }
            ptr = ptr->next;
        }
    }
    
    /*   for (i = 0; i < mod->grid->nnodes_sur; i++) {
     if(surf_unref_flags[i] == YES)printf("MYID %d grid UNREFINING!\n");
     }
     */
    /* set the unref flags for all nodes below surface nodes that are marked to be removed. */
    column_adpt_set_flags(mod->grid, surf_unref_flags, node_unref_flags);
    
#ifdef _MESSG
    
    comm_adj_nodes(node_unref_flags, node_hashtab, mod);
    
    
    int num_adj_nodes_added;
    if(mod->grid->nnodes_sur != old_nnodes_sur){
        surf_unref_flags = (int *) tl_realloc(sizeof(int), mod->grid->nnodes_sur, old_nnodes_sur, surf_unref_flags);
        old_nnodes_sur = mod->grid->nnodes_sur;
    }
    /* Update surf_unref_flags by peeking at node_unref_flags */
    for (i = 0; i < mod->grid->nnodes_sur; i++) {
        
        if (node_unref_flags[mod->grid->nodeID_2d_to_3d_sur[i]] == NO) {
            surf_unref_flags[i] = NO;
        }
        else {
            surf_unref_flags[i] = YES;
        }
    }
    new_size_node_unref_flags = mod->grid->nnodes;
    num_adj_nodes_added = new_size_node_unref_flags - old_nnodes;
    if (num_adj_nodes_added != 0)
        node_unref_flags = (int *) tl_realloc(sizeof(int), new_size_node_unref_flags, old_nnodes, node_unref_flags);
#endif
    //for(i=0;i<mod->grid->nnodes;i++)if (node_unref_flags[i] == YES) printf("MYID %d removing node %d gid %d x y z %10.5e %10.5e %10.5e \n", mod->grid->smpi->myid, i, mod->grid->node[i].gid,mod->grid->node[i].x,mod->grid->node[i].y,mod->grid->node[i].z);
    /* Now, determine which columns can be unrefined */
    
    for (i = 0; i < mod->grid->ncolumns; i++) {
        ie = mod->grid->elem2d_sur[i];
        for (jj = 0; jj < 3; jj++) {
            column_id[jj] = find_vertical_segment(mod->grid, mod->grid->elem2d[ie].nodes[jj], mod->grid->vertical_hash);
        }
        
        if (surf_unref_flags[column_id[0]] == YES || surf_unref_flags[column_id[1]] == YES || surf_unref_flags[column_id[2]] == YES) {
            /* Unrefine the column */
            /* First, find which node is to be removed */
            hi_level = UNSET_INT;
            for (jj=0; jj<mod->grid->elem2d[ie].nnodes; jj++){
                current_level = mod->grid->elem2d[ie].levels[jj];
                if (current_level > hi_level) {
                    hi_level = current_level;
                    new_local = jj;
                }
            }
            /* Gajanan gkc possible bugfix - I think this should have been hi_level... Please crosscheck, someone! */
            if (hi_level == 0) {
                //if (current_level == UNSET_INT) {
                tl_error("Tried to unref original element in column_adpt_unref.");
            }
            
            /* new_node is the node that is going to be removed */
            new_node = mod->grid->elem2d[ie].nodes[new_local];
#ifdef _MESSG
            column_adpt_merge_elem(mod->grid, i, new_node, node_unref_flags, send_2d_elem, node_hashtab, elem3d_hashtab, elem2d_hashtab);
#else
            column_adpt_merge_elem(mod->grid, i, new_node, node_unref_flags, NULL, NULL, elem3d_hashtab, elem2d_hashtab);
#endif
        }
    }
    if(mod->grid->nnodes_sur != old_nnodes_sur){
        surf_unref_flags = (int *) tl_realloc(sizeof(int), mod->grid->nnodes_sur, old_nnodes_sur, surf_unref_flags);
        old_nnodes_sur = mod->grid->nnodes_sur;
    }
#ifdef _MESSG
    comm_elem_levels(send_3d_elem, send_2d_elem, send_1d_elem, node_hashtab,elem3d_hashtab,  elem2d_hashtab, elem1d_hashtab, mod->grid);
#endif
    
    if(mod->grid->nnodes_sur != old_nnodes_sur){
        surf_unref_flags = (int *) tl_realloc(sizeof(int), mod->grid->nnodes_sur, old_nnodes_sur, surf_unref_flags);
        old_nnodes_sur = mod->grid->nnodes_sur;
    }
    /* fix the adjacent nodes */
#ifdef _MESSG
    adpt_fix_adj(mod->grid, node_unref_flags, node_hashtab);
    if(mod->grid->nnodes_sur != old_nnodes_sur){
        surf_unref_flags = (int *) tl_realloc(sizeof(int), mod->grid->nnodes_sur, old_nnodes_sur, surf_unref_flags);
        old_nnodes_sur = mod->grid->nnodes_sur;
    }
    partition_cleanup(mod->grid, mod->grid->smpi->myid);
    if(mod->grid->nnodes_sur != old_nnodes_sur){
        surf_unref_flags = (int *) tl_realloc(sizeof(int), mod->grid->nnodes_sur, old_nnodes_sur, surf_unref_flags);
        old_nnodes_sur = mod->grid->nnodes_sur;
    }
#else
    adpt_fix_adj(mod->grid, node_unref_flags);
#endif
    
    node_renumber(mod, 0);
    elem3d_renumber(mod->grid);
    elem2d_renumber(mod->grid);
    elem1d_renumber(mod->grid);
    if(mod->grid->nnodes_sur != old_nnodes_sur){
        surf_unref_flags = (int *) tl_realloc(sizeof(int), mod->grid->nnodes_sur, old_nnodes_sur, surf_unref_flags);
        old_nnodes_sur = mod->grid->nnodes_sur;
    }
    
    /* flag to reset the matrix if nodes were removed */
    if(check != mod->grid->nnodes) {
        mod->flag.UNREFINE = YES;
    } else {
        mod->flag.UNREFINE = NO;
    }
    
    /* re-build columns */
    build_columns(mod->grid, NO);
    if(mod->grid->nnodes_sur != old_nnodes_sur){
        surf_unref_flags = (int *) tl_realloc(sizeof(int), mod->grid->nnodes_sur, old_nnodes_sur, surf_unref_flags);
        old_nnodes_sur = mod->grid->nnodes_sur;
    }
    node_renumber_surface(mod);
    if(mod->grid->nnodes_sur != old_nnodes_sur){
        surf_unref_flags = (int *) tl_realloc(sizeof(int), mod->grid->nnodes_sur, old_nnodes_sur, surf_unref_flags);
        old_nnodes_sur = mod->grid->nnodes_sur;
    }
    
    adpt_fix_global_ids(mod->grid);
    if(mod->grid->nnodes_sur != old_nnodes_sur){
        surf_unref_flags = (int *) tl_realloc(sizeof(int), mod->grid->nnodes_sur, old_nnodes_sur, surf_unref_flags);
        old_nnodes_sur = mod->grid->nnodes_sur;
    }
    /* reset maxes for adaption */
    //mod->grid->old_global_bed = tl_realloc(sizeof(int), mod->grid->nnodes_bed, mod->grid->max_nnodes_bed, mod->grid->old_global_bed);
    //mod->grid->old_global_surf = tl_realloc(sizeof(int), mod->grid->nnodes_sur, mod->grid->max_nnodes_sur, mod->grid->old_global_surf);
    
    
#ifdef _MESSG
    /* updates the message keys for the repartitioned mesh */
    comm_set_keys(mod->grid);
    
    /* update ghost nodes and node adjacencies */
    comm_update_GN(1, mod->grid);
    
    /* update ghost nodes and node adjacencies */
    comm_update_snode(mod->grid);
    
#endif
    
    /* flag to reset the matrix if nodes were removed */
    if(check != mod->grid->nnodes) {
        mod->flag.UNREFINE = YES;
    } else {
        mod->flag.UNREFINE = NO;
    }
    
    /* free memory */
    if (surf_unref_flags != NULL) {
        tl_free(sizeof(int), mod->grid->nnodes_sur, (void *) surf_unref_flags);
    }
    
    if (mod->grid->nelems1d != 0)
        tl_list_free_all(ELEM1D_LIST);
    if (mod->grid->nelems2d != 0)
        tl_list_free_all(ELEM2D_LIST);
    if (mod->grid->nelems3d != 0)
        tl_list_free_all(ELEM3D_LIST);
    
#ifdef _MESSG
    tl_free(sizeof(int), new_size_node_unref_flags, (void *) node_unref_flags);
    tl_list_free_all(NODE_LIST);
    if(send_3d_elem != NULL)
        send_3d_elem = (ELEM_REF_LIST_ITEM * *)tl_free(sizeof(ELEM_REF_LIST_ITEM *), mod->grid->smpi->npes, send_3d_elem);
    if(send_2d_elem != NULL)
        send_2d_elem = (ELEM_REF_LIST_ITEM * *)tl_free(sizeof(ELEM_REF_LIST_ITEM *), mod->grid->smpi->npes, send_2d_elem);
    if(send_1d_elem != NULL)
        send_1d_elem = (ELEM_REF_LIST_ITEM * *)tl_free(sizeof(ELEM_REF_LIST_ITEM *), mod->grid->smpi->npes, send_1d_elem);
#else
    tl_free(sizeof(int), check, (void *) node_unref_flags);
#endif
}

/*void column_adpt_unref_clean(SMODEL *) {
 
 int i;
 
 #ifdef _MESSG
 if (send_3d_elem != NULL)
 tl_free(sizeof(ELEM_REF_LIST_ITEM *), mod->grid->smpi->npes, (void *) send_3d_elem);
 if (send_2d_elem != NULL)
 tl_free(sizeof(ELEM_REF_LIST_ITEM *), mod->grid->smpi->npes, (void *) send_2d_elem);
 if (send_1d_elem != NULL)
 tl_free(sizeof(ELEM_REF_LIST_ITEM *), mod->grid->smpi->npes, (void *) send_1d_elem);
 
 #endif
 }
 */
