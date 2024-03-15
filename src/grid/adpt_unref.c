/* Unrefines the mesh never coarsening past the initial grid. */
/* JRC: add functionality to conserve mass when unrefined nodes exist */
#include "global_header.h"

#ifdef _MESSG
static ELEM_REF_LIST_ITEM **send_3d_elem = NULL;	/* linked list of 3D elements to send */
static ELEM_REF_LIST_ITEM **send_2d_elem = NULL;	/* linked list of 2D elements to send */
static ELEM_REF_LIST_ITEM **send_1d_elem = NULL;	/* linked list of 1D elements to send */
#endif

void adpt_unref(SMODEL *mod) {
    
    assert(mod->grid->type != COLUMNAR);
    
    int check;			/* flag for reseting matrix */
    int i;                /* loop counter */
    int ie;               /* loop counter */
    int *node_unref_flags;/* flags nodes to be eliminated */
    
#ifdef _MESSG
    int isd;			/* loop counter for the processor number */
    int new_size_node_unref_flags;	/* new size for realloc of node_unref_flags */
    NODE_LIST_ITEM *node_hashtab[HASHSIZE];	/* node hash table */
#endif
    
    ELEM1D_LIST_ITEM *elem1d_hashtab[HASHSIZE];	/* 1D element hash table */
    ELEM2D_LIST_ITEM *elem2d_hashtab[HASHSIZE];	/* 2D element hash table */
    ELEM3D_LIST_ITEM *elem3d_hashtab[HASHSIZE];	/* 3D element hash table */
    
    /* set check = nnode to determine whether nodes were removed */
    check = mod->grid->nnodes;
    
#ifdef _MESSG
    if(send_3d_elem == NULL) send_3d_elem = (ELEM_REF_LIST_ITEM **) tl_alloc(sizeof(ELEM_REF_LIST_ITEM *), mod->grid->smpi->npes);
    if(send_2d_elem == NULL) send_2d_elem = (ELEM_REF_LIST_ITEM **) tl_alloc(sizeof(ELEM_REF_LIST_ITEM *), mod->grid->smpi->npes);
    if(send_1d_elem == NULL) send_1d_elem = (ELEM_REF_LIST_ITEM **) tl_alloc(sizeof(ELEM_REF_LIST_ITEM *), mod->grid->smpi->npes);
    for(isd = 0; isd < mod->grid->smpi->npes; isd++) {
        send_3d_elem[isd] = NULL;
        send_2d_elem[isd] = NULL;
        send_1d_elem[isd] = NULL;
    }
#endif
    
    /* initialize hash table */
    for(i = 0; i < HASHSIZE; i++) {
#ifdef _MESSG
        node_hashtab[i] = NULL;
#endif
        elem3d_hashtab[i] = NULL;
        elem2d_hashtab[i] = NULL;
        elem1d_hashtab[i] = NULL;
    }
    
    /* fill hash tables */
#ifdef _MESSG
    for(i = 0; i < mod->grid->nnodes; i++)
        node_hash_add_entry(mod->grid->node[i], node_hashtab, i, mod->grid->smpi->npes);
#endif
    for(ie = 0; ie < mod->grid->nelems3d; ie++) {
        elem3d_hash_add_entry(mod->grid->elem3d[ie].nodes[0],
                              mod->grid->elem3d[ie].nodes[1],
                              mod->grid->elem3d[ie].nodes[2],
                              mod->grid->elem3d[ie].nodes[3], elem3d_hashtab, ie);
    }
    for(ie = 0; ie < mod->grid->nelems2d; ie++) {
        elem2d_hash_add_entry(mod->grid->elem2d[ie].nodes[0],
                              mod->grid->elem2d[ie].nodes[1],
                              mod->grid->elem2d[ie].nodes[2], elem2d_hashtab, ie);
    }
    for(ie = 0; ie < mod->grid->nelems1d; ie++) {
        elem1d_hash_add_entry(mod->grid->elem1d[ie].nodes[0],
                              mod->grid->elem1d[ie].nodes[1], elem1d_hashtab, ie);
    }
    
    /* allocate unref flags */
    node_unref_flags = (int *)tl_alloc(sizeof(int), mod->grid->nnodes);
    
    /* set the unrefinement flags */
    adpt_set_flags(mod->grid, node_unref_flags);
    
#ifdef _MESSG
    int old_nnodes = mod->grid->nnodes;
    int num_adj_nodes_added;
    
    comm_adj_nodes(node_unref_flags, node_hashtab, mod);
    
    //comm_update_snode(grid);
    //set_parents_local(grid);
    
    new_size_node_unref_flags = mod->grid->nnodes;
    num_adj_nodes_added = new_size_node_unref_flags - old_nnodes;
    if(num_adj_nodes_added != 0) {
        node_unref_flags = (int *)tl_realloc(sizeof(int), new_size_node_unref_flags, old_nnodes, node_unref_flags);
    }
#endif
    
    /* merge the elements */
    adpt_merge_elem (mod->grid, node_unref_flags,
#ifdef _MESSG
                     send_3d_elem, send_2d_elem, send_1d_elem, node_hashtab,
#endif
                     elem3d_hashtab, elem2d_hashtab, elem1d_hashtab);
    
    /* communicate element-wise node levels */
#ifdef _MESSG
    /* post the receives */
    comm_elem_levels(send_3d_elem, send_2d_elem, send_1d_elem, node_hashtab, elem3d_hashtab, elem2d_hashtab, elem1d_hashtab, mod->grid);
#endif
    
    //print_grid_to_file(grid,"POST_ELEM_LEVELS");int elemerror=0;
    
    /* error traps */
    for(ie = 0; ie < mod->grid->nelems3d; ie++)
        if(mod->grid->elem3d[ie].string != UNSET_INT) {
            if(mod->grid->elem3d[ie].levels[0] == UNSET_INT || mod->grid->elem3d[ie].levels[1] == UNSET_INT ||
               mod->grid->elem3d[ie].levels[2] == UNSET_INT || mod->grid->elem3d[ie].levels[3] == UNSET_INT) {
                tl_error("3D element levels are not set in adpt_unref.");
            }
            if(node_unref_flags[mod->grid->elem3d[ie].nodes[0]] == YES || node_unref_flags[mod->grid->elem3d[ie].nodes[1]] == YES ||
               node_unref_flags[mod->grid->elem3d[ie].nodes[2]] == YES || node_unref_flags[mod->grid->elem3d[ie].nodes[3]] == YES) {
                tl_error("3D element nodes being pulled are still in use in adpt_unref.");
            }
        }
    for(ie = 0; ie < mod->grid->nelems2d; ie++)
        if(mod->grid->elem2d[ie].mat != UNSET_INT) {
            if(mod->grid->elem2d[ie].levels[0] == UNSET_INT ||
               mod->grid->elem2d[ie].levels[1] == UNSET_INT ||
               mod->grid->elem2d[ie].levels[2] == UNSET_INT) {
                tl_error("2D element levels are not set in adpt_unref.");
            }
            if(node_unref_flags[mod->grid->elem2d[ie].nodes[0]] == YES ||
               node_unref_flags[mod->grid->elem2d[ie].nodes[1]] == YES ||
               node_unref_flags[mod->grid->elem2d[ie].nodes[2]] == YES) {
                tl_error("2D element nodes being pulled are still in use in adpt_unref.");
            }
        }
    for(ie = 0; ie < mod->grid->nelems1d; ie++)
        if(mod->grid->elem1d[ie].string != UNSET_INT) {
            if(mod->grid->elem1d[ie].levels[0] == UNSET_INT || mod->grid->elem1d[ie].levels[1] == UNSET_INT) {
                tl_error("1D element levels are not set in adpt_unref.");
            }
            if(node_unref_flags[mod->grid->elem1d[ie].nodes[0]] == YES ||
               node_unref_flags[mod->grid->elem1d[ie].nodes[1]] == YES) {
                tl_error("1D element nodes being pulled are still in use in adpt_unref.");
            }
        }
    
#ifdef _MESSG
    /* fix the adjacent nodes */
    adpt_fix_adj(mod->grid, node_unref_flags, node_hashtab);
    
    /* clean up unneeded nodes and elements */
    //  partition_cleanup();
#else
    adpt_fix_adj(mod->grid, node_unref_flags);
#endif
    
    
    //print_grid_to_file(grid,"POST_UNREF_FIX");
#ifdef _MESSG
    partition_cleanup(mod->grid, mod->grid->smpi->myid);
#endif
    
    //print_grid_to_file(grid,"POST_UNREF_CLEAN");
    /* renumbers the mesh */
    node_renumber(mod, 0);//print_grid_to_file(grid,"POST_UNREF_NODE");
    elem3d_renumber(mod->grid);
    elem2d_renumber(mod->grid);//print_grid_to_file(grid,"POST_UNREF_ELEM2D");
    elem1d_renumber(mod->grid);//print_grid_to_file(grid,"POST_UNREF_ELEM1D");
    adpt_fix_global_ids(mod->grid);//print_grid_to_file(grid,"POST_UNREF_GID");
    if (mod->grid->ndim == 3) {
        classify_2d_elements(&mod->grid);
    }
    
    
    //print_grid_to_file(grid,"POST_UNREF_RENUM");
    
#ifdef _MESSG
    /* updates the message keys for the repartitioned mesh */
    comm_set_keys(mod->grid);
    
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
    if(mod->grid->nelems1d != 0) {
        tl_list_free_all(ELEM1D_LIST);
    }
    if(mod->grid->nelems2d != 0) {
        tl_list_free_all(ELEM2D_LIST);
    }
    if(mod->grid->nelems3d != 0) {
        tl_list_free_all(ELEM3D_LIST);
    }
    
#ifdef _MESSG
    node_unref_flags = (int *)tl_free(sizeof(int), new_size_node_unref_flags, node_unref_flags);
    tl_list_free_all(NODE_LIST);
    if(send_3d_elem != NULL)
        send_3d_elem = (ELEM_REF_LIST_ITEM * *)tl_free(sizeof(ELEM_REF_LIST_ITEM *), mod->grid->smpi->npes, send_3d_elem);
    if(send_2d_elem != NULL)
        send_2d_elem = (ELEM_REF_LIST_ITEM * *)tl_free(sizeof(ELEM_REF_LIST_ITEM *), mod->grid->smpi->npes, send_2d_elem);
    if(send_1d_elem != NULL)
        send_1d_elem = (ELEM_REF_LIST_ITEM * *)tl_free(sizeof(ELEM_REF_LIST_ITEM *), mod->grid->smpi->npes, send_1d_elem);
#else
    node_unref_flags = (int *)tl_free(sizeof(int), check, node_unref_flags);
#endif
    
}

//*******************************************************************************************************//
//*******************************************************************************************************//
//*******************************************************************************************************//
//*******************************************************************************************************//

void adpt_unref_clean() {
#ifdef _MESSG_ADPT
    if(send_3d_elem != NULL)
        send_3d_elem = (ELEM_REF_LIST_ITEM * *)tl_free(sizeof(ELEM_REF_LIST_ITEM *), npes, send_3d_elem);
    if(send_2d_elem != NULL)
        send_2d_elem = (ELEM_REF_LIST_ITEM * *)tl_free(sizeof(ELEM_REF_LIST_ITEM *), npes, send_2d_elem);
    if(send_1d_elem != NULL)
        send_1d_elem = (ELEM_REF_LIST_ITEM * *)tl_free(sizeof(ELEM_REF_LIST_ITEM *), npes, send_1d_elem);
#endif
}
