/* ADH Version 2.0.0 6-04 */
/* refines the mesh */

#include "global_header.h"

// cycle = adaption counter
long column_adpt_ref(SMODEL *mod, long cycle) {
    
    
    int ie;                       /* loop counter over the elements */
    int i;                        /* loop counters */
    int iconform_flag;            /* flag for whether or not a mesh is conforming */
#ifdef _MESSG
    int isd;
    int *nrecv_edge, *nsend_edge;
    int npes = mod->grid->smpi->npes;
#endif
    EDGE_LIST_ITEM *edge_hashtab[HASHSIZE];   /* the hash table of edges - this serves two purposes:
                                               the first is to hold the integer rank of an edge for
                                               length comparisons, and the second is to hold the number
                                               of a node imbedded in the edge.  This allows us to
                                               determine if an edge has been split */
    EDGE_LIST_ITEM *e0, *e1, *e2, *e3, *e4, *e5;  /* the edges of an element */
    
    ID_LIST_ITEM *ptr;
    int refine_col;
    int update_cols;
    int ival;
    int count;
    int nd1, nd2, iedge, nedges;
    int **edges;
    SVECT2D center;
    FILE *fp;
    
    // save the original nnodes and max_nnodes for physics reallocation
    assert(mod->grid->nnodes_sur == mod->grid->nnodes_bed); // sanity
    int old_nnodes = mod->grid->nnodes;
    int old_max_nnodes = mod->grid->max_nnodes;
    int old_nnodes2d = mod->grid->nnodes_sur;
    int old_max_nnodes2d = mod->grid->max_nnodes_sur;
    
    
    /* initialize the hash table */
    for (i = 0; i < HASHSIZE; i++)
        edge_hashtab[i] = NULL;
    
    /* load the edges into the hash table */
    for (ie = 0; ie < mod->grid->nelems3d; ie++) {
        for (iedge = 0; iedge < mod->grid->elem3d[ie].nedges; iedge++) {
            nd1 = mod->grid->elem3d[ie].edges[iedge][0]; nd2 = mod->grid->elem3d[ie].edges[iedge][1];
            edge_hash_add_entry(mod->grid->elem3d[ie].nodes[nd1], mod->grid->elem3d[ie].nodes[nd2], edge_hashtab);
            
        }
        /* Gajanan gkc adding for adaption, to remove edge-not-found error in adpt_rank_edges */
        /* Doing this based on a previous version of AdH that I have */
        for (iedge = 0; iedge < mod->grid->elem3d[ie].nedges; iedge++) {
            nd1 = mod->grid->elem3d[ie].edges[iedge][0]; nd2 = mod->grid->elem3d[ie].edges[iedge][1];
            edge_hash_add_entry(mod->grid->elem3d[ie].nodes[nd2], mod->grid->elem3d[ie].nodes[nd1], edge_hashtab);
            
        }
    }
    for (ie = 0; ie < mod->grid->nelems2d; ie++) {
        for (iedge = 0; iedge < mod->grid->elem2d[ie].nedges; iedge++) {
            nd1 = mod->grid->elem2d[ie].edges[iedge][0]; nd2 = mod->grid->elem2d[ie].edges[iedge][1];
            edge_hash_add_entry(mod->grid->elem2d[ie].nodes[nd1], mod->grid->elem2d[ie].nodes[nd2], edge_hashtab);
            
        }
        /* Gajanan gkc adding for adaption, to remove edge-not-found error in adpt_rank_edges */
        /* Doing this based on a previous version of AdH that I have */
        for (iedge = 0; iedge < mod->grid->elem2d[ie].nedges; iedge++) {
            nd1 = mod->grid->elem2d[ie].edges[iedge][0]; nd2 = mod->grid->elem2d[ie].edges[iedge][1];
            edge_hash_add_entry(mod->grid->elem2d[ie].nodes[nd2], mod->grid->elem2d[ie].nodes[nd1], edge_hashtab);
            
        }
    }
    
    for (ie = 0; ie < mod->grid->nelems1d; ie++) {
        edge_hash_add_entry(mod->grid->elem1d[ie].nodes[0], mod->grid->elem1d[ie].nodes[1], edge_hashtab);
        edge_hash_add_entry(mod->grid->elem1d[ie].nodes[1], mod->grid->elem1d[ie].nodes[0], edge_hashtab);
    }
    
    /* set the edge keys */
#ifdef _MESSG
    nrecv_edge = (int *) tl_alloc(sizeof(int), npes);
    nsend_edge = (int *) tl_alloc(sizeof(int), npes);
    comm_edge_setup(edge_hashtab, mod->grid, nrecv_edge, nsend_edge);
#endif
    
    /* rank the edges */
#ifdef _MESSG
    adpt_rank_edges(edge_hashtab, mod->grid, nrecv_edge, nsend_edge);
#else
    adpt_rank_edges(edge_hashtab, mod->grid);
#endif
    
    /* Add resolution */
    update_cols = 0;  /* flag to determine if we have refined the mesh and need to rebuild the column structures */
    
    //******************************************************************************//
    //******************************************************************************//
    
    /* loop through the columns and see if we have to refine a column */
    for (i = 0; i < mod->grid->ncolumns; i++) {
        refine_col = 0;
        ptr = mod->grid->column_list[i];
        while (ptr->next != NULL) {
            ie = ptr->id;
            if (mod->grid->elem3d[ie].mat != UNSET_INT && mod->grid->elem_error[ie] > REF_TOL && elem3d_level(mod->grid->elem3d[ie]) < mod->mat[mod->grid->elem3d[ie].mat].sw->max_lev && mod->grid->elem3d[ie].interface==0) {
                refine_col = 1;
                break;
            }
            ptr = ptr->next;
        }
        if (refine_col > 0) {
            column_elem3d_split(mod, i, mod->grid->column_list, edge_hashtab);
            update_cols = 1;
        }
        
    }
    
    //******************************************************************************//
    //******************************************************************************//
    
    /* increment the adaption counter */
    cycle++;
    
    /* Make the grid conforming :: loop until the mesh is conforming */
    do {
        if (update_cols > 0) {
            
            build_column_hash(&mod->grid);
            build_column_list(&mod->grid);
            
            
            update_cols = 0;
            elem2d_renumber(mod->grid); /* need to renumber 2d elements so that the refinement of the 2d elements works */
            classify_2d_elements(&mod->grid); // re-classify elements based on bflag (must occur after elem2d_renumber)
            
        }
        
#ifdef _MESSG
        /* notify other pe's of edges that belong to them that I have split, and they notify me likewise */
        comm_flag_edges(edge_hashtab, mod->grid, nrecv_edge, nsend_edge, mod);
        
        /* update new nodes for the edges in the hash table that I own */
        comm_update_edges(edge_hashtab, mod->grid, nrecv_edge, nsend_edge, mod);
        
        /* correct the adjacencies temporarily set in adpt_get_node */
        for (i = 0; i < mod->grid->nnodes; i++) {
            if ((mod->grid->node[i].parent_res_pe[0] != UNSET_INT) &&
                (mod->grid->node[i].parent_res_pe[0] == mod->grid->smpi->myid) &&
                (mod->grid->node[mod->grid->node[i].parent_res_id[0]].resident_pe != mod->grid->smpi->myid)) {
                mod->grid->node[i].parent_res_pe[0] = mod->grid->node[mod->grid->node[i].parent_res_id[0]].resident_pe;
                mod->grid->node[i].parent_res_id[0] = mod->grid->node[mod->grid->node[i].parent_res_id[0]].resident_id;
                
            }
            if ((mod->grid->node[i].parent_res_pe[1] != UNSET_INT) &&
                (mod->grid->node[i].parent_res_pe[1] == mod->grid->smpi->myid) &&
                (mod->grid->node[mod->grid->node[i].parent_res_id[1]].resident_pe != mod->grid->smpi->myid)) {
                mod->grid->node[i].parent_res_pe[1] = mod->grid->node[mod->grid->node[i].parent_res_id[1]].resident_pe;
                mod->grid->node[i].parent_res_id[1] = mod->grid->node[mod->grid->node[i].parent_res_id[1]].resident_id;
                
            }
        }
#endif
        
        /* initialize the mesh flag to conforming */
        iconform_flag = YES;
        
        /* For columns, see if we have a non conforming element */
        for (i = 0; i < mod->grid->ncolumns; i++) {
            /* Just take first element in column list;  they all need to cut an edge, or none of them need to be refined */
            ptr = mod->grid->column_list[i];
            ie = ptr->id;
            for (iedge = 0; iedge < mod->grid->elem3d[ie].nedges; iedge++) {
                nd1 = mod->grid->elem3d[ie].edges[iedge][0]; nd2 = mod->grid->elem3d[ie].edges[iedge][1];
                e0 = edge_hash_lookup(mod->grid->elem3d[ie].nodes[nd1], mod->grid->elem3d[ie].nodes[nd2], edge_hashtab);
                if (e0 != NULL && e0->new_node != UNSET_INT) {
                    column_elem3d_split(mod, i, mod->grid->column_list, edge_hashtab);
                    update_cols = 1;
                    iconform_flag = NO;
                    break;
                }
            }
        }
        
        /* synchronize the conform flag */
#ifdef _MESSG
        iconform_flag = messg_imin(iconform_flag, mod->grid->smpi->ADH_COMM);
#endif
        
        /* increment the adaption counter */
        cycle++;
    }
    while (iconform_flag == NO);
    
    /* clean up unneeded nodes and elements */
    
#ifdef _MESSG
    partition_cleanup(mod->grid, mod->grid->smpi->myid);
#endif
    
    
    /* reallocate and average physics arrays */
    //adpt_physics_update(mod, old_nnodes, old_nnodes2d, old_max_nnodes, old_max_nnodes2d);
    
    /* renumbers the mesh */
    node_renumber(mod,0);
    elem3d_renumber(mod->grid);
    elem2d_renumber(mod->grid);
    elem1d_renumber(mod->grid);
    
    /* re-build columns */
    build_columns(mod->grid, NO);
    node_renumber_surface(mod);
    /* reset maxes for adaption */
    //mod->grid->old_global_bed = tl_realloc(sizeof(int), mod->grid->nnodes_bed, mod->grid->max_nnodes_bed, mod->grid->old_global_bed);
    //mod->grid->old_global_surf = tl_realloc(sizeof(int), mod->grid->nnodes_sur, mod->grid->max_nnodes_sur, mod->grid->old_global_surf);
    
    // update displacements using nodal coordinates instead of avg dpl (which shouldn't matter since linear)
    tl_vertical_adapt(mod->grid, mod->sw->d3->displacement);
    tl_vertical_adapt(mod->grid, mod->sw->d3->old_displacement);
    tl_vertical_adapt(mod->grid, mod->sw->d3->older_displacement);
    
#ifdef _MESSG
    /* updates the message keys for the repartitioned mesh */
    comm_set_keys(mod->grid);
    
    /* update ghost nodes and node adjacencies */
    comm_update_GN(1, mod->grid);
#endif
    adpt_fix_global_ids(mod->grid);
#ifdef _MESSG
    comm_update_snode(mod->grid);
    for (isd = 0; isd < npes; isd++){
        if (nsend_edge[isd] != 0) {
            mod->grid->smpi->send_edge_key[isd] = (EDGE_LIST_ITEM **) tl_free(sizeof(EDGE_LIST_ITEM *), nsend_edge[isd], mod->grid->smpi->send_edge_key[isd]);
        }
        if (nrecv_edge[isd] != 0) {
            mod->grid->smpi->recv_edge_key[isd] = (EDGE_LIST_ITEM **) tl_free(sizeof(EDGE_LIST_ITEM *), nrecv_edge[isd], mod->grid->smpi->recv_edge_key[isd]);
        }
    }
    nrecv_edge = (int *) tl_free(sizeof(int), npes, nrecv_edge);
    nsend_edge = (int *) tl_free(sizeof(int), npes, nsend_edge);
#endif
    /*update Wind and Wave Series */
    if (mod->flag.WIND) {
        if (mod->flag.WIND_STATION) {
            SSERIES *series = mod->series_wind_head;
            while(series != NULL) {
                
                smeteor_station_realloc(&(series->station), mod->grid->nnodes_sur, series->nnodes);
                series->nnodes = mod->grid->nnodes_sur;
                series = series->next;
            }
            mod->series_wind_head->nnodes = mod->grid->nnodes_sur;
            sseries_set_meteor_stations(mod->series_wind_head, mod->grid, WIND_SERIES);
        }
    }
    if (mod->flag.WAVE) {
        if (mod->flag.WAVE_STATION) {
            SSERIES *series = mod->series_wave_head;
            while(series != NULL) {
                
                smeteor_station_realloc(&(series->station), mod->grid->nnodes_sur, series->nnodes);
                series->nnodes = mod->grid->nnodes_sur;
                series = series->next;
            }
            if(mod->series_wave_head != NULL) mod->series_wave_head->nnodes = mod->grid->nnodes_sur;
            
            sseries_set_meteor_stations(mod->series_wave_head, mod->grid, WAVE_SERIES);
        }
    } 
    /* free the hash table */
    tl_list_free_all(EDGE_LIST);
    
    
    /* return the adaption counter */
    return (cycle);
}

