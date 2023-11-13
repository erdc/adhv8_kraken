/* refines the mesh */
#include "global_header.h"

long adpt_ref(SMODEL *mod, long cycle) {
    
    int ie;      /* loop counter over the elements */
    int i, isd;      /* loop counters */
    int iconform_flag;    /* flag for whether or not a mesh is conforming */
    int nd1, nd2, iedge, nedges;
    int **edges;
#ifdef _MESSG
    int *nrecv_edge, *nsend_edge;
    int npes = mod->grid->smpi->npes;
#endif
    
    SGRID *grid = mod->grid;
    assert(grid->type != COLUMNAR);
    
    //print_parents_to_file(grid,"PARENTS_PREREF");
    EDGE_LIST_ITEM *edge_hashtab[HASHSIZE];  /* the hash table of edges - this serves two purposes:
                                              the first is to hold the integer rank of an edge for
                                              length comparisons, and the second is to hold the number
                                              of a node imbedded in the edge.  This allows us to
                                              determine if an edge has been split */
    
    EDGE_LIST_ITEM *e0, *e1, *e2, *e3, *e4, *e5;  /* the edges of an element */
    
    /* initialize the hash table */
    for (i = 0; i < HASHSIZE; i++) {
        edge_hashtab[i] = NULL;
    }
    
    /* load the edges into the hash table */
    for (ie = 0; ie < grid->nelems3d; ie++) {
        for (iedge = 0; iedge < grid->elem3d[ie].nedges; iedge++) {
            nd1 = grid->elem3d[ie].edges[iedge][0]; nd2 = grid->elem3d[ie].edges[iedge][1];
            edge_hash_add_entry(grid->elem3d[ie].nodes[nd1], grid->elem3d[ie].nodes[nd2], edge_hashtab);
            edge_hash_add_entry(grid->elem3d[ie].nodes[nd2], grid->elem3d[ie].nodes[nd1], edge_hashtab); /* gkc adding to remove edge not found error in adpt_rank_edges */
        }
    }
    
    for (ie = 0; ie < grid->nelems2d; ie++) {
        for (iedge = 0; iedge < grid->elem2d[ie].nedges; iedge++) {
            nd1 = grid->elem2d[ie].edges[iedge][0]; nd2 = grid->elem2d[ie].edges[iedge][1];
            edge_hash_add_entry(grid->elem2d[ie].nodes[nd1], grid->elem2d[ie].nodes[nd2], edge_hashtab);
            edge_hash_add_entry(grid->elem2d[ie].nodes[nd2], grid->elem2d[ie].nodes[nd1], edge_hashtab); /* gkc adding to remove edge not found error in adpt_rank_edges */
        }
    }
    
    for (ie = 0; ie < grid->nelems1d; ie++) {
        edge_hash_add_entry(grid->elem1d[ie].nodes[0], grid->elem1d[ie].nodes[1], edge_hashtab);
        edge_hash_add_entry(grid->elem1d[ie].nodes[1], grid->elem1d[ie].nodes[0], edge_hashtab);
    }
    
    /* set the edge keys */
#ifdef _MESSG
    nrecv_edge = (int *) tl_alloc(sizeof(int), npes);
    nsend_edge = (int *) tl_alloc(sizeof(int), npes);
    comm_edge_setup(edge_hashtab, grid, nrecv_edge, nsend_edge);
#endif
    
    /* rank the edges */
#ifdef _MESSG
    adpt_rank_edges(edge_hashtab, grid, nrecv_edge, nsend_edge);
#else
    adpt_rank_edges(edge_hashtab, grid);
#endif
    
   
    int max_lev = UNSET_INT, k=0;
    
    if (grid->ndim == 3) {  // Add 3D resolution here only if the grid is fully unstructured (note :: also splits 2D if external boundary)
        for (ie = 0; ie < grid->nelems3d; ie++) {
#ifdef _ADH_GROUNDWATER
            if (mod->flag.NS3_FLOW) max_lev = mod->mat[grid->elem3d[ie].mat].ns->max_lev;
            else if (mod->flag.GW_FLOW) max_lev = mod->mat[grid->elem3d[ie].mat].gw->max_lev;
#else
            max_lev = mod->mat[grid->elem3d[ie].mat].ns->max_lev;
#endif
            if (grid->elem3d[ie].mat != UNSET_INT && grid->elem_error[ie] > REF_TOL && elem3d_level(grid->elem3d[ie]) < max_lev) {
                int force_adaption = NO;
                elem3d_split(mod, ie, edge_hashtab, force_adaption);
                k++;
            }
        }
    } else if (grid->ndim == 2) {  // Add 2D resolution here (note :: also splits 1D if external boundary)
        for (ie = 0; ie < grid->nelems2d; ie++) {
            if (grid->elem2d[ie].mat != UNSET_INT && grid->elem_error[ie] > REF_TOL && elem2d_level(grid->elem2d[ie]) < max_lev && grid->elem2d[ie].interface==0) {
                elem2d_split(mod, ie, edge_hashtab);
            }
        }
    } else {
        //for (ie = 0; ie < grid->nelems1d; ie++) {
        //    if (elem1d[ie].string != UNSET_INT && grid->elem_error[ie] > REF_TOL && elem1d_level(elem1d[ie]) < max_lev) {
        //        elem1d_split(mod, ie, edge_hashtab);
        //    }
        //}
    }
    
    
    /* increment the adaption counter */
    cycle++;
    
    //   print_grid_to_file(grid,"POST_SPLIT");
    //   print_parents_to_file(grid,"PARENTS_SPLIT");
    /* Make the grid conforming */
    do {
        
#ifdef _MESSG
        /* notify other pe's of edges that belong to them that I have split,
         and they notify me likewise */
        
        comm_flag_edges(edge_hashtab, grid, nrecv_edge, nsend_edge, mod);
        
        /* update new nodes for the edges in the hash table that I own */
        comm_update_edges(edge_hashtab, grid, nrecv_edge, nsend_edge, mod);
        
        /* correct the adjacencies temporarily set in adpt_get_node */
        
        for (i = 0; i < grid->nnodes; i++) {
            if ((grid->node[i].parent_res_pe[0] != UNSET_INT) &&
                (grid->node[i].parent_res_pe[0] == grid->smpi->myid) &&
                (grid->node[grid->node[i].parent_res_id[0]].resident_pe != grid->smpi->myid)) {
                grid->node[i].parent_res_pe[0] = grid->node[grid->node[i].parent_res_id[0]].resident_pe;
                grid->node[i].parent_res_id[0] = grid->node[grid->node[i].parent_res_id[0]].resident_id;
                
            }
            if ((grid->node[i].parent_res_pe[1] != UNSET_INT) &&
                (grid->node[i].parent_res_pe[1] == grid->smpi->myid) &&
                (grid->node[grid->node[i].parent_res_id[1]].resident_pe != grid->smpi->myid)) {
                grid->node[i].parent_res_pe[1] = grid->node[grid->node[i].parent_res_id[1]].resident_pe;
                grid->node[i].parent_res_id[1] = grid->node[grid->node[i].parent_res_id[1]].resident_id;
                
            }
        }
#endif
        
        /* initialize the mesh flag to conforming */
        iconform_flag = YES;
        
        /* loop over the 3d elements and split the non conforming elements */
        for (ie = 0; ie < grid->nelems3d; ie++) {
            if (grid->elem3d[ie].mat != UNSET_INT) {
                for (iedge = 0; iedge < grid->elem3d[ie].nedges; iedge++) {
                    nd1 = grid->elem3d[ie].edges[iedge][0]; nd2 = grid->elem3d[ie].edges[iedge][1];
                    e0 = edge_hash_lookup(grid->elem3d[ie].nodes[nd1], grid->elem3d[ie].nodes[nd2], edge_hashtab);
                    if (e0 != NULL && e0->new_node != UNSET_INT) {
                        elem3d_split(mod, ie, edge_hashtab, NO);
                        iconform_flag = NO;
                        break;
                    }
                }
            }
        }
        
        //// all new nodes == UNSET_INT!!!!  NOW THAT I LOOK, IT MAY NOT NEED TO SPLIT.  IT looks like it works fine for level 2
        
        /* loop over the 2d elements and split the non conforming elements */
        for (ie = 0; ie < grid->nelems2d; ie++) {
            if (grid->elem2d[ie].mat != UNSET_INT) {
                for (iedge = 0; iedge < grid->elem2d[ie].nedges; iedge++) {
                    nd1 = grid->elem2d[ie].edges[iedge][0]; nd2 = grid->elem2d[ie].edges[iedge][1];
                    e0 = edge_hash_lookup(grid->elem2d[ie].nodes[nd1], grid->elem2d[ie].nodes[nd2], edge_hashtab);
                    if (e0 != NULL && e0->new_node != UNSET_INT) {
                        elem2d_split(mod, ie, edge_hashtab);
                        iconform_flag = NO;
                        break;
                    }
                }
            }
        }
        
        /* loop over the 1d elements and split the non conforming elements */
        for (ie = 0; ie < grid->nelems1d; ie++) {
            if (grid->elem1d[ie].string != UNSET_INT) {
                /* look up the edges */
                e0 = edge_hash_lookup(grid->elem1d[ie].nodes[0], grid->elem1d[ie].nodes[1], edge_hashtab);
                
                /* checks that the edges are conforming */
                if (e0 != NULL && e0->new_node != UNSET_INT) {
                    elem1d_split(mod, ie, edge_hashtab);
                    iconform_flag = NO;
                }
            }
        }
        
        /* synchronize the conform flag */
#ifdef _MESSG
        iconform_flag = messg_imin(iconform_flag, grid->smpi->ADH_COMM);
#endif
        
        /* increment the adaption counter */
        cycle++;
    }
    while (iconform_flag == NO);
    
    //  print_parents_to_file(grid,"PARENTS_CONFORM");
    /* clean up unneeded nodes and elements */
#ifdef _MESSG
    partition_cleanup(grid, grid->smpi->myid);
#endif
    
    
    //   print_parents_to_file(grid,"PARENTS_CLEAN");
    //    print_grid_to_file(grid,"POST_CLEAN");
    /* renumbers the mesh */
    node_renumber(mod, 0);
    elem3d_renumber(grid);
    elem2d_renumber(grid);
    elem1d_renumber(grid);
    if (grid->ndim == 3) classify_2d_elements(&grid);
    
#ifdef _MESSG
    /* updates the message keys for the repartitioned mesh */
    comm_set_keys(grid);
    comm_update_GN(1, grid);
    //comm_update_snode(grid);messg_barrier(MPI_COMM_WORLD);
    adpt_fix_global_ids(grid);
    comm_update_snode(grid);
    
    /* update ghost nodes and node adjacencies */
    //comm_update_GLOBAL_NODE(1, node_pair);
    //comm_update_GLOBAL_NODE(0, node_ladj);
    //comm_update_GLOBAL_NODE(0, node_radj);
    for (isd = 0; isd < npes; isd++){
        if (nsend_edge[isd] != 0) {
            grid->smpi->send_edge_key[isd] = (EDGE_LIST_ITEM **) tl_free(sizeof(EDGE_LIST_ITEM *), nsend_edge[isd], grid->smpi->send_edge_key[isd]);
        }
        if (nrecv_edge[isd] != 0) {
            grid->smpi->recv_edge_key[isd] = (EDGE_LIST_ITEM **) tl_free(sizeof(EDGE_LIST_ITEM *), nrecv_edge[isd], grid->smpi->recv_edge_key[isd]);
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
                
                smeteor_station_realloc(&(series->station), mod->grid->nnodes, series->nnodes);
                series->nnodes = grid->nnodes_sur;
                series = series->next;
            }
            mod->series_wind_head->nnodes = mod->grid->nnodes;
            sseries_set_meteor_stations(mod->series_wind_head, mod->grid, WIND_SERIES);
        }
    }
    if (mod->flag.WAVE) {
        if (mod->flag.WAVE_STATION) {
            SSERIES *series = mod->series_wave_head;
            while(series != NULL) {
                
                smeteor_station_realloc(&(series->station), mod->grid->nnodes, series->nnodes);
                series->nnodes = mod->grid->nnodes;
                series = series->next;
            }
            if(mod->series_wave_head != NULL) mod->series_wave_head->nnodes = mod->grid->nnodes;
            
            sseries_set_meteor_stations(mod->series_wave_head, mod->grid, WAVE_SERIES);
        }
    }
    /* free the hash table */
    tl_list_free_all(EDGE_LIST);
    //print_grid_to_file(grid,"POST_REF");
    /* return the adaption counter */
    
    
    //printf("\nnelem3d: %d \t k: %d\n",grid->nelems3d,grid->nelems3d_old);
    //printf("nelem2d: %d \t k: %d\n",grid->nelems2d,grid->nelems2d_old);
    //for (i=0; i<grid->nnodes; i++) {
    //    if (grid->node[i].original_id == UNSET_INT && (grid->node[i].z == 0 || fabs(grid->node[i].z + 10) < 1e-6)) printf("z: %20.10f\n",grid->node[i].z);
    //    if (grid->node[i].original_id == UNSET_INT) printf("z: %20.10f\n",grid->node[i].z);
    //}
    //exit(-1);
    
    return (cycle);
}
