/* ADH Version 2.0.0 6-04 */
/* this routine splits a 3d elem */

// CJT :: unlikely to work for mixed grids .... 

#include "global_header.h"

#define DIFF_TOL 0.0000001

int split_column_elem(SMODEL *, int, int *, int, int, EDGE_LIST_ITEM **, int);

void column_elem3d_split(
                         SMODEL *mod, int icol,  /* the column we are currently in */
                         ID_LIST_ITEM ** column_list,   /* the linked list for the columns */
                         EDGE_LIST_ITEM ** edge_hashtab /* the hash table of edges */
)
{
    int lo_nd;                    /* the low node on the split edge */
    int hi_nd;                    /* the high node on the split edge */
    int new_node;                 /* the new node */
    int new_node_level;           /* level for the new node */
    int my_new_level;             /* the level for the new nodes should all be identical */
    /* We are going down through the verical prism and splitting any edge that
     * lies below one of the surface edges. Some of the elements will then have
     * two edges that need to be cut */
    int num_splits;               /* counter: 1 or 2 depending on whether there are
                                   1 or 2 edges for this element to be split */
    int split_edge[2];            /* the edges tc be split */
    int new_elem;                 /* the new 3D element */
    int new_elem2d;               /* new 2D element */
    int nd1, nd2;                 /* the nodes on an edge */
    double z_tmp;
    double z_rank;
    EDGE_LIST_ITEM *edge_pntr;    /* pointer to the edge */
    int i, j;
    int istr;
    int isurf;
    int ie;
    int count;
    double length, max_length;
    int is_surf_node[NDPRELM];
    SVECT2D center, midpoint;
    SVECT midpt;
    ID_LIST_ITEM *ptr;
    int ie2;
    int ival;
    
    ie = mod->grid->elem2d_sur[icol];
    
    int **edges = mod->grid->elem2d[ie].edges;
    
    /* selects the edge to be split */
    if (tl_long_edge(mod->grid->node, mod->grid->elem2d[ie].nodes[edges[0][0]], mod->grid->elem2d[ie].nodes[edges[0][1]], mod->grid->elem2d[ie].nodes[edges[1][0]], mod->grid->elem2d[ie].nodes[edges[1][1]], edge_hashtab) == 1)
        split_edge[0] = 0;
    else
        split_edge[0] = 1;
    
    if (tl_long_edge(mod->grid->node, mod->grid->elem2d[ie].nodes[edges[split_edge[0]][0]], mod->grid->elem2d[ie].nodes[edges[split_edge[0]][1]], mod->grid->elem2d[ie].nodes[edges[2][0]], mod->grid->elem2d[ie].nodes[edges[2][1]], edge_hashtab) == 2)
        split_edge[0] = 2;
    
    if (split_edge[0] < 0) {
        printf("\n split_edge[0]: %d\n",split_edge[0]);
        tl_error("ERROR: Unable to split element.\n");
    }
    
    /* sets the hi and lo nodes on the edge */
    lo_nd = edges[split_edge[0]][0];
    hi_nd = edges[split_edge[0]][1];
    
    /* refines the edge if needed */
    nd1 = mod->grid->elem2d[ie].nodes[hi_nd];
    nd2 = mod->grid->elem2d[ie].nodes[lo_nd];
    
    center.x = (mod->grid->node[nd1].x + mod->grid->node[nd2].x) / 2.0;
    center.y = (mod->grid->node[nd1].y + mod->grid->node[nd2].y) / 2.0;
    
    my_new_level = elem2d_level(mod->grid->elem2d[ie]) + 1;
    
    
    /* use the center above to locate the edge to cut for each element in the column */
    ptr = mod->grid->column_list[icol];
    while (ptr->next != NULL) {
        split_edge[0] = -1;
        split_edge[1] = -1;
        num_splits = 0;
        
        ie = ptr->id;
        edges = mod->grid->elem3d[ie].edges;
        for (i = 0; i < mod->grid->elem3d[ie].nedges/* NEDGEPRELM */; i++) {
            /* find the midpoint of each edge */
            nd1 = mod->grid->elem3d[ie].nodes[edges[i][0]];
            nd2 = mod->grid->elem3d[ie].nodes[edges[i][1]];
            midpoint.x = (mod->grid->node[nd1].x + mod->grid->node[nd2].x) / 2.0;
            midpoint.y = (mod->grid->node[nd1].y + mod->grid->node[nd2].y) / 2.0;
            z_tmp = (mod->grid->node[nd1].z + mod->grid->node[nd2].z) / 2.0;
            
            /* check to see if the edge midpoint lies in vertical line with the center from above */
            if ((fabs(center.x - midpoint.x) < DIFF_TOL) && (fabs(center.y - midpoint.y) < DIFF_TOL)) {
                /* If there are two edges to be split, cut the bottom one first */
                if (num_splits == 0) {
                    /* found the first edge */
                    split_edge[0] = i;
                    z_rank = z_tmp;
                }
                else {
                    /* found a second edge */
                    if (z_tmp < z_rank) {
                        /* This edge lies below the previous edge */
                        split_edge[1] = split_edge[0];
                        split_edge[0] = i;
                    }
                    else {
                        /* This edge lies above the previously found edge */
                        split_edge[1] = i;
                    }
                }
                num_splits++;
            }
        }
        ptr = ptr->next;
        
        if (split_edge[0] < 0) {
            tl_error("ERROR: Unable to split element.\n");
        }
        
        new_elem = split_column_elem(mod, ie, &(split_edge[0]), my_new_level, 0, edge_hashtab, NO);
        
        if (mod->grid->elem3d[ie].nnodes == NDONTET){
            /* If num_splits > 1, then we need to further refine the two elements */
            if (num_splits > 1) {
                nd1 = mod->grid->elem3d[ie].nodes[edges[split_edge[1]][0]];
                nd2 = mod->grid->elem3d[ie].nodes[edges[split_edge[1]][1]];
                if (edge_hash_lookup(nd1, nd2, edge_hashtab) != NULL) {
                    split_column_elem(mod, ie, &(split_edge[1]), UNSET_INT, 0, edge_hashtab, YES); // cjt :: force split
                }
                nd1 = mod->grid->elem3d[new_elem].nodes[edges[split_edge[1]][0]];
                nd2 = mod->grid->elem3d[new_elem].nodes[edges[split_edge[1]][1]];
                if (edge_hash_lookup(nd1, nd2, edge_hashtab) != NULL) {
                    split_column_elem(mod, new_elem, &(split_edge[1]), UNSET_INT, 0, edge_hashtab, YES); // cjt :: force split
                }
            }
        }
    }
    
    /* Loop over the 2D elements in this column and split them */
    ptr = mod->grid->column_list2d[icol];
    while (ptr->next != NULL) {
        split_edge[0] = -1;
        split_edge[1] = -1;
        num_splits = 0;
        
        ie = ptr->id;
        edges = mod->grid->elem2d[ie].edges; /* Gajanan gkc bug fix */
        if (mod->grid->elem2d[ie].nedges == NDONQUAD) {
            tl_error("\n    Found a quadrilateral element where a triangle was expected.");
        }
        for (i = 0; i < mod->grid->elem2d[ie].nedges /* NEDGEPRFC */; i++) {
            /* find the midpoint of each edge */
            nd1 = mod->grid->elem2d[ie].nodes[edges[i][0]];
            nd2 = mod->grid->elem2d[ie].nodes[edges[i][1]];
            midpoint.x = (mod->grid->node[nd1].x + mod->grid->node[nd2].x) / 2.0;
            midpoint.y = (mod->grid->node[nd1].y + mod->grid->node[nd2].y) / 2.0;
            z_tmp = (mod->grid->node[nd1].z + mod->grid->node[nd2].z) / 2.0;
            
            /* check to see if the edge midpoint lies in vertical line
             with the center from above */
            if ((fabs(center.x - midpoint.x) < DIFF_TOL) && (fabs(center.y - midpoint.y) < DIFF_TOL)) {
                /* If there are two edges to be split, put the one on top first */
                if (num_splits == 0) {
                    split_edge[0] = i;
                    z_rank = z_tmp;
                }
                else {
                    if (z_tmp > z_rank) {
                        /* This edge lies on top of previous edge */
                        split_edge[1] = split_edge[0];
                        split_edge[0] = i;
                    }
                    else {
                        /* This edge lies below the previously found edge */
                        split_edge[1] = i;
                    }
                }
                num_splits++;
            }
        }
        ptr = ptr->next;
        
        /* Found an edge to split. This should be true, unless we have an element on a side wall, and that wall is not being split. */
        if (split_edge[0] >= 0) {
            new_elem2d = column_elem2d_split(mod, ie, &(split_edge[0]), edge_hashtab);
        }
        
        if (mod->grid->elem2d[ie].nnodes == NDONTRI){
            /* If num_splits > 1, then we need to further refine the two elements */
            if (num_splits > 1) {
                nd1 = mod->grid->elem2d[ie].nodes[edges[split_edge[1]][0]];
                nd2 = mod->grid->elem2d[ie].nodes[edges[split_edge[1]][1]];
                if (edge_hash_lookup(nd1, nd2, edge_hashtab) != NULL) {
                    /* Gajanan gkc bugfix - looks like there might have been a bug here - replacing split_edge[0] with splid_edge[1]. */
                    column_elem2d_split(mod, ie, &(split_edge[1]), edge_hashtab);
                    //column_elem2d_split(mod, ie, split_edge[0], edge_hashtab);
                }
                else {
                    nd1 = mod->grid->elem2d[new_elem2d].nodes[edges[split_edge[1]][0]];
                    nd2 = mod->grid->elem2d[new_elem2d].nodes[edges[split_edge[1]][1]];
                    if (edge_hash_lookup(nd1, nd2, edge_hashtab) != NULL) {
                        /* Gajanan gkc bugfix - looks like there might have been a bug here - replacing split_edge[0] with splid_edge[1]. */
                        column_elem2d_split(mod, new_elem2d, &(split_edge[1]), edge_hashtab);
                        //column_elem2d_split(mod, new_elem2d, split_edge[0], edge_hashtab);
                    }
                }
            }
        }
    }
    
    /* Now check to see if we need to split any sidewall elements */
    midpt.x = center.x;
    midpt.y = center.y;
    midpt.z = 0.0;
    
    /* Find the midpt segment corresponding to the new node */
    ival = find_midpt_segment(mod->grid, midpt, mod->grid->midpt_hash);
    ptr = mod->grid->sidewall_list[ival];
    while (ptr->next != NULL) {
        split_edge[0] = -1;
        split_edge[1] = -1;
        num_splits = 0;
        
        ie = ptr->id;
        edges = mod->grid->elem2d[ie].edges; /* Gajanan gkc bug fix */
        for (i = 0; i < mod->grid->elem2d[ie].nedges /* NEDGEPRFC */; i++) {
            /* find the midpoint of each edge */
            nd1 = mod->grid->elem2d[ie].nodes[edges[i][0]];
            nd2 = mod->grid->elem2d[ie].nodes[edges[i][1]];
            midpoint.x = (mod->grid->node[nd1].x + mod->grid->node[nd2].x) / 2.0;
            midpoint.y = (mod->grid->node[nd1].y + mod->grid->node[nd2].y) / 2.0;
            z_tmp = (mod->grid->node[nd1].z + mod->grid->node[nd2].z) / 2.0;
            
            /* check to see if the edge midpoint lies in vertical line with the center from above */
            if ((fabs(center.x - midpoint.x) < DIFF_TOL) && (fabs(center.y - midpoint.y) < DIFF_TOL)) {
                /* If there are two edges to be split, put the one on top first */
                if (num_splits == 0) {
                    split_edge[0] = i;
                    z_rank = z_tmp;
                }
                else {
                    if (z_tmp < z_rank) {
                        /* This edge lies on top of previous edge */
                        split_edge[1] = split_edge[0];
                        split_edge[0] = i;
                    }
                    else {
                        /* This edge lies below the previously found edge */
                        split_edge[1] = i;
                    }
                }
                num_splits++;
            }
        }
        ptr = ptr->next;
        
        /* Found an edge to split. This should be true, unless we have an element on a side
         * wall, and that wall is not being split.
         */
        if (split_edge[0] >= 0) {
            new_elem2d = column_elem2d_split(mod, ie, &(split_edge[0]), edge_hashtab);
        }
        
        if (mod->grid->elem2d[ie].nnodes == NDONTRI){
            /* If num_splits > 1, then we need to further refine the two elements */
            if (num_splits > 1) {
                nd1 = mod->grid->elem2d[ie].nodes[edges[split_edge[1]][0]];
                nd2 = mod->grid->elem2d[ie].nodes[edges[split_edge[1]][1]];
                if (edge_hash_lookup(nd1, nd2, edge_hashtab) != NULL) {
                    column_elem2d_split(mod, ie, &(split_edge[1]), edge_hashtab);
                }
                nd1 = mod->grid->elem2d[new_elem2d].nodes[edges[split_edge[1]][0]];
                nd2 = mod->grid->elem2d[new_elem2d].nodes[edges[split_edge[1]][1]];
                if (edge_hash_lookup(nd1, nd2, edge_hashtab) != NULL) {
                    column_elem2d_split(mod, new_elem2d, &(split_edge[1]), edge_hashtab);
                }
            }
        }
    }
    
}

int split_column_elem(SMODEL *mod,
                      int ie,    /* the element to be split */
                      int *split_edge,  /* the edge(s) to be split */
                      int node_level,   /* the node_level to be assigned to the new node */
                      int force_level,  /* indicates if this is a secondary refinement in column */
                      EDGE_LIST_ITEM ** edge_hashtab,    /* the hash table of edges */
                      int forced  /* cjt :: force splitting */
)
{
    int i;
    int lo_nd1, lo_nd2;    /* the low node on the split edge */
    int hi_nd1, hi_nd2;    /* the high node on the split edge */
    int new_node1, new_node2;                 /* the new node */
    int nd1, nd2;                 /* the nodes on an edge */
    EDGE_LIST_ITEM *edge_pntr;    /* pointer to the edge */
    int new_elem;                 /* the new element */
    SNODE nodes[NDONPRISM];
    
    if (mod->grid->elem3d[ie].nnodes == NDONTET) {
        /* sets the hi and lo nodes on the edge */
        hi_nd1 = mod->grid->elem3d[ie].edges[split_edge[0]][0];
        lo_nd1 = mod->grid->elem3d[ie].edges[split_edge[0]][1];
        
        /* refines the edge if needed */
        nd1 = mod->grid->elem3d[ie].nodes[hi_nd1];
        nd2 = mod->grid->elem3d[ie].nodes[lo_nd1];
        edge_pntr = edge_hash_lookup(nd1, nd2, edge_hashtab);
        new_node1 = adpt_get_node(mod, edge_pntr, edge_hashtab);
        
        /* gets the new element */
        new_elem = elem3d_new(mod->grid, mod->nalloc_inc, mod->grid->elem3d[ie].nnodes);
        
        /* copies the current element to the new element */
        selem3d_copy(&(mod->grid->elem3d[new_elem]), mod->grid->elem3d[ie]);
        
        /* corrects the current elements data */
        mod->grid->elem3d[ie].nodes[lo_nd1] = new_node1;
        mod->grid->elem_error[ie] = 0.0;
        for (i=0; i<NDONTET; i++) {
            snode_copy(&(nodes[i]), mod->grid->node[mod->grid->elem3d[ie].nodes[i]]);
        }
        get_tet_linear_djac_gradPhi(&(mod->grid->elem3d[ie]), nodes, NULL);
        
        /* corrects the new elements data */
        mod->grid->elem3d[new_elem].nodes[hi_nd1] = new_node1;
        mod->grid->elem_error[new_elem] = 0.0;
        for (i=0; i<NDONTET; i++) {
            snode_copy(&(nodes[i]), mod->grid->node[mod->grid->elem3d[new_elem].nodes[i]]);
        }
        get_tet_linear_djac_gradPhi(&(mod->grid->elem3d[new_elem]), nodes, NULL);
        
        int new_node_level = UNSET_INT;
        if (forced == YES) {  // leave the level alone
            new_node_level = elem3d_level(mod->grid->elem3d[ie]);
        } else {
            new_node_level = node_level;
        }
        
        mod->grid->elem3d[ie].levels[lo_nd1] = new_node_level;
        mod->grid->elem3d[new_elem].levels[hi_nd1] = new_node_level;
    }
    else if (mod->grid->elem3d[new_elem].nnodes == NDONPRISM) {
        /* sets the hi and lo nodes on the edge */
        hi_nd1 = mod->grid->elem3d[ie].edges[split_edge[0]][0];
        lo_nd1 = mod->grid->elem3d[ie].edges[split_edge[0]][1];
        hi_nd2 = mod->grid->elem3d[ie].edges[split_edge[1]][0];
        lo_nd2 = mod->grid->elem3d[ie].edges[split_edge[1]][1];
        
        /* refines the edge if needed */
        nd1 = mod->grid->elem3d[ie].nodes[hi_nd1];
        nd2 = mod->grid->elem3d[ie].nodes[lo_nd1];
        edge_pntr = edge_hash_lookup(nd1, nd2, edge_hashtab);
        new_node1 = adpt_get_node(mod, edge_pntr, edge_hashtab);
        nd1 = mod->grid->elem3d[ie].nodes[hi_nd2];
        nd2 = mod->grid->elem3d[ie].nodes[lo_nd2];
        edge_pntr = edge_hash_lookup(nd1, nd2, edge_hashtab);
        new_node2 = adpt_get_node(mod, edge_pntr, edge_hashtab);
        
        /* gets the new element */
        new_elem = elem3d_new(mod->grid, mod->nalloc_inc, mod->grid->elem3d[ie].nnodes);
        
        /* copies the current element to the new element */
        selem3d_copy(&(mod->grid->elem3d[new_elem]), mod->grid->elem3d[ie]);
        
        /* corrects the current elements data */
        mod->grid->elem3d[ie].nodes[lo_nd1] = new_node1;
        mod->grid->elem3d[ie].nodes[lo_nd2] = new_node2;
        mod->grid->elem_error[ie] = 0.0;
        for (i=0; i<NDONPRISM; i++) {
            snode_copy(&(nodes[i]), mod->grid->node[mod->grid->elem3d[ie].nodes[i]]);
        }
        /* Gajanan gkc - I don't know what to do here. */
        //get_triprism_linear_djac_gradPhi(&(mod->grid->elem3d[ie]), nodes, NULL);
        
        /* corrects the new elements data */
        mod->grid->elem3d[new_elem].nodes[hi_nd1] = new_node1;
        mod->grid->elem3d[new_elem].nodes[hi_nd2] = new_node2;
        mod->grid->elem_error[new_elem] = 0.0;
        for (i=0; i<NDONPRISM; i++) {
            snode_copy(&(nodes[i]), mod->grid->node[mod->grid->elem3d[new_elem].nodes[i]]);
        }
        /* Gajanan gkc - I don't know what to do here. */
        //get_triprism_linear_djac_gradPhi(&(mod->grid->elem3d[new_elem]), nodes, NULL);
        
        /* This is likely not going to happen for any prism. */
        int new_node_level = UNSET_INT;
        if (forced == YES) {  // leave the level alone
            new_node_level = elem3d_level(mod->grid->elem3d[ie]);
        } else {
            new_node_level = node_level;
        }
        
        mod->grid->elem3d[ie].levels[lo_nd1] = new_node_level;
        mod->grid->elem3d[ie].levels[lo_nd2] = new_node_level;
        mod->grid->elem3d[new_elem].levels[hi_nd1] = new_node_level;
        mod->grid->elem3d[new_elem].levels[hi_nd2] = new_node_level;
    }
    
    return (new_elem);
    
}
