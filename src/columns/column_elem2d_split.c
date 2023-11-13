/* ADH Version 2.0.0 6-04 */
/* this routine splits a 2d elem  which is being split after during a call to elem_3d_split */

#include "global_header.h"

int column_elem2d_split(SMODEL *mod,
                        int ie,  /* the element to be split */
                        int *split_edge, /* the edge to split */
                        EDGE_LIST_ITEM ** edge_hashtab  /* the hash table of edges */
)
{
    int i;
    int lo_nd1, lo_nd2;           /* the low node on the split edge */
    int hi_nd1, hi_nd2;           /* the high node on the split edge */
    int new_node1, new_node2;     /* the new node */
    int new_node_level;           /* level for the new node */
    int new_elem;                 /* the new element */
    int nd1, nd2;                 /* the nodes on an edge */
    EDGE_LIST_ITEM *edge_pntr;    /* pointer to the edge */
    
    int **edges;
    if (mod->grid->elem2d[ie].nnodes == NDONTRI) {
        edges = mod->grid->nd_on_TriEdge;

        /* sets the hi and lo nodes on the edge */
        lo_nd1 = edges[split_edge[0]][0];
        hi_nd1 = edges[split_edge[0]][1];
        
        /* refines the edge if needed */
        nd1 = mod->grid->elem2d[ie].nodes[hi_nd1];
        nd2 = mod->grid->elem2d[ie].nodes[lo_nd1];
        edge_pntr = edge_hash_lookup(nd1, nd2, edge_hashtab);
        new_node1 = adpt_get_node(mod, edge_pntr, edge_hashtab);
        
        /* gets the new element */
        new_elem = elem2d_new(mod->grid, mod->nalloc_inc, mod->grid->elem2d[ie].nnodes);
        
        /* copies the current element to the new element */
        selem2d_copy(&(mod->grid->elem2d[new_elem]), mod->grid->elem2d[ie]);
        
        /* calculate level for the new node */
        new_node_level = elem2d_level(mod->grid->elem2d[ie]) + 1;
        
        /* corrects the current elements data */
        mod->grid->elem2d[ie].nodes[lo_nd1] = new_node1;
        mod->grid->elem2d[ie].levels[lo_nd1] = new_node_level;
        mod->grid->elem2d[ie].djac3d_fixed = 0.; // zero this so that lingrad recalculates
        get_triangle_linear_djac_nrml_gradPhi(&(mod->grid->elem2d[ie]), mod->grid->node, NULL);
        
        /* corrects the new elements data */
        mod->grid->elem2d[new_elem].nodes[hi_nd1] = new_node1;
        mod->grid->elem2d[new_elem].levels[hi_nd1] = new_node_level;
        mod->grid->elem2d[new_elem].djac3d_fixed = 0.; // zero this so that lingrad recalculates
        get_triangle_linear_djac_nrml_gradPhi(&(mod->grid->elem2d[new_elem]), mod->grid->node, NULL);
    
    } else {
        edges = mod->grid->nd_on_QuadEdge;
        
        /* sets the hi and lo nodes on the edge */
        lo_nd1 = edges[split_edge[0]][0];
        hi_nd1 = edges[split_edge[0]][1];
        lo_nd2 = edges[split_edge[1]][1]; /* Gajanan gkc - Note the numbers!!! This needs to be looked at more rigorously. */
        hi_nd2 = edges[split_edge[1]][0];
        
        /* refines the edge if needed */
        nd1 = mod->grid->elem2d[ie].nodes[hi_nd1];
        nd2 = mod->grid->elem2d[ie].nodes[lo_nd1];
        edge_pntr = edge_hash_lookup(nd1, nd2, edge_hashtab);
        new_node1 = adpt_get_node(mod, edge_pntr, edge_hashtab);
        
        nd1 = mod->grid->elem2d[ie].nodes[hi_nd2];
        nd2 = mod->grid->elem2d[ie].nodes[lo_nd2];
        edge_pntr = edge_hash_lookup(nd1, nd2, edge_hashtab);
        new_node2 = adpt_get_node(mod, edge_pntr, edge_hashtab);
        
        /* gets the new element */
        new_elem = elem2d_new(mod->grid, mod->nalloc_inc, mod->grid->elem2d[ie].nnodes);
        
        /* copies the current element to the new element */
        selem2d_copy(&(mod->grid->elem2d[new_elem]), mod->grid->elem2d[ie]);
        
        /* calculate level for the new node */
        new_node_level = elem2d_level(mod->grid->elem2d[ie]) + 1;
        
        /* corrects the current elements data */
        mod->grid->elem2d[ie].nodes[lo_nd1] = new_node1;
        mod->grid->elem2d[ie].nodes[lo_nd2] = new_node2;
        mod->grid->elem2d[ie].levels[lo_nd1] = new_node_level;
        mod->grid->elem2d[ie].levels[lo_nd2] = new_node_level;
        mod->grid->elem2d[ie].djac3d_fixed = 0.; // zero this so that lingrad recalculates
        SVECT nds[mod->grid->elem2d[ie].nnodes];
        for (i=0; i<mod->grid->elem2d[ie].nnodes; i++) {
            nds[i].x = mod->grid->node[ mod->grid->elem2d[ie].nodes[i] ].x;
            nds[i].y = mod->grid->node[ mod->grid->elem2d[ie].nodes[i] ].y;
            nds[i].z = mod->grid->node[ mod->grid->elem2d[ie].nodes[i] ].z;  // initial displacement should be added here
        }
        mod->grid->elem2d[ie].nrml = get_elem2d_normals(nds);
        
        /* corrects the new elements data */
        mod->grid->elem2d[new_elem].nodes[hi_nd1] = new_node1;
        mod->grid->elem2d[new_elem].nodes[hi_nd2] = new_node2;
        mod->grid->elem2d[new_elem].levels[hi_nd1] = new_node_level;
        mod->grid->elem2d[new_elem].levels[hi_nd2] = new_node_level;
        mod->grid->elem2d[new_elem].djac3d_fixed = 0.; // zero this so that lingrad recalculates
        for (i=0; i<mod->grid->elem2d[new_elem].nnodes; i++) {
            nds[i].x = mod->grid->node[ mod->grid->elem2d[new_elem].nodes[i] ].x;
            nds[i].y = mod->grid->node[ mod->grid->elem2d[new_elem].nodes[i] ].y;
            nds[i].z = mod->grid->node[ mod->grid->elem2d[new_elem].nodes[i] ].z;  // initial displacement should be added here
        }
        mod->grid->elem2d[new_elem].nrml = get_elem2d_normals(nds);
    
    }

    return new_elem;
}
