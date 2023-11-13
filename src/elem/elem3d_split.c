// CJT :: Splits a 3D element
// notes :: currently only for tetrahedrons although the algorithm is general
// notes :: will probably need some work for prisms, check column_elem3d_split
#include "global_header.h"

void elem3d_split(SMODEL *mod,
                  int ie,                           /* the element to be split */
                  EDGE_LIST_ITEM ** edge_hashtab,   /* the hash table of edges */
                  int forced  /* cjt :: force splitting */
)
{
    int lo_nd;              /* the low node on the split edge */
    int hi_nd;              /* the high node on the split edge */
    int new_node;			/* the new node */
    int split_edge;         /* the split edge */
    int new_elem;			/* the new element */
    int nd1, nd2;			/* the nodes on an edge */
    int i, iedge;
    EDGE_LIST_ITEM *edge_pntr;	/* pointer to the edge */
    SNODE nodes[NDONPRISM];
    
    // currently only for tetrahedrons although the algorithm is general
    assert(mod->grid->elem3d[ie].nnodes == NDONTET);
    
    // selects the edge to be split
    if (tl_long_edge(mod->grid->node,
                     mod->grid->elem3d[ie].nodes[mod->grid->elem3d[ie].edges[0][0]],
                     mod->grid->elem3d[ie].nodes[mod->grid->elem3d[ie].edges[0][1]],
                     mod->grid->elem3d[ie].nodes[mod->grid->elem3d[ie].edges[1][0]],
                     mod->grid->elem3d[ie].nodes[mod->grid->elem3d[ie].edges[1][1]],
                     edge_hashtab) == 1) {
        split_edge = 0;
    } else {split_edge = 1;}
    
    for (iedge=2;iedge<mod->grid->elem3d[ie].nedges; iedge++) {
        if (tl_long_edge(mod->grid->node,
                         mod->grid->elem3d[ie].nodes[mod->grid->elem3d[ie].edges[split_edge][0]],
                         mod->grid->elem3d[ie].nodes[mod->grid->elem3d[ie].edges[split_edge][1]],
                         mod->grid->elem3d[ie].nodes[mod->grid->elem3d[ie].edges[iedge][0]],
                         mod->grid->elem3d[ie].nodes[mod->grid->elem3d[ie].edges[iedge][1]],
                         edge_hashtab) == 2) {
            split_edge = iedge;
        }
    }

    
    // sets the hi and lo nodes on the edge
    lo_nd = mod->grid->elem3d[ie].edges[split_edge][0];
    hi_nd = mod->grid->elem3d[ie].edges[split_edge][1];

    // refines the edge if needed
    nd1 = mod->grid->elem3d[ie].nodes[hi_nd];
    nd2 = mod->grid->elem3d[ie].nodes[lo_nd];
    edge_pntr = edge_hash_lookup(nd1, nd2, edge_hashtab);
    new_node = adpt_get_node(mod, edge_pntr, edge_hashtab);
    
#ifdef _ADH_GROUNDWATER
    //printf("\n\n old : %i, new : %i, max : %i", mod->grid->nelems3d_old, mod->grid->nelems3d, mod->grid->max_nelems3d);
    sgw_3d_elem_realloc_init(mod->sgw, mod->grid->max_nelems3d, mod->grid->nelems3d, mod->nalloc_inc);
#endif
    // gets the new element
    new_elem = elem3d_new(mod->grid, mod->nalloc_inc, mod->grid->elem3d[ie].nnodes);
    
    // copies the current element to the new element
    selem3d_copy(&(mod->grid->elem3d[new_elem]), mod->grid->elem3d[ie]);
    
    // corrects the current elemental data
    mod->grid->elem3d[ie].nodes[lo_nd] = new_node;
    mod->grid->elem_error[ie] = 0.0;
    for (i=0; i<mod->grid->elem3d[ie].nnodes; i++) {
        snode_copy(&(nodes[i]), mod->grid->node[mod->grid->elem3d[ie].nodes[i]]);
    }
    // CJT :: below only for tets!
    get_tet_linear_djac_gradPhi(&(mod->grid->elem3d[ie]), nodes, NULL);
    
    // corrects the new elements data
    mod->grid->elem3d[new_elem].nodes[hi_nd] = new_node;
    mod->grid->elem_error[new_elem] = 0.0;
    for (i=0; i<mod->grid->elem3d[ie].nnodes; i++) {
        snode_copy(&(nodes[i]), mod->grid->node[mod->grid->elem3d[new_elem].nodes[i]]);
    }
    // CJT :: below only for tets!
    get_tet_linear_djac_gradPhi(&(mod->grid->elem3d[new_elem]), nodes, NULL);

    // get new element adaption level
    int new_node_level = UNSET_INT;
    if (forced == YES) {  // leave the level alone
        new_node_level = elem3d_level(mod->grid->elem3d[ie]);
    } else {
        new_node_level = elem3d_level(mod->grid->elem3d[ie]) + 1;
    }
    
    mod->grid->elem3d[ie].levels[lo_nd] = new_node_level;
    mod->grid->elem3d[new_elem].levels[hi_nd] = new_node_level;
    
#ifdef _ADH_GROUNDWATER
    /* Note that grid->elem3d is realloc'd in elem3d_new() above. Only after that can we do the following: */
    /* Make sure lo_nd and hi_nd are correctly used below! */
    int new_node_in_orig_elem = lo_nd;
    int new_node_in_new_elem = hi_nd;
    //sgw_3d_elem_realloc_init(mod->sgw, mod->grid->max_nelems3d, mod->grid->nelems3d_old);
    sgw_3d_elem_avg(mod->sgw, ie, new_elem, new_node_in_orig_elem, new_node_in_new_elem, mod->grid);
#endif
}
