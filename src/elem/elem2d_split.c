/* this routine splits a 2d elem */

#include "global_header.h"

#ifdef UNREF_CONSV
void average_conserved_values(int,int,int,int,int,int,int);
#endif

void elem2d_split(
                  SMODEL *mod,
                  int ielem,			/* the element to be split */
                  EDGE_LIST_ITEM ** edge_hashtab	/* the hash table of edges */
)
{
    int i;                        /*loop counter*/
    int lo_nd;			/* the low node on the split edge */
    int hi_nd;			/* the high node on the split edge */
    int new_node;			/* the new node */
    int new_node_level;		/* level for the new node */
    int split_edge;		/* the split edge */
    int new_elem;			/* the new element */
    int nd1, nd2;			/* the nodes on an edge */
    EDGE_LIST_ITEM *edge_pntr;	/* pointer to the edge */
#ifdef UNREF_CONSV
    int j, k, ierr;
    UNREF_elem *elemunref = unref_obj.elemunrefPtr;
    CONSV_elem *e12 = NULL, *e23 = NULL, *elemconsv = consv_obj.elemconsvPtr;
#endif
    SELEM_2D *elem2d = mod->grid->elem2d;
    
    /* selects the edge to be split */
    if (elem2d[ielem].nnodes == 3) {
        
        if( tl_long_edge(mod->grid->node,
                         elem2d[ielem].nodes[elem2d[ielem].edges[0][0]],
                         elem2d[ielem].nodes[elem2d[ielem].edges[0][1]],
                         elem2d[ielem].nodes[elem2d[ielem].edges[1][0]],
                         elem2d[ielem].nodes[elem2d[ielem].edges[1][1]], edge_hashtab) == 1) {
            split_edge = 0;
        } else {
            split_edge = 1;
        }
        
        if(tl_long_edge(mod->grid->node,
                        elem2d[ielem].nodes[elem2d[ielem].edges[split_edge][0]],
                        elem2d[ielem].nodes[elem2d[ielem].edges[split_edge][1]],
                        elem2d[ielem].nodes[elem2d[ielem].edges[2][0]],
                        elem2d[ielem].nodes[elem2d[ielem].edges[2][1]], edge_hashtab) == 2) {
            split_edge = 2;
        }
        
        /* sets the hi and lo nodes on the edge */
        lo_nd = elem2d[ielem].edges[split_edge][0];
        hi_nd = elem2d[ielem].edges[split_edge][1];
        
    } else {
        
        if(tl_long_edge(mod->grid->node,
                        elem2d[ielem].nodes[elem2d[ielem].edges[0][0]],
                        elem2d[ielem].nodes[elem2d[ielem].edges[0][1]],
                        elem2d[ielem].nodes[elem2d[ielem].edges[1][0]],
                        elem2d[ielem].nodes[elem2d[ielem].edges[1][1]], edge_hashtab) == 1) {
            split_edge = 0;
        } else {
            split_edge = 1;
        }
        
        if(tl_long_edge(mod->grid->node,
                        elem2d[ielem].nodes[elem2d[ielem].edges[split_edge][0]],
                        elem2d[ielem].nodes[elem2d[ielem].edges[split_edge][1]],
                        elem2d[ielem].nodes[elem2d[ielem].edges[2][0]],
                        elem2d[ielem].nodes[elem2d[ielem].edges[2][1]], edge_hashtab) == 2) {
            split_edge = 2;
        }
        
        if(tl_long_edge(mod->grid->node,
                        elem2d[ielem].nodes[elem2d[ielem].edges[split_edge][0]],
                        elem2d[ielem].nodes[elem2d[ielem].edges[split_edge][1]],
                        elem2d[ielem].nodes[elem2d[ielem].edges[3][0]],
                        elem2d[ielem].nodes[elem2d[ielem].edges[3][1]], edge_hashtab) == 2) {
            split_edge = 3;
        }
        
        /* sets the hi and lo nodes on the edge */
        lo_nd = elem2d[ielem].edges[split_edge][0];
        hi_nd = elem2d[ielem].edges[split_edge][1];
        
    }
    
    
    
    /* refines the edge if needed */
    nd1 = elem2d[ielem].nodes[hi_nd];
    nd2 = elem2d[ielem].nodes[lo_nd];
    edge_pntr = edge_hash_lookup(nd1, nd2, edge_hashtab);
    new_node = adpt_get_node(mod, edge_pntr, edge_hashtab);
    
    /* gets the new element */
    new_elem = elem2d_new(mod->grid, mod->nalloc_inc, mod->grid->elem2d[ielem].nnodes);
    
    /* since we reallocated, must redefine pointer .. weird .. */
    elem2d = mod->grid->elem2d;
    
    /* copies the current element to the new element */
    selem2d_copy(&(elem2d[new_elem]), elem2d[ielem]);
    
    /* calculate level for the new node */
    new_node_level = elem2d_level(elem2d[ielem]) + 1;
    
    /* corrects the current elements data */
    mod->grid->elem_error[ielem] = 0.;
    elem2d[ielem].nodes[lo_nd] = new_node;
    elem2d[ielem].levels[lo_nd] = new_node_level;
    // get 2D element djacs, normals and basis function gradients (elemental constants)
    if (mod->grid->elem2d[ielem].nnodes == NDONTRI) {
        get_triangle_linear_djac_nrml_gradPhi(&(mod->grid->elem2d[ielem]), mod->grid->node, NULL);
    } else {
        SVECT nds[mod->grid->elem2d[ielem].nnodes];
        for (i=0; i<mod->grid->elem2d[ielem].nnodes; i++) {
            nds[i].x = mod->grid->node[ mod->grid->elem2d[ielem].nodes[i] ].x;
            nds[i].y = mod->grid->node[ mod->grid->elem2d[ielem].nodes[i] ].y;
            nds[i].z = mod->grid->node[ mod->grid->elem2d[ielem].nodes[i] ].z;  // initial displacement should be added here
        }
        mod->grid->elem2d[ielem].nrml = get_elem2d_normals(nds);
    }
    
    /* corrects the new elements data */
    mod->grid->elem_error[new_elem] = 0.;
    elem2d[new_elem].id = new_elem;
    elem2d[new_elem].nodes[hi_nd] = new_node;
    elem2d[new_elem].levels[hi_nd] = new_node_level;
    // get 2D element djacs, normals and basis function gradients (elemental constants)
    if (mod->grid->elem2d[new_elem].nnodes == NDONTRI) {
        get_triangle_linear_djac_nrml_gradPhi(&(mod->grid->elem2d[new_elem]), mod->grid->node, NULL);
    } else {
        SVECT nds[mod->grid->elem2d[new_elem].nnodes];
        for (i=0; i<mod->grid->elem2d[new_elem].nnodes; i++) {
            nds[i].x = mod->grid->node[ mod->grid->elem2d[new_elem].nodes[i] ].x;
            nds[i].y = mod->grid->node[ mod->grid->elem2d[new_elem].nodes[i] ].y;
            nds[i].z = mod->grid->node[ mod->grid->elem2d[new_elem].nodes[i] ].z;  // initial displacement should be added here
        }
        mod->grid->elem2d[new_elem].nrml = get_elem2d_normals(nds);
    }
    
    
#ifdef UNREF_CONSV
    if (unref_obj.UNREFhash_tablePtr != NULL)
    {
        i = Wash_Hash_find_data(unref_obj.UNREFhash_tablePtr,ielem);
        if (i >= 0)
        {
            Wash_Hash_delete (unref_obj.UNREFhash_tablePtr, ielem);
            elemunref[i].ielem = UNSET_INT;
            average_conserved_values(lo_nd,hi_nd,i,new_node,ielem,new_elem,UNREF);
        }
    }
    if (consv_obj.CONSVhash_tablePtr != NULL)
    {
        i = Wash_Hash_find_data(consv_obj.CONSVhash_tablePtr, ielem);
        if (i >= 0)
        {
            Wash_Hash_delete (consv_obj.CONSVhash_tablePtr,ielem);
            elemconsv[i].ielem = UNSET_INT;
            average_conserved_values(lo_nd,hi_nd,i,new_node,ielem,new_elem,CONSV);
        }
    }
#endif
}

#ifdef UNREF_CONSV
void average_conserved_values(int node1,int node2,int unref_entry,int new_node,int elem1,int elem2,int elem2d_type)
{
    int ie, i, j, k;
    int tot_comp = 1 + 3*VEL_CONSV + ntransport;
    int iem[NDPRFC];
    VECT2D elem_vel[NDPRFC];
    double *ustar;
    CONSV_elem *elemconsv = consv_obj.elemconsvPtr;
    UNREF_elem *elemunref = unref_obj.elemunrefPtr;
    
    if (SW2_FLOW && elem2d_type == UNREF) {
        if (unref_obj.orig_nelunref == unref_obj.nelem2dunref) {
            elemunref = (UNREF_elem *)tl_realloc(sizeof(UNREF_elem),
                                                 unref_obj.orig_nelunref+1,unref_obj.orig_nelunref,elemunref);
            unref_obj.elemunrefPtr = elemunref;
            unref_obj.orig_nelunref += 1;
        }
        for(i=0;i<NDPRFC;i++)
            iem[i] = i;
        
        elemunref[unref_entry].ielem = elem1;
        if (Wash_Hash_insert (&unref_obj.UNREFhash_tablePtr,unref_entry,elemunref[unref_entry].ielem) == -105)
            printf("WASH_HASH_DUPLICATE key returned in elem2d_split.c\n");
        
        elemunref[unref_obj.nelem2dunref].ielem = elem2;
        elemunref[unref_obj.nelem2dunref].consv_ielem[0] = UNSET_INT;
        elemunref[unref_obj.nelem2dunref].consv_ielem[1] = UNSET_INT;
        elemunref[unref_obj.nelem2dunref].unref_local[0] = UNSET_INT;
        elemunref[unref_obj.nelem2dunref].unref_local[1] = UNSET_INT;
        elemunref[unref_obj.nelem2dunref].edge_num = UNSET_INT;
        elemunref[unref_obj.nelem2dunref].unref_value = NULL;
        elemunref[unref_obj.nelem2dunref].old_unref_value = NULL;
        elemunref[unref_obj.nelem2dunref].ustar[0] = (double *)tl_alloc(sizeof(double),tot_comp*NDPRFC);
        memcpy((elemunref + unref_obj.nelem2dunref)->ustar[0], (elemunref + unref_entry)->ustar[0],
               sizeof(double)*tot_comp*NDPRFC);
        elemunref[unref_obj.nelem2dunref].ustar[1] = (double *)tl_alloc(sizeof(double),tot_comp*NDPRFC);
        memcpy((elemunref + unref_obj.nelem2dunref)->ustar[1], (elemunref + unref_entry)->ustar[1],
               sizeof(double)*tot_comp*NDPRFC);
        if (Wash_Hash_insert (&unref_obj.UNREFhash_tablePtr,unref_obj.nelem2dunref,elemunref[unref_obj.nelem2dunref].ielem) == -105)
            printf("WASH_HASH_DUPLICATE key returned in elem2d_split.c\n");
        
        for (j = 0; j < 2; j++) {
            for(k = 0,ustar=elemunref[unref_entry].ustar[j]; k < 1+1*VEL_CONSV+ntransport; k++,ustar+=NDPRFC) {
                if (k == 0) {
                    ustar[node1] = (ustar[node1] + ustar[node2]) * 0.5;
                } else if (k == 1) {
                    ustar[2*node1]     = (ustar[2*node1] + ustar[2*node2]) * 0.5;
                    ustar[2*node1 + 1] = (ustar[2*node1 + 1] + ustar[2*node2 + 1]) * 0.5;
                } else {
                    if (k == 2) ustar += 2*NDPRFC; /** skip 3 double for vz **/
                    ustar[node1] = (ustar[node1] + ustar[node2]) * 0.5;
                }
            }
        }
        
        for (j = 0; j < 2; j++) {
            for(k = 0,ustar=elemunref[unref_obj.nelem2dunref].ustar[j]; k < 1+1*VEL_CONSV+ntransport; k++,ustar+=NDPRFC) {
                if (k == 0) {
                    ustar[node2] = (ustar[node1] + ustar[node2]) * 0.5;
                } else if (k == 1) {
                    ustar[2*node2]     = (ustar[2*node1] + ustar[2*node2]) * 0.5;
                    ustar[2*node2 + 1] = (ustar[2*node1 + 1] + ustar[2*node2 + 1]) * 0.5;
                } else {
                    if (k == 2) ustar += 2*NDPRFC; /** skip 3 double for vz **/
                    ustar[node2] = (ustar[node1] + ustar[node2]) * 0.5;
                }
            }
        }
        unref_obj.nelem2dunref += 1;
    } else if (SW2_FLOW && elem2d_type == CONSV) {
        if (consv_obj.orig_nelconsv == consv_obj.nelem2dconsv) {
            elemconsv = (CONSV_elem *)tl_realloc(sizeof(CONSV_elem),
                                                 consv_obj.orig_nelconsv+1,consv_obj.orig_nelconsv,elemconsv);
            consv_obj.elemconsvPtr = elemconsv;
            consv_obj.orig_nelconsv += 1;
        }
        for(i=0;i<NDPRFC;i++)
            iem[i] = i;
        
        elemconsv[unref_entry].ielem = elem1;
        if (Wash_Hash_insert (&consv_obj.CONSVhash_tablePtr,unref_entry,elemconsv[unref_entry].ielem) == -105)
            printf("WASH_HASH_DUPLICATE key returned in elem2d_split.c\n");
        
        elemconsv[consv_obj.nelem2dconsv].ielem = elem2;
        elemconsv[consv_obj.nelem2dconsv].ustar = (double *)tl_alloc(sizeof(double),tot_comp*NDPRFC);
        memcpy((elemconsv + consv_obj.nelem2dconsv)->ustar, (elemconsv + unref_entry)->ustar,
               sizeof(double)*tot_comp*NDPRFC);
        if (Wash_Hash_insert (&consv_obj.CONSVhash_tablePtr,consv_obj.nelem2dconsv,elemconsv[consv_obj.nelem2dconsv].ielem) == -105)
            printf("WASH_HASH_DUPLICATE key returned in elem2d_split.c\n");
        
        for(k = 0,ustar=elemconsv[unref_entry].ustar; k < 1+1*VEL_CONSV+ntransport; k++,ustar+=NDPRFC) {
            if (k == 0) {
                ustar[node1] = (ustar[node1] + ustar[node2]) * 0.5;
            } else if (k == 1) {
                ustar[2*node1]     = (ustar[2*node1] + ustar[2*node2]) * 0.5;
                ustar[2*node1 + 1] = (ustar[2*node1 + 1] + ustar[2*node2 + 1]) * 0.5;
            } else {
                if (k == 2) ustar += 2*NDPRFC; /** skip 3 double for vz **/
                ustar[node1] = (ustar[node1] + ustar[node2]) * 0.5;
            }
        }
        
        for(k = 0,ustar=elemconsv[consv_obj.nelem2dconsv].ustar; k < 1+1*VEL_CONSV+ntransport; k++,ustar+=NDPRFC) {
            if (k == 0) {
                ustar[node2] = (ustar[node1] + ustar[node2]) * 0.5;
            } else if (k == 1) {
                ustar[2*node2]     = (ustar[2*node1] + ustar[2*node2]) * 0.5;
                ustar[2*node2 + 1] = (ustar[2*node1 + 1] + ustar[2*node2 + 1]) * 0.5;
            } else {
                if (k == 2) ustar += 2*NDPRFC; /** skip 3 double for vz **/
                ustar[node2] = (ustar[node1] + ustar[node2]) * 0.5;
            }
        }
        consv_obj.nelem2dconsv += 1;
    }
}
#endif
