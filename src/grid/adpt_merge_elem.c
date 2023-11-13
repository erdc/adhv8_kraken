/* unrefines the mesh */
/* JRC: merge element and save the unrefined nodal number and the edge num it sits on */

#include "global_header.h"

void adpt_merge_elem(
                     SGRID *grid,                           /* model grid */
                     int *node_unref_flags,                 /* flags nodes to be eliminated */
#ifdef _MESSG
                     ELEM_REF_LIST_ITEM **send_3d_elem,     /* linked list of 3D elements to send */
                     ELEM_REF_LIST_ITEM **send_2d_elem,     /* linked list of 2D elements to send */
                     ELEM_REF_LIST_ITEM **send_1d_elem,     /* linked list of 1D elements to send */
                     NODE_LIST_ITEM **node_hashtab,        /* node hash table */
#endif
                     ELEM3D_LIST_ITEM **elem3d_hashtab,	    /* 3D element hash table */
                     ELEM2D_LIST_ITEM **elem2d_hashtab,	    /* 2D element hash table */
                     ELEM1D_LIST_ITEM **elem1d_hashtab      /* 1D element hash table */

) {
    
    int new_elem_num;		/* the merged element number */
    int ie;			/* loop counter over the elements */
#ifdef _MESSG
    int isd;			/* loop counter over the processors */
    ELEM_REF_LIST_ITEM *elem_ref_list_pntr;	/* the pointer to the element linked lists */
#endif
    
    /* unrefinement for 3D elements */
    for(ie = 0; ie < grid->nelems3d; ie++) {
        if(grid->elem3d[ie].mat != UNSET_INT) {
            if(node_unref_flags[grid->elem3d[ie].nodes[0]] == YES ||
               node_unref_flags[grid->elem3d[ie].nodes[1]] == YES ||
               node_unref_flags[grid->elem3d[ie].nodes[2]] == YES ||
               node_unref_flags[grid->elem3d[ie].nodes[3]] == YES) {
#ifdef _MESSG
                isd = elem3d_merge(grid, ie, &new_elem_num, node_hashtab, elem3d_hashtab, node_unref_flags);
                if(isd != UNSET_INT) {
                    elem_ref_list_pntr = (ELEM_REF_LIST_ITEM *) tl_list_alloc(ELEM3D_LIST);
                    elem_ref_list_pntr->next = send_3d_elem[isd];
                    elem_ref_list_pntr->ielem = new_elem_num;
                    send_3d_elem[isd] = elem_ref_list_pntr;
                }
#else
                elem3d_merge(grid, ie, &new_elem_num, elem3d_hashtab, node_unref_flags);
#endif
            }
        }
    }
    
    /* unrefinement for 2D elements */
    for(ie = 0; ie < grid->nelems2d; ie++) {
        if(grid->elem2d[ie].mat != UNSET_INT) {
            if(node_unref_flags[grid->elem2d[ie].nodes[0]] == YES ||
               node_unref_flags[grid->elem2d[ie].nodes[1]] == YES ||
               node_unref_flags[grid->elem2d[ie].nodes[2]] == YES) {
#ifdef _MESSG
                isd = elem2d_merge(grid, ie, &new_elem_num, node_hashtab, elem2d_hashtab, node_unref_flags);
                if (isd != UNSET_INT) {
                    elem_ref_list_pntr = (ELEM_REF_LIST_ITEM *) tl_list_alloc(ELEM2D_LIST);
                    elem_ref_list_pntr->next = send_2d_elem[isd];
                    elem_ref_list_pntr->ielem = new_elem_num;
                    send_2d_elem[isd] = elem_ref_list_pntr;
                }
#else
                elem2d_merge(grid, ie, &new_elem_num, elem2d_hashtab, node_unref_flags);
#endif
            }
        }
    }  /*** end of 2d unrefined element loop ***/
    
    /* unrefinement for 1D elements */
    for(ie = 0; ie < grid->nelems1d; ie++) {
        if(grid->elem1d[ie].string != UNSET_INT) {
            if(node_unref_flags[grid->elem1d[ie].nodes[0]] == YES ||
               node_unref_flags[grid->elem1d[ie].nodes[1]] == YES) {
#ifdef _MESSG
                isd = elem1d_merge(grid, ie, &new_elem_num, node_hashtab, elem1d_hashtab, node_unref_flags);
                if(isd != UNSET_INT) {
                    elem_ref_list_pntr = (ELEM_REF_LIST_ITEM *) tl_list_alloc(ELEM1D_LIST);
                    elem_ref_list_pntr->next = send_1d_elem[isd];
                    elem_ref_list_pntr->ielem = new_elem_num;
                    send_1d_elem[isd] = elem_ref_list_pntr;
                }
#else
                elem1d_merge(grid, ie, &new_elem_num, elem1d_hashtab, node_unref_flags);
#endif
            }
        }
    }
}
