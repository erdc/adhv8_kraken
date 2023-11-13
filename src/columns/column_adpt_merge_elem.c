/* ADH Version 2.0.0 6-04 */
/* unrefines the mesh */

#include "global_header.h"

void column_adpt_merge_elem(SGRID *grid,
                            int icol,   /* the column that is being unrefined */
                            int inode,  /* the surface node that is being removed */
                            int *node_unref_flags, ELEM_REF_LIST_ITEM ** send_2d_elem,  /* linked list of 2D elements to send */
                            NODE_LIST_ITEM ** node_hashtab, /* node hash table */
                            ELEM3D_LIST_ITEM ** elem3d_hashtab, /* 3D element hash table */
                            ELEM2D_LIST_ITEM ** elem2d_hashtab  /* 2D element hash table */
)
{
    int new_elem_num;             /* the merged element number */
    int ie;                       /* loop counter over the elements */
    int current_vseg;             /* the id of the vertical list that inode is in */
    int node_count;               /* count of number of nodes in an element that lie below the
                                   surface node */
    int new_local[2];            /* the node(s) of the element that is(are) to be removed */
    int ii, jj;
    ID_LIST_ITEM *ptr;
    int isd;                      /* loop counter over the processors */
#ifdef _MESSG
    ELEM_REF_LIST_ITEM *elem_ref_list_pntr;   /* the pointer to the element linked lists */
    SNODE nd_pntr;
#endif
    int adj1, adj2;
    int ival;
    SVECT midpt[2];
    
    /* determine which vertical line of nodes we are dealing with */
    current_vseg = find_vertical_segment(grid, inode, grid->vertical_hash);
    
    /* We first loop through the column and find any elements that have 2 nodes
     * to be removed. These elements need to be merged first. */
    
    ptr = grid->column_list[icol];
    while (ptr->next != NULL) {
        ie = ptr->id;
        node_count = 0;
        new_local[0] = UNSET_INT;
        new_local[1] = UNSET_INT;
        for (ii = 0; ii < grid->elem3d[ie].nnodes; ii++) {
            if (find_vertical_segment(grid, grid->elem3d[ie].nodes[ii], grid->vertical_hash) == current_vseg) {
                /* If we have already found node, figure out which one is on top */
                if (node_count > 0) {
                    if (grid->node[grid->elem3d[ie].nodes[ii]].z > grid->node[grid->elem3d[ie].nodes[new_local[0]]].z) {
                        new_local[1] = new_local[0];
                        new_local[0] = ii;
                    } else {
                        new_local[1] = ii;
                    }
                } else {
                    new_local[0] = ii;
                }
                node_count++;
            }
        }
        if (node_count == 0) {
            tl_error("Could not find vertical string of nodes for given node.");
        }
        column_elem3d_merge(grid, ie, new_local, node_unref_flags, node_hashtab, elem3d_hashtab);
        ptr = ptr->next;
    }
    
    /* Now need to loop again and remove all of the elements
     * Note that the column list may contain references to elements that have
     * been removed (by calling elem3d_init). This isn't a problem: the column
     * list still contains a reference to at least one of the two elements to be
     * merged and that's enough to recover the entire unrefined mesh.
     */
    
    /* Gajanan gkc - prisms. I don't understand the need for this loop. Particularly if it is applicable for prisms or not. */
    
    ptr = grid->column_list[icol];
    while (ptr->next != NULL) {
        ie = ptr->id;
        /* Check to see if this element has been removed already by checking if material ID is defined */
        if (grid->elem3d[ie].mat == UNSET_INT) {
            break;
        }
        new_local[0] = UNSET_INT;
        new_local[1] = UNSET_INT;
        for (ii = 0; ii < NDONTET; ii++) {
            if (find_vertical_segment(grid, grid->elem3d[ie].nodes[ii], grid->vertical_hash) == current_vseg) {
                new_local[0] = ii;
                break;
            }
        }
        if (new_local[0] == UNSET_INT) {
            /* OOPS. This may not be an error condition, but just indicate that
             * we are dealing with an element that we have already unrefined. */
            /* tl_error("Could not find vertical string of nodes for given node."); */
            // printf("Unrefine: passing on.\n");
            ptr = ptr->next;
            continue;
        }
        column_elem3d_merge(grid, ie, new_local, node_unref_flags, node_hashtab, elem3d_hashtab);
        ptr = ptr->next;
    }
    
    /* Now, merge the 2D elements. Since there may be 2D elements on side walls, we may have
     * the case here as well that an element was cut twice. */
    
    ptr = grid->column_list2d[icol];
    while (ptr->next != NULL) {
        ie = ptr->id;
        node_count = 0;
        new_local[0] = UNSET_INT;
        new_local[1] = UNSET_INT;
        for (ii = 0; ii < grid->elem2d[ie].nnodes; ii++) {
            if (find_vertical_segment(grid, grid->elem2d[ie].nodes[ii], grid->vertical_hash) == current_vseg) {
                /* If we have already found node, figure out which one is on bottom */
                if (node_count > 0) {
                    if (grid->node[grid->elem2d[ie].nodes[ii]].z < grid->node[grid->elem2d[ie].nodes[new_local[0]]].z) {
                        new_local[1] = new_local[0];
                        new_local[0] = ii;
                    } else {
                        new_local[1] = ii;
                    }
                } else {
                    new_local[0] = ii;
                }
                node_count++;
            }
        }
        /* If we didn't find any nodes, then this is 2D element on side wall which wasn't cut */
        if (node_count > 0) {
#ifdef _MESSG
            isd = column_elem2d_merge(grid, ie, new_local, &new_elem_num, node_hashtab, elem2d_hashtab);
            if (isd != UNSET_INT) {
                elem_ref_list_pntr = (ELEM_REF_LIST_ITEM *) tl_list_alloc(ELEM2D_LIST);
                elem_ref_list_pntr->next = send_2d_elem[isd];
                elem_ref_list_pntr->ielem = new_elem_num;
                send_2d_elem[isd] = elem_ref_list_pntr;
            }
#else
            isd = column_elem2d_merge(grid, ie, new_local, &new_elem_num, elem2d_hashtab);
#endif
        }
        ptr = ptr->next;
    }
    
    /* NOTE: THIS LOOP IS ALMOST CERTAINLY UNECESSARY. The above is a loop over
     * the 2d elements 'inside' the column; that is, the surface face and the face
     * on the bottom. There can't be a second instance to remove ....
     */
    
    /* Now need to loop again and remove the second instances */
    
    /* Gajanan gkc - prisms. I don't understand the need for this loop. Particularly if it is applicable for prisms or not. */
    
    ptr = grid->column_list2d[icol];
    while (ptr->next != NULL) {
        ie = ptr->id;
        /* Check to see if this element has been removed already */
        if (grid->elem2d[ie].string == UNSET_INT) {
            break;
        }
        new_local[0] = UNSET_INT;
        new_local[1] = UNSET_INT;
        for (ii = 0; ii < NDONTRI; ii++) {
            if (find_vertical_segment(grid, grid->elem2d[ie].nodes[ii], grid->vertical_hash) == current_vseg) {
                new_local[0] = ii;
                break;
            }
        }
        if (new_local[0] == UNSET_INT) {
            /* OOPS. This may not be an error condition, but just indicate that
             * we are dealing with an element that we have already unrefined.
             */
            /* tl_error("Could not find vertical string of nodes for given node."); */
            // printf("Unrefine: passing on.\n");
            ptr = ptr->next;
            continue;
        }
#ifdef _MESSG
        isd = column_elem2d_merge(grid, ie, new_local, &new_elem_num, node_hashtab, elem2d_hashtab);
        if (isd != UNSET_INT) {
            elem_ref_list_pntr = (ELEM_REF_LIST_ITEM *) tl_list_alloc(ELEM2D_LIST);
            elem_ref_list_pntr->next = send_2d_elem[isd];
            elem_ref_list_pntr->ielem = new_elem_num;
            send_2d_elem[isd] = elem_ref_list_pntr;
        }
#else
        isd = column_elem2d_merge(grid, ie, new_local, &new_elem_num, elem2d_hashtab);
#endif
        ptr = ptr->next;
    }
    
    /* Merge the sidewall elements.
     * First, locate the edges and midpts of the segments to be merged */
    
#ifdef _MESSG
    nd_pntr.resident_pe=grid->node[inode].parent_res_pe[0];
    nd_pntr.resident_id=grid->node[inode].parent_res_id[0];
    adj1 = -1;
    adj1 = node_hash_lookup(nd_pntr, node_hashtab, grid->smpi->npes);
    if (adj1 < 0) tl_error("Adjacent node parent0 lookup failed in elem2d_merge!");
    nd_pntr.resident_pe=grid->node[inode].parent_res_pe[1];
    nd_pntr.resident_id=grid->node[inode].parent_res_id[1];
    adj2 = -1;
    adj2 = node_hash_lookup(nd_pntr, node_hashtab, grid->smpi->npes);
    if (adj2 < 0) tl_error("Adjacent node parent0 lookup failed in elem2d_merge!");
#else
    adj1 = grid->node[inode].parent[0];
    adj2 = grid->node[inode].parent[1];
#endif
    
    midpt[0].x = (grid->node[inode].x + grid->node[adj1].x) * 0.5;
    midpt[0].y = (grid->node[inode].y + grid->node[adj1].y) * 0.5;
    midpt[1].x = (grid->node[inode].x + grid->node[adj2].x) * 0.5;
    midpt[1].y = (grid->node[inode].y + grid->node[adj2].y) * 0.5;
    
    for (ii = 0; ii < 2; ii++) {
        ival = find_midpt_segment(grid, midpt[ii], grid->midpt_hash);
        if (ival >= 0) {
            ptr = grid->sidewall_list[ival];
            while (ptr->next != NULL) {
                ie = ptr->id;
                if (grid->elem2d[ie].string == UNSET_INT) {
                    break;
                }
                node_count = 0;
                new_local[0] = UNSET_INT;
                new_local[1] = UNSET_INT;
                for (jj = 0; jj < grid->elem2d[ie].nnodes; jj++) {
                    if (find_vertical_segment(grid, grid->elem2d[ie].nodes[jj], grid->vertical_hash) == current_vseg) {
                        /* If we have already found node, figure out which one is on bottom */
                        if (node_count > 0) {
                            if (grid->node[grid->elem2d[ie].nodes[jj]].z > grid->node[grid->elem2d[ie].nodes[new_local[0]]].z) {
                                new_local[1] = new_local[0];
                                new_local[0] = jj;
                            }
                            else{
                                new_local[1] = jj;
                            }
                        }
                        else {
                            new_local[0] = jj;
                        }
                        node_count++;
                    }
                }
                /* If we didn't find any nodes, then this is 2D element on side wall which wasn't cut */
                if (node_count > 0) {
#ifdef _MESSG
                    isd = column_elem2d_merge(grid, ie, new_local, &new_elem_num, node_hashtab, elem2d_hashtab);
                    /* For the sidewall elements, we don't have to worry about
                     * updating levels. The levels info is used only for the
                     * surface 2d elements
                     */
                    /*if(isd != UNSET_INT)
                     {
                     elem_ref_list_pntr = (ELEM_REF_LIST_ITEM *) tl_list_alloc(ELEM2D_LIST);
                     elem_ref_list_pntr->next = send_2d_elem[isd];
                     elem_ref_list_pntr->ielem = new_elem_num;
                     send_2d_elem[isd] = elem_ref_list_pntr;
                     }
                     */
#else
                    isd = column_elem2d_merge(grid, ie, new_local, &new_elem_num, elem2d_hashtab);
#endif
                }
                ptr = ptr->next;
            }
            
            /* Now need to loop again and remove the second instances */
            
            /* Gajanan gkc - prisms. I don't understand the need for this loop. Particularly if it is applicable for prisms or not. */
            
            ptr = grid->sidewall_list[ival];
            while (ptr->next != NULL) {
                ie = ptr->id;
                /* Check to see if this element has been removed already */
                if (grid->elem2d[ie].string == UNSET_INT) {
                    break;
                }
                new_local[0] = UNSET_INT;
                new_local[1] = UNSET_INT;
                for (jj = 0; jj < NDONTRI; jj++) {
                    if (find_vertical_segment(grid, grid->elem2d[ie].nodes[jj], grid->vertical_hash) == current_vseg) {
                        new_local[0] = jj;
                        break;
                    }
                }
                if (new_local[0] == UNSET_INT) {
                    /* OOPS. This may not be an error condition, but just indicate that
                     * we are dealing with an element that we have already unrefined.
                     */
                    /* tl_error("Could not find vertical string of nodes for given node."); */
                    // printf("Unrefine: passing on.\n");
                    ptr = ptr->next;
                    continue;
                }
#ifdef _MESSG
                isd = column_elem2d_merge(grid, ie, new_local, &new_elem_num, node_hashtab, elem2d_hashtab);
                /* For the sidewall elements, we don't have to worry about
                 * updating levels. The levels info is used only for the
                 * surface 2d elements
                 */
                /* if(isd != UNSET_INT)
                 {
                 elem_ref_list_pntr = (ELEM_REF_LIST_ITEM *) tl_list_alloc(ELEM2D_LIST);
                 elem_ref_list_pntr->next = send_2d_elem[isd];
                 elem_ref_list_pntr->ielem = new_elem_num;
                 send_2d_elem[isd] = elem_ref_list_pntr;
                 }
                 */
#else
                isd = column_elem2d_merge(grid, ie, new_local, &new_elem_num, elem2d_hashtab);
#endif
                ptr = ptr->next;
            }
        }
    }
}
