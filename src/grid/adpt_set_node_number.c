#include "global_header.h"

void adpt_set_node_number(SGRID *grid) {
    int n, n_old, n_new, n_strt, n_orig;
    int nnode_tot, nelem2d_tot;
    int mec, ie, nelem2d_loc;
    int myelem;
    
    /* Find Maximum Array Sizes And Allocate Temporary Arrays */
    
    nnode_tot = grid->my_nnodes;
/*#ifdef _MESSG
    ierr = MPI_Allreduce(&grid->my_nnodes, &my_nnode_max, 1, MPI_INT, MPI_MAX, ADH_COMM);
    if(ierr != MPI_SUCCESS) {
        messg_err(ierr);
    }
    ierr = MPI_Allreduce(&grid->my_nnode, &nnode_tot, 1, MPI_INT, MPI_SUM, ADH_COMM);
    if(ierr != MPI_SUCCESS) {
        messg_err(ierr);
    }
#endif*/
    
    /* Calculate Node Numbers For Adapted Nodes */
    for(n = 0, n_old = 0, n_new = 0; n < grid->my_nnodes; n++) {
        if(grid->node[n].original_id != UNSET_INT) {
            n_old++;
        } else {
            n_new++;
        }
    }
    
    n_strt = n_new;
    n_orig = n_old;
/*#ifdef _MESSG
    ierr = MPI_Scan(&n_new, &n_strt, 1, MPI_INT, MPI_SUM, ADH_COMM);
    if(ierr != MPI_SUCCESS) {
        messg_err(ierr);
    }
    ierr = MPI_Allreduce(&n_old, &n_orig, 1, MPI_INT, MPI_SUM, ADH_COMM);
    if(ierr != MPI_SUCCESS) {
        messg_err(ierr);
    }
#endif*/
    
    n_strt = n_strt - n_new;
    
    for(n = 0; n < grid->my_nnodes; n++) {
        if(grid->node[n].original_id == UNSET_INT) {
            grid->node[n].original_id  = n_orig + n_strt;
            n_strt++;
        }
    }
    
#ifdef _MESSG
  //  comm_update_int(orig_nd_number, 1); // DANG IT ... original ID needs to be either taken out of node struct or a dummy array made
    
/*    for(n = 0; n < HASHSIZE; n++)
        node_hashtab[n] = NULL;
    for(n = 0; n < nnode; n++)
        node_hash_add_entry(node_pair[n], node_hashtab, orig_nd_number[n]);
    for(n = 0; n < my_nnode; n++) {
        if(node_ladj[n].rnode != UNSET_INT && node_radj[n].rnode != UNSET_INT) {
            nd_pntr = node_hash_lookup(node_ladj[n], node_hashtab);
            if(nd_pntr == NULL)
                tl_error("Left adjacent node lookup failed in print_adpt_mesh.");
            node_ladj[n].global_num = nd_pntr->local;
            
            nd_pntr = node_hash_lookup(node_radj[n], node_hashtab);
            if(nd_pntr == NULL)
                tl_error("Right adjacent node lookup failed in print_adpt_mesh.");
            node_radj[n].global_num = nd_pntr->local;
        }
    }
    tl_list_free_all(NODE_LIST);
    */
   // comm_update_GLOBAL_NODE(2, node_ladj);
  //  comm_update_GLOBAL_NODE(2, node_radj);
#endif
    
    /* Determine Processor Owning Element (Processor Owning Minimum Node Number Of Element) */
    mec = 0;
    for(ie = 0; ie < grid->nelems2d; ie++) {
/*#ifdef _MESSG
        myelem = 1;
        for(n = 0; n < NDPRFC; n++) {
            if(node_pair[elem2d[ie].nodes[n]].sd < myid) {
                myelem = 0;
                break;
            }
        }
        if(myelem == 1) {
            grid->elem2d[ie].sd = myid;
        }
#else*/
        myelem = 1;
//#endif
        if(myelem == 1) {
            mec++;
        }
    }
    
    nelem2d_loc = mec;
    nelem2d_tot = nelem2d_loc;
/*#ifdef _MESSG
    ierr = MPI_Allreduce(&nelem2d_loc, &nelem2d_tot, 1, MPI_INT, MPI_SUM, ADH_COMM);
    if(ierr != MPI_SUCCESS)
        messg_err(ierr);
#endif*/
    
    grid->orig_initial_nnodes = grid->initial_nnodes;
    grid->orig_initial_nelems = grid->initial_nelems;
    grid->initial_nnodes = nnode_tot;
    grid->initial_nelems = nelem2d_tot;
}

//************************************************************************//
//************************************************************************//
//************************************************************************//
void adpt_unset_node_number(SGRID *grid) {
  int n;

  /* Reset values for non-adaptive output */
/*#ifdef _MESSG
  for(n = 0; n < grid->nelems2d; n++) {
      grid->elem2d[n].sd = UNSET_INT;
  }
#endif*/

  for(n = 0; n < grid->nnodes; n++) {
      if(grid->node[n].level > 0)
          grid->node[n].original_id = UNSET_INT;
/*#ifdef _MESSG
      node_ladj[n].global_num = UNSET_INT;
      node_radj[n].global_num = UNSET_INT;
#endif*/
  }
  grid->initial_nnodes = grid->orig_initial_nnodes;
  grid->initial_nelems = grid->orig_initial_nelems;
}
