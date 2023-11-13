/* moves nodes between processors during repartitioning */

#include "global_header.h"

#ifdef _MESSG
void partition_transfer(SMODEL *mod) {
    
    SGRID *g = mod->grid;
    int i;			/* loop counters */
    int isd;			/* loop counter over the subdomains */
    int ierr;
    int *nnode_out;		/* the number of nodes being sent to each processor */
    int **nodes_out;		/* the outgoing nodes */
    int *messg_flag;		/* yes/no flag to indicate if a node message is going
                             to the given processor */
    SNODE *old_node_pair;	/* the node numbers before repartitioning */
    int nnode = g->nnodes;
    int old_nnode;		/* number of nodes before repartitioning */
    int npes = g->smpi->npes;
    int myid = g->smpi->myid;
    int *part = g->smpi->partition_info;
    
    /* calculate auxiliary partition arrays */
    nnode_out = (int *)tl_alloc(sizeof(int), npes);
    nodes_out = (int **)tl_alloc(sizeof(int *), npes);
    
    for(isd = 0; isd < npes; isd++) nnode_out[isd] = 0;
    //tag(MPI_COMM_WORLD);
    for(i = 0; i < g->my_nnodes; i++){
        if(part[i] != myid) nnode_out[part[i]]++;
    }
    //tag(MPI_COMM_WORLD);
    for(isd = 0; isd < npes; isd++) {
        if(nnode_out[isd] != 0) nodes_out[isd] = (int *)tl_alloc(sizeof(int), nnode_out[isd]);
    }
    //tag(MPI_COMM_WORLD);
    for(isd = 0; isd < npes; isd++) nnode_out[isd] = 0;
    for(i = 0; i < g->my_nnodes; i++) {
        if(part[i] != myid){
            nodes_out[part[i]][nnode_out[part[i]]++] = i;
        }
    }
    //tag(MPI_COMM_WORLD);
    old_node_pair = (SNODE *) tl_alloc(sizeof(SNODE), g->nnodes);
    
    // char str[25];
    //  sprintf(str,"prenums_%10.5f", mod->t_prev);
    //  print_parents_to_file(mod->grid, str);
    /* save node numbers before renumbering
     (needed for storage in hash table in comm_node_data so
     that local numbers for adjacent nodes can be retrieved
     when resetting adjacent nodes after renumbering) */
    for(i = 0; i < nnode; i++) {
        old_node_pair[i].resident_id = g->node[i].resident_id;
        old_node_pair[i].resident_pe = g->node[i].resident_pe;
    }
    //tag(MPI_COMM_WORLD);
    //print_grid_to_file(mod->grid,"PRE_Comm_Nums");
    /* get new node numbers (node_pair) for departing nodes */
    old_nnode = nnode;
    comm_new_node_nums(nnode_out, nodes_out, npes, myid, mod);
    //tag(MPI_COMM_WORLD);
    //sprintf(str,"preupdate_%10.5f", mod->t_prev);
    //    print_parents_to_file(mod->grid, str);
    //print_grid_to_file(mod->grid,"PRE_UPDATE_GN");
    /* update ghost locations with the new node numbers */
    comm_update_GN(0, g);
    //tag(MPI_COMM_WORLD);
    //print_grid_to_file(mod->grid,"PRE_NODE_DATA");
    /* move the nodal data */
    //sprintf(str,"prenodedata_%10.5f", mod->t_prev);
    //    print_parents_to_file(mod->grid, str);
    
    comm_node_data(nnode_out, nodes_out, old_node_pair, old_nnode, mod);
    //tag(MPI_COMM_WORLD);
    //sprintf(str,"postnodedata_%10.5f", mod->t_prev);
    //    print_parents_to_file(mod->grid, str);
    /* free memory */
    for(isd = 0; isd < npes; isd++) {
        if(nnode_out[isd] != 0) {
            nodes_out[isd] = (int *)tl_free(sizeof(int), nnode_out[isd], nodes_out[isd]);
        }
    }
    //tag(MPI_COMM_WORLD);
    nnode_out = (int *)tl_free(sizeof(int), npes, nnode_out);
    nodes_out = (int **)tl_free(sizeof(int *), npes, nodes_out);
    old_node_pair = (SNODE *) tl_free(sizeof(SNODE), old_nnode, old_node_pair);
    
}
#else
void partition_transfer(
                        void
                        )
{
}
#endif
