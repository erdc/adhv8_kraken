/* this routine renumbers the nodes */
/* NOTE:  THE GLOBAL NODE NUMBERS NEED TO BE UPDATED AFTER CALLING THIS ROUTINE !! */

#include "global_header.h"

void node_renumber(SMODEL *mod, int flag) {
    
    int i, j, ie;
    int num_alt_nodes;              /* the total number of nodes that have been changed */
    int old_node;
    SGRID *grid = mod->grid; // alias
    int max_nnode = mod->grid->max_nnodes;
    
    // temporary storage arrays
    SVECT *vtmp = (SVECT *) tl_alloc(sizeof(SVECT), max_nnode);
    SVECT2D *v2dtmp = (SVECT2D *) tl_alloc(sizeof(SVECT2D), max_nnode);
    int *itmp = (int *) tl_alloc(sizeof(int), max_nnode);
    double *dtmp = (double *) tl_alloc(sizeof(double), max_nnode);
    int *order_tmp = (int *) tl_alloc(sizeof(int), max_nnode);
    int *new_node_number = (int *) tl_alloc(sizeof(int), max_nnode);
    
    /* form the new order */
    node_order(grid, mod->nblock, new_node_number);
    grid = mod->grid;
    
    /* Initialize the order array so we don't do it every sub call */
    for (i = 0; i < grid->max_nnodes; i++) {
        order_tmp[i]=0;
    }
    //assert(grid->type == COLUMNAR);
    if(grid->type == COLUMNAR){
        for(i=0;i<grid->max_nnodes_sur;i++){
            if(grid->nodeID_2d_to_3d_sur[i] != UNSET_INT){
                grid->old_global_surf[i]=grid->node[grid->nodeID_2d_to_3d_sur[i]].global_surf_id;
            }else{
                grid->old_global_surf[i]=UNSET_INT;
            }
        }
        
        for(i=0;i<grid->max_nnodes_bed;i++){
            if(grid->nodeID_2d_to_3d_bed[i] != UNSET_INT){
                grid->old_global_bed[i]=grid->node[grid->nodeID_2d_to_3d_bed[i]].global_bed_id;
            }else{
                grid->old_global_bed[i]=UNSET_INT;
            }
        }
    }
    
    /* check to see if any nodes have changed numbers */
    num_alt_nodes = 0;
    for (i = 0; i < max_nnode; i++) {
        if (i != new_node_number[i]) {
            num_alt_nodes++;
            if (new_node_number[i] > i) order_tmp[new_node_number[i]] = 1;
        }
    }
    
    
    //printf("myid %d num_alt_nodes: %d\n",grid->smpi->myid, num_alt_nodes);
    
    if (num_alt_nodes > 0) {
        
        /* renumber physics variables */
        if (mod->flag.SW2_FLOW)  {ssw_2d_renumber(mod->sw->d2, max_nnode, new_node_number, order_tmp, dtmp, v2dtmp);}
        if (mod->flag.SW2_FLOW && mod->flag.WAVE) {swave_renumber(mod->sw->d2->waves, max_nnode, new_node_number, order_tmp);}
        if (mod->flag.SW2_FLOW && mod->flag.WIND) {swind_renumber(mod->sw->d2->winds, max_nnode, new_node_number, order_tmp, v2dtmp);}
        if (mod->flag.SW3_FLOW)  {ssw_3d_renumber(mod->sw->d3, max_nnode, new_node_number, order_tmp, dtmp, v2dtmp, vtmp);}
        if (mod->flag.NS3_FLOW)  {sns_3d_renumber(mod->ns->d3, max_nnode, new_node_number, order_tmp, dtmp, v2dtmp, vtmp);}
#ifdef _ADH_GROUNDWATER
        if (mod->flag.GW_FLOW)  {sgw_3d_renumber(mod->sgw, max_nnode, new_node_number, order_tmp, dtmp, v2dtmp, vtmp);}
#endif
        if (mod->flag.TRANSPORT) {scon_renumber(mod->con, mod->ntransport, max_nnode, new_node_number, order_tmp, dtmp, v2dtmp);}
#ifdef _SEDIMENT
        if (mod->flag.SW2_FLOW && mod->flag.SEDIMENT)  {ssediment_renumber(mod->sed, max_nnode, new_node_number, order_tmp, itmp, dtmp, v2dtmp);}
        if (mod->flag.SW3_FLOW && mod->flag.SEDIMENT)  {ssediment_renumber_3d_sus(mod->sed, max_nnode, new_node_number, order_tmp, itmp, dtmp, v2dtmp);}
        if (mod->flag.NS3_FLOW && mod->flag.SEDIMENT)  {ssediment_renumber_3d_sus(mod->sed, max_nnode, new_node_number, order_tmp, itmp, dtmp, v2dtmp);}
#endif
        
#ifdef _ADH_ICM
        if (mod->flag.ICM) {}
#endif
        
        /* renumber anything in node structs that needs it (parents, etc) */
        node_renumber_snode(grid->max_nnodes, grid->node, new_node_number, order_tmp);
        /* renumber element connectivities */
        for (ie = 0; ie < grid->nelems3d; ie++) {
            for (i = 0; i < grid->elem3d[ie].nnodes; i++){
                if (grid->elem3d[ie].nodes[i] != UNSET_INT)
                    grid->elem3d[ie].nodes[i] = new_node_number[grid->elem3d[ie].nodes[i]];
            }
        }
        for (ie = 0; ie < grid->nelems2d; ie++) {
            for (i = 0; i < grid->elem2d[ie].nnodes; i++){
                if (grid->elem2d[ie].nodes[i] != UNSET_INT)
                    grid->elem2d[ie].nodes[i] = new_node_number[grid->elem2d[ie].nodes[i]];
            }
        }
        for (ie = 0; ie < grid->nelems1d; ie++) {
            for (i = 0; i < grid->elem1d[ie].nnodes; i++){
                if (grid->elem1d[ie].nodes[i] != UNSET_INT)
                    grid->elem1d[ie].nodes[i] = new_node_number[grid->elem1d[ie].nodes[i]];
            }
        }
    }
    
#ifdef _MESSG
    /* fixes node_info.rnode for the nodes I own */
    for (i = 0; i < grid->my_nnodes; i++) {
        grid->node[i].resident_id = i;
    }
#endif
    
    if (mod->flag.GRID_ADAPTION ) {
#ifdef _MESSG
        for (i = 0; i < grid->my_nnodes; i++) {
            if(grid->node[i].parent_res_pe[0] == grid->smpi->myid){
                old_node = grid->node[i].parent_res_id[0];
                grid->node[i].parent_res_id[0] = new_node_number[old_node];
                
            }
            if(grid->node[i].parent_res_pe[1] == grid->smpi->myid){
                old_node = grid->node[i].parent_res_id[1];
                grid->node[i].parent_res_id[1] = new_node_number[old_node];
            }
        }
#else
        for (i = 0; i < grid->nnodes; i++) {
            /* fix the first parent (adjacent) node */
            if (grid->node[i].original_id == UNSET_INT) {
                j = grid->node[i].parent[0];
                if(new_node_number[j]>=grid->nnodes || j < 0){
                    grid->node[i].parent[0] = UNSET_INT;
                }else{
                    grid->node[i].parent[0] = new_node_number[j];
                }
                j = grid->node[i].parent[1];
                if(new_node_number[j]>=grid->nnodes || j < 0){
                    grid->node[i].parent[1] = UNSET_INT;
                }else{
                    grid->node[i].parent[1] = new_node_number[j];
                }
            }
        }
#endif
    }
    //print_int_array_to_file("new_nums",new_node_number,max_nnode);
    //messg_barrier(MPI_COMM_WORLD);
    // free temporary storage arrays
    vtmp = (SVECT *) tl_free(sizeof(SVECT), max_nnode, vtmp);
    v2dtmp = (SVECT2D *) tl_free(sizeof(SVECT2D), max_nnode, v2dtmp);
    itmp = (int *) tl_free(sizeof(int), max_nnode, itmp);
    dtmp = (double *) tl_free(sizeof(double), max_nnode, dtmp);
    order_tmp = (int *) tl_free(sizeof(int), max_nnode, order_tmp);
    new_node_number = (int *) tl_free(sizeof(int), max_nnode, new_node_number);
    
}
