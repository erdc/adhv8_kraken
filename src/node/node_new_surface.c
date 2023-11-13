/* this routine gets a new node */

#include "global_header.h"

int node_new_surface(SMODEL *mod) {
    
    int the_node;                 /* the node */
    int i;               /* loop counters */
    
   
    /* checks to see if it is time to allocate more nodes */
    if (mod->grid->nnodes_sur == mod->grid->max_nnodes_sur) {
   
        mod->grid->max_nnodes_sur += mod->nalloc_inc;
        if (mod->grid->type == COLUMNAR) mod->grid->old_global_surf = tl_realloc(sizeof(int), mod->grid->max_nnodes_sur, mod->grid->old_max_nnodes_sur, mod->grid->old_global_surf);
        mod->grid->nodeID_2d_to_3d_sur = tl_realloc(sizeof(int), mod->grid->max_nnodes_sur, mod->grid->old_max_nnodes_sur, mod->grid->nodeID_2d_to_3d_sur);
#ifdef _MESSG
        if (mod->grid->type == COLUMNAR)  mod->grid->smpi->surface_partition_info = tl_realloc(sizeof(int), mod->grid->max_nnodes_sur, mod->grid->old_max_nnodes_sur, mod->grid->smpi->surface_partition_info);
#endif
        /* Gajanan gkc - possible bug fix - comparing with node_new_bed.c. Corey, please check this which one of node_new_bed and node_new_surface is correct! */
        //for(i=mod->grid->nnodes_sur;i<mod->grid->max_nnodes_sur;i++) { /* gkc */
        for(i=mod->grid->old_max_nnodes_sur;i<mod->grid->max_nnodes_sur;i++) { /* Original */
          if (mod->grid->type == COLUMNAR) mod->grid->old_global_surf[i]=UNSET_INT;
          mod->grid->nodeID_2d_to_3d_sur[i]=UNSET_INT;
        }
        mod->grid->old_max_nnodes_sur = mod->grid->max_nnodes_sur;

        //printf("to  %d\n ", mod->grid->max_nnodes_sur);   
        /* reallocate and initalize physics arrays for more nodes */
        if (mod->flag.SW3_FLOW) {ssw_3d_realloc_init_surface(mod->sw->d3, mod->grid->nnodes_sur, mod->grid->max_nnodes_sur);}
        //if (mod->flag.NS3_FLOW) {sns_3d_realloc_init_surface(mod->ns->d3, mod->grid->nnodes_sur, mod->grid->max_nnodes_sur);}
//#ifdef _ADH_GROUNDWATER
//        if (mod->flag.GW_FLOW) {sgw_3d_realloc_init_surface(mod->sgw, mod->grid->nnodes_sur, mod->grid->max_nnodes_sur);}
//#endif
#ifdef _ADH_ICM
        if (mod->flag.ICM) {}
#endif
    }
    
    the_node = mod->grid->nnodes_sur;
    mod->grid->nnodes_sur++;
    return (the_node);
}
