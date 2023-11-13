/* this routine gets a new node */

#include "global_header.h"

int node_new(SMODEL *mod) {
    
    int the_node;                 /* the node */
    int i;               /* loop counters */
    
    /* checks to see if it is time to allocate more nodes */
    if (mod->grid->nnodes == mod->grid->max_nnodes) {
      //printf("************myid %d adding nodes from %d", mod->grid->smpi->myid, mod->grid->max_nnodes);
        mod->grid->max_nnodes += mod->nalloc_inc;
        //printf("to  %d\n ", mod->grid->max_nnodes);
        mod->grid->node = (SNODE *) tl_realloc(sizeof(SNODE), mod->grid->max_nnodes, mod->grid->nnodes, mod->grid->node);
        if(mod->grid->ndim == 3){
          mod->grid->nodeID_3d_to_2d_sur = (int*) tl_realloc(sizeof(int), mod->grid->max_nnodes, mod->grid->nnodes, mod->grid->nodeID_3d_to_2d_sur);
          mod->grid->nodeID_3d_to_2d_bed = (int*) tl_realloc(sizeof(int), mod->grid->max_nnodes, mod->grid->nnodes, mod->grid->nodeID_3d_to_2d_bed);
        }
       
#ifdef _MESSG
       mod->grid->smpi->partition_info = (int *) tl_realloc(sizeof(int), mod->grid->max_nnodes, mod->grid->nnodes, mod->grid->smpi->partition_info);
#endif
        // initializes new nodes
        for (i = mod->grid->nnodes; i < mod->grid->max_nnodes; i++) {
            snode_init(&mod->grid->node[i]);
            mod->grid->node[i].id = i;
            if(mod->grid->ndim == 3){
              mod->grid->nodeID_3d_to_2d_sur[i]=UNSET_INT;
              mod->grid->nodeID_3d_to_2d_bed[i]=UNSET_INT;
            }
        }
        
        /* reallocate and initalize physics arrays for more nodes */
        if (mod->flag.SW2_FLOW) {ssw_2d_realloc_init(mod->sw->d2, mod->grid->nnodes, mod->grid->max_nnodes);}
        if (mod->flag.SW3_FLOW) {ssw_3d_realloc_init(mod->sw->d3, mod->grid->nnodes, mod->grid->max_nnodes);}
        if (mod->flag.NS3_FLOW) {sns_3d_realloc_init(mod->ns->d3, mod->grid->nnodes, mod->grid->max_nnodes);}
#ifdef _ADH_GROUNDWATER
        if (mod->flag.GW_FLOW) {sgw_3d_realloc_init(mod->sgw, mod->grid->nnodes, mod->grid->max_nnodes);}
#endif
        if (mod->flag.TRANSPORT) {scon_realloc_init(mod->con, mod->ntransport, mod->grid->nnodes, mod->grid->max_nnodes);}
#ifdef _SEDIMENT
        if (mod->flag.SW2_FLOW && mod->flag.SEDIMENT) {
          ssediment_realloc_init_bed(mod->sed, mod->grid->nnodes, mod->grid->max_nnodes);
          ssediment_realloc_init_sus(mod->sed, mod->grid->nnodes, mod->grid->max_nnodes);
        }
        if (mod->flag.SW3_FLOW && mod->flag.SEDIMENT) {
          ssediment_realloc_init_sus(mod->sed, mod->grid->nnodes, mod->grid->max_nnodes);
        }
        if (mod->flag.NS3_FLOW && mod->flag.SEDIMENT) {
            ssediment_realloc_init_sus(mod->sed, mod->grid->nnodes, mod->grid->max_nnodes);
        }
#endif

#ifdef _ADH_ICM
        if (mod->flag.ICM) {}
#endif
    }
    
    the_node = mod->grid->nnodes;
    mod->grid->nnodes++;
    
    /* returns the new node */
    return (the_node);
}
