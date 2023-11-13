/* this routine gets a new node */

#include "global_header.h"

int node_new_bed(SMODEL *mod) {

    int the_node;                 /* the node */
    int i;               /* loop counters */


    /* checks to see if it is time to allocate more nodes */
    if (mod->grid->nnodes_bed >= mod->grid->max_nnodes_bed) {
        mod->grid->max_nnodes_bed += mod->nalloc_inc;
        mod->grid->old_global_bed = tl_realloc(sizeof(int), mod->grid->max_nnodes_bed, mod->grid->old_max_nnodes_bed, mod->grid->old_global_bed);
        mod->grid->nodeID_2d_to_3d_bed = tl_realloc(sizeof(int), mod->grid->max_nnodes_bed, mod->grid->old_max_nnodes_bed, mod->grid->nodeID_2d_to_3d_bed);
        /* Gajanan gkc - possible bug fix - comparing with node_new_surface.c.*/
        for(i=mod->grid->old_max_nnodes_bed;i<mod->grid->max_nnodes_bed;i++) { /* gkc */
        //for(i=mod->grid->nnodes_bed;i<mod->grid->max_nnodes_bed;i++) { /* Original */
          mod->grid->old_global_bed[i]=UNSET_INT;
          mod->grid->nodeID_2d_to_3d_bed[i]=UNSET_INT;
        }
        mod->grid->old_max_nnodes_bed = mod->grid->max_nnodes_bed;
#ifdef _SEDIMENT
        if (mod->flag.SW3_FLOW && mod->flag.SEDIMENT) {
          ssediment_realloc_init_bed(mod->sed, mod->grid->nnodes_bed, mod->grid->max_nnodes_bed);
        }
#endif
        if (mod->flag.SW3_FLOW) {ssw_3d_realloc_init_bed(mod->sw->d3, mod->grid->nnodes_bed, mod->grid->max_nnodes_bed);}
        //if (mod->flag.NS3_FLOW) {sns_3d_realloc_init_bed(mod->ns->d3, mod->grid->nnodes_bed, mod->grid->max_nnodes_bed);}
//#ifdef _ADH_GROUNDWATER
//        if (mod->flag.GW_FLOW) {sgw_3d_realloc_init_bed(mod->sgw, mod->grid->nnodes_bed, mod->grid->max_nnodes_bed);}
//#endif

    }
    the_node = mod->grid->nnodes_bed;
    mod->grid->nnodes_bed++;
    
    return (the_node);
}
