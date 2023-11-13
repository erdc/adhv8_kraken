#include "global_header.h"

/* This updates the GW struct. It needs to be called when ghost nodes have been added after element renumbering */

void comm_update_gw(SGW *gw, SGRID *grid) {
    
    int i,j,k; //counters
    double *packed_gw;
    int num_doubles_gw = 6;
    
    packed_gw = (double *)tl_alloc(sizeof(double), grid->nnodes * num_doubles_gw);
    
    /* pack the array */
    for(j=0,i=0;i<grid->nnodes;i++){
        packed_gw[j++] = gw->gw_phead[i];
        packed_gw[j++] = gw->old_gw_phead[i];
        packed_gw[j++] = gw->older_gw_phead[i];
        packed_gw[j++] = gw->predict_gw_phead[i];
        packed_gw[j++] = gw->gw_density[i];
        packed_gw[j++] = gw->error[i];
    }
    
    /*update the array */
    comm_update_double(packed_gw, num_doubles_gw, grid->smpi);
    
    /* update the gw struct (only ghost nodes to save time and be safe)*/
    for(j=0,i=0;i<grid->nnodes;i++){
        gw->gw_phead[i] = packed_gw[j++];
        gw->old_gw_phead[i] = packed_gw[j++];
        gw->older_gw_phead[i] = packed_gw[j++];
        gw->predict_gw_phead[i] = packed_gw[j++];
        gw->gw_density[i] = packed_gw[j++];
        gw->error[i] = packed_gw[j++];
    }
    packed_gw = (double *)tl_free(sizeof(double), grid->nnodes * num_doubles_gw, packed_gw);
    return;
}
