#include "global_header.h"

void sinterface_alloc_init(SINTERFACE *ifce, int NumNodeColumns, int m2d_id, int m3d_id) {
    int i;
    ifce->NumNodeColumns = NumNodeColumns;
    ifce->model_id[0] = m2d_id;
    ifce->model_id[1] = m3d_id;
    ifce->nodelist = (SINTERFACE_NODELIST *) tl_alloc(sizeof(SINTERFACE_NODELIST), NumNodeColumns);
    for (i=0; i<NumNodeColumns; i++){
        ifce->nodelist[i].surfnode[0] = UNSET_INT;
        ifce->nodelist[i].surfnode[1] = UNSET_INT;
        ifce->nodelist[i].size[0] = 0;
        ifce->nodelist[i].size[1] = 0;
        ifce->nodelist[i].save_wvel_residual = NULL;
    }
}

void sinterface_free(SINTERFACE *ifce) {
    int i;
    if (ifce != NULL){
        for (i=0; i<ifce->NumNodeColumns; i++) {
            if(ifce->nodelist[i].size[0]>0){
                ifce->nodelist[i].couplednodes[0] = (int *) tl_free(sizeof(int), ifce->nodelist[i].size[0], ifce->nodelist[i].couplednodes[0]);
            }
            if(ifce->nodelist[i].size[1]>0){
                ifce->nodelist[i].couplednodes[1] = (int *) tl_free(sizeof(int), ifce->nodelist[i].size[1], ifce->nodelist[i].couplednodes[1]);
            }
            if (ifce->nodelist[i].save_wvel_residual != NULL){
                ifce->nodelist[i].save_wvel_residual = (double *) tl_free(sizeof(double), ifce->nodelist[i].wvel_residual_length, ifce->nodelist[i].save_wvel_residual);
            }
        }
        ifce->nodelist = (SINTERFACE_NODELIST *) tl_free(sizeof(SINTERFACE_NODELIST), ifce->NumNodeColumns, ifce->nodelist);
    }
}
