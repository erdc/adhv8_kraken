/* This routine loops over the submodels to calculate the Jacobian of the shallow water supermodel */
#include "global_header.h"

/***********************************************************/
/***********************************************************/
static int DEBUG_MATRIX = OFF;

/* Gajanan gkc WARNING :*/
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* KEEP THESE VALUES SAME AS IN FE_SW_HYBRID_WVEL_RESID.C ! */
/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
static int is_wvel_dirichlet  =  YES; /* YES enforces linearly varying wvel along a node column */

/***********************************************************/
/***********************************************************/

void fe_sw_hybrid_wvel_load(SSUPER_MODEL *sm, int idum) {
    
    int i, j, k, l, iend;
    int model_3d_id=UNSET_INT;
    SMODEL *mod;
    
    if (debug.matrix_hvel == ON) {
        DEBUG_MATRIX = ON;
    }
    
    // initialize matrix
    init_adh_matrix(sm->total_nnodes, sm->max_nsys_sq, sm->matrix, sm->diagonal);
    for (i = 0; i < sm->total_nnodes; i++){
        spv_reset(&(sm->matrix[i]));
        sm->matrix[i].size = 0;
    }
    
    // GAJANAN gkc: loop over all models
    for (i = 0; i<sm->nsubmodels; i++){
        mod = &(sm->submodel[i]);
        //if(mod->proc_flag==1){
            if(mod->flag.SW3_FLOW){
                fe_wvel_load(sm,i);
            }
        //}
    }
    
#ifdef _MESSG
    //    comm_update_double(sm->diagonal, sm->nsys_sq, sm->supersmpi_wvel);
#endif
    
    for (j=0; j<sm->NumInterfaces; j++){
        if (sm->submodel[sm->interface[j].model_id[0]].flag.SW2_FLOW &&
            sm->submodel[sm->interface[j].model_id[1]].flag.SW3_FLOW){
            model_3d_id = 1;
        }
        else if (sm->submodel[sm->interface[j].model_id[0]].flag.SW3_FLOW &&
                 sm->submodel[sm->interface[j].model_id[1]].flag.SW2_FLOW){
            model_3d_id = 0;
        }
        else{  /* i.e., if this interface is 2D-2D or 3D-3D. */
            continue;
        }
        // What follows is only executed for 2D-3D coupling.
#ifdef _DEBUG
        if (DEBUG_MATRIX){
            printf("\n******************************\nModifying WVEL Jacobian for interface %i :: model_3d_id: %d\n",j+1,model_3d_id);
        }
#endif
        mod = &(sm->submodel[sm->interface[j].model_id[model_3d_id]]); /* 3d model alias */
        int *fmap = mod->fmap_wvel;
        SGRID *grid = mod->grid;
        SINTERFACE_NODELIST *ndlist;
        
        
        if (is_wvel_dirichlet == YES){
            for (k=0; k<sm->interface[j].NumNodeColumns; k++){
                ndlist = &(sm->interface[j].nodelist[k]);
                if (ndlist->size[model_3d_id] > 0){
                    for (i = 0; i<ndlist->size[model_3d_id]; i++){
                        int inode = ndlist->couplednodes[model_3d_id][i];
                        int irow = inode; //fmap[inode];

                        // Reset the matrix rows
                        spv_init(sm->matrix[irow], sm->max_nsys_sq);
                        spv_reset(&(sm->matrix[irow]));
                        sm->matrix[irow].size = 0;

                        // Set the diagonal to -1 consistent with the residual.
                        sm->diagonal[irow*sm->nsys_sq] = -1.0;

                    }
                }
            }
            mod = NULL;
        } else {
        }
    }
    
    for (i = 0, iend = sm->wvel_nnodes; i < iend; i++){     // This loop relies on the fact that nsys = 1
        if (fabs(sm->diagonal[i]) < NOT_QUITE_SMALL){
            sm->diagonal[i] = NOT_QUITE_SMALL;
        }
    }
#ifdef _MESSG
        //comm_update_double(sm->diagonal, sm->nsys_sq, sm->supersmpi_wvel);
#endif
    
    //printf("sm->supersmpi_wvel->npes: %d \t sm->wvel_nnodes: %d \t  sm->nsys_sq: %d sm->max_nsys_sq: %d\n",sm->supersmpi_wvel->npes,sm->wvel_nnodes, sm->nsys_sq, sm->max_nsys_sq);
    
#ifdef _DEBUG
    if (DEBUG_MATRIX) {
#ifdef _MESSG
        int proc_count;
        for(proc_count=0;proc_count<sm->supersmpi_wvel->npes;proc_count++){
            if(proc_count==sm->supersmpi_wvel->myid){
                printf("***********myid %d",sm->supersmpi_wvel->myid);
#endif
                printScreen_matrix("Hybrid WVEL matrix", sm->diagonal, sm->matrix, sm->wvel_nnodes, sm->nsys_sq, __LINE__, __FILE__);
#ifdef _MESSG
            }
            fflush(stdout);
            messg_barrier(sm->supersmpi_wvel->ADH_COMM);
        }
        tl_error("exit");
#endif
    }
#endif
}
