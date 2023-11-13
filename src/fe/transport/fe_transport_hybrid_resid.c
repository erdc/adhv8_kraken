/* This routine loops over submodels to calculate residuals for the transport problem */
#include "global_header.h"

static int DEBUG = OFF;

//*****************************************************************************************//
//*****************************************************************************************//
//*****************************************************************************************//

void fe_transport_hybrid_resid(SSUPER_MODEL *sm, int dum) {
#ifndef _PETSC
    int i, idof, ndof;
    
    for (idof=0, ndof = sm->nnodes * sm->max_nsys; idof<ndof; idof++){
        sm->residual[idof] = 0.0;
    }
    
    for (i=0; i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            fe_transport_resid(sm,i);
        }
    }
    
    /* checks the residual for very small entries */
    for (i = 0; i < sm->nnodes; i++) {   /* sm->nsys = 1 here! */
        if (fabs(sm->residual[i]) < SMALL) {
            sm->residual[i] = 0.0;
        }
    }
    
    
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
        int proc_count;
        for(proc_count=0;proc_count<sm->supersmpi->npes;proc_count++){
            if(proc_count==sm->supersmpi->myid){
                printf("***********myid %d",sm->supersmpi->myid);
#endif
                printScreen_resid("Hybrid TRANSPORT residual", sm->residual, sm->nnodes, sm->nsys, __LINE__, __FILE__);
#ifdef _MESSG
            }
            messg_barrier(MPI_COMM_WORLD);
        }
        //exit(-1);
#endif
    }
#endif
#endif
}
