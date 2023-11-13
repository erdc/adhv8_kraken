/* This routine calculates the DW-GW coupled residuals */
#include "global_header.h"

/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/

void fe_dwgw_hybrid_resid(SSUPER_MODEL *sm, int idum) {
    int DEBUG = OFF;
    int idof,ndof;
    
    for (idof=0, ndof = sm->nnodes * sm->max_nsys; idof<ndof; idof++){
        sm->residual[idof] = 0.0;
    }
    
    /* Process the groundwater part */
    fe_gw_resid(sm,0);
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
        int proc_count;
        for(proc_count=0;proc_count<sm->supersmpi->npes;proc_count++){
            if(proc_count==sm->supersmpi->myid){
                printf("***********myid %d",sm->supersmpi->myid);
#endif
                printScreen_resid("Final DW-GW residual", sm->residual, sm->nnodes, sm->nsys, __LINE__, __FILE__);
#ifdef _MESSG
            }
            messg_barrier(MPI_COMM_WORLD);
        }
        //exit(-1);
#endif
    }
#endif

    /* Process the diffusive wave part */
    sm->submodel[0].grid->ndim=2;    /* Hacked to 2 */
    fe_diffusive_resid(sm,0);
    sm->submodel[0].grid->ndim=3;    /* Reset to 3 */

    
    /*****************************************************************************************/
    /*****************************************************************************************/
    
#ifdef _DEBUG
    if (DEBUG) {
#ifdef _MESSG
        int proc_count;
        for(proc_count=0;proc_count<sm->supersmpi->npes;proc_count++){
            if(proc_count==sm->supersmpi->myid){
                printf("***********myid %d",sm->supersmpi->myid);
#endif
                printScreen_resid("Final DW-GW residual", sm->residual, sm->nnodes, sm->nsys, __LINE__, __FILE__);
#ifdef _MESSG
            }
            messg_barrier(MPI_COMM_WORLD);
        }
#endif
        //exit(-1);
    }
#endif
    
}
