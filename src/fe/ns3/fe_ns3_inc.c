/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Increments the Navier Stokes solutions after a Newton iterate.
 * \author    Corey Trahan, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    David Smith, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[inout]  mod (SMODEL *) pointer to the model struct
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_ns3_inc(SSUPER_MODEL *sm, int imod) {
#ifndef _PETSC

    SMODEL *mod = &(sm->submodel[imod]);
    
    int i;
    int i4;
    double factor = 1.;
    
    for (i = 0; i < mod->grid->nnodes; i++) {
        i4 = mod->fmap[i] * 4;
        mod->ns->d3->vel[i].x += factor * sm->sol[i4];
        mod->ns->d3->vel[i].y += factor * sm->sol[i4 + 1];
        mod->ns->d3->vel[i].z += factor * sm->sol[i4 + 2];
        
        /*
        if (mod->flag.MG == ON && mod->grid->node[i].bflag == 0) {
            mod->ns->d3->displacement[i] += mod->sol[i4 + 3];
            //printf("prs: %20.10f \t dpl: %20.10f\n",mod->ns->d3->prs[i],mod->ns->d3->displacement[i]);
        } else {
            mod->ns->d3->prs[i] += mod->sol[i4 + 3];
        }
        */

        ///*
        mod->ns->d3->prs[i] += factor * sm->sol[i4 + 3];
        if (mod->flag.MG == ON && mod->grid->node[i].bflag == 0) { 
            exit(-1);
            mod->ns->d3->displacement[i] += factor * sm->sol[i4 + 2] * mod->dt;
            printf("prs: %20.10f \t dpl: %20.10f\n",mod->ns->d3->prs[i],mod->ns->d3->displacement[i]);
        }
        
        
        //*/ CJT :: TEMPORARY
        if (i==345 &&  mod->ns->d3->vel[i].x > 0.004) {printf("front moved too fast! x: %20.10f u: %20.10f\n",mod->grid->node[i].x, mod->ns->d3->vel[i].x); exit(-1);}
    }
    
    //printf("\n");
    //for (i = 0; i < mod->grid->nnodes; i++) printf("i: %d \t u: %20.10e \t v: %20.10e \t w: %20.10e \t p: %20.10e\n", i,mod->ns->d3->vel[i].x,mod->ns->d3->vel[i].y,mod->ns->d3->vel[i].z,mod->ns->d3->prs[i]);
    //exit(-1);
#endif
}
