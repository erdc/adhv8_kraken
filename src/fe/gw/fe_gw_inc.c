/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Increments the GW solutions after a Newton iterate.
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
void fe_gw_inc(SSUPER_MODEL *sm, int imod) {
#ifndef _PETSC

    SMODEL *mod = &(sm->submodel[imod]);

    int i, i1;
    double factor = 1.;
    /*mwf debug 
    printf("in fe_gw_inc\n");
    */
    for (i = 0; i < mod->grid->nnodes; i++) {
      i1 = mod->fmap[i]; /* mod->fmap[i] * 1  */
      /*printf("Node %d phead: %30.20e inc %30.20e\n",i,mod->sgw->gw_phead[i],mod->sol[i]);*/
      mod->sgw->gw_phead[i]  += factor * sm->sol[i1];
    }
    
    ///exit(-1);
#endif
}
