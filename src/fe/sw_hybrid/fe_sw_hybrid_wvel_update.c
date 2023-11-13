// This is routine updates the horizontal velocity and displacement solution variables for the 3D SW equations
#include "global_header.h"

/***********************************************************/
/***********************************************************/
static int DEBUG = OFF;
/***********************************************************/
/***********************************************************/

void fe_sw_hybrid_wvel_update(SSUPER_MODEL *sm, int idum) {
    int i;
    SMODEL *mod;
    for (i=0; i<sm->nsubmodels; i++){
        mod = &(sm->submodel[i]);
        if (mod->flag.SW3_FLOW){
            fe_wvel_update(sm,i);
        }
    }
}
