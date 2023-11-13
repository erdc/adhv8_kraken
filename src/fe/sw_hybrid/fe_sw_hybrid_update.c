/* 3D : This is routine updates the horizontal velocity and displacement solution variables for the 3D shallow water equations */
/* 2D : This is routine updates the shallow water variables */

#include "global_header.h"

/***********************************************************/
/***********************************************************/

static int DEBUG = OFF;

/***********************************************************/
/***********************************************************/

void fe_sw_hybrid_update(SSUPER_MODEL *sm, int idum) {
    int i;	    		/* loop counter */
    // aliases
    SMODEL *mod;
    for (i=0; i<sm->nsubmodels; i++){
        mod = &(sm->submodel[i]);
      if(mod->proc_flag==1){
        if (mod->flag.SW2_FLOW){
            fe_sw2_update(sm,i);
        }
        else if (mod->flag.SW3_FLOW){
            fe_hvel_update(sm,i);
        }
      }
    }
}
