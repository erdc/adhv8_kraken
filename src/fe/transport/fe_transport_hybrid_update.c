/* This routine loops over all submodels to call the respective transport update functions */
#include "global_header.h"

/***********************************************************/
/***********************************************************/

static int DEBUG = OFF;

/***********************************************************/
/***********************************************************/

void fe_transport_hybrid_update(SSUPER_MODEL *sm, int dum) {
    int i;	    		/* loop counter */
    // aliases
    SMODEL *mod;
    for (i=0; i<sm->nsubmodels; i++){
        if(sm->submodel[i].proc_flag==1){
            fe_transport_update(sm,i);
        }
    }
}
