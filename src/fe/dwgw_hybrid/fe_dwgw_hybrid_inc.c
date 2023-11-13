/* This routine loops over submodels to increment the shallow water solution */
#include "global_header.h"

/***********************************************************/
/***********************************************************/

void fe_dwgw_hybrid_inc(SSUPER_MODEL *sm, int idum) {

    int DEBUG_LOCAL = OFF;
    
    /* Process the groundwater part */
    fe_gw_inc(sm,0);

    /* Process the diffusive wave part */
    fe_diffusive_inc(sm,0);
    
#ifdef _DEBUG
    /*Do nothing for now*/
#endif
}
