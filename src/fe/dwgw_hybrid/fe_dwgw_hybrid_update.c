/* This routine calls the parallel DW and GW update functions */
#include "global_header.h"

/***********************************************************/
/***********************************************************/

void fe_dwgw_hybrid_update(SSUPER_MODEL *sm, int idum) {

    int DEBUG_LOCAL = OFF;

    /* Process the groundwater part */
    fe_gw_update(sm,0);

    /* Process the diffusive wave part */
    fe_diffusive_update(sm,0);
}
