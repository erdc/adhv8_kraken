/* this routine gives the count of double items per node to be packed */

#include "global_header.h"


int node_pack_dcnt(SMODEL *mod)
{
    int nd_ddatum = 3;            /* the number of doubles */
    
    /* set the size of the node double datum - look at node_packd to count the data */
    if (mod->flag.SW2_FLOW) {
        nd_ddatum += 11;
        if (mod->flag.TRANSPORT) {
            nd_ddatum += (11 * mod->ntransport);
        }
        if (mod->flag.WAVE) { /* cjt 7/2012 */
            nd_ddatum += 5;   /* stress.xx, stress.xy, stress.yy, height, period, angle, number, break, speed, edb, edf (cjt) */
        }
        if (mod->flag.WIND) { /* cjt 7/2012 */
            nd_ddatum += 2;   /* VECT2D winds.x, winds.y (cjt) */
        }
#ifdef _SEDIMENT
        if(mod->sed != NULL){
            nd_ddatum += 28;
            nd_ddatum += 5*mod->nconti;
            nd_ddatum += 28*(mod->nsed + (mod->nlayers*2));
            nd_ddatum += 16*mod->nlayers;
        }
#endif
        
    }
    if (mod->flag.SW3_FLOW) {
        nd_ddatum += 29;
        if (mod->flag.TRANSPORT) {    /* cjt */
            nd_ddatum += (11 * mod->ntransport);    /* nodal_source, nodal_decay_coeff, concX3 */
        }
#ifdef _SEDIMENT
        if(mod->sed != NULL){
            nd_ddatum += 2;
            nd_ddatum += 6*(mod->nsed);
        }
#endif
    }
    
    if (mod->flag.GW_FLOW) {
        nd_ddatum += 6;
    }
    
    /* returns the number of items */
    return nd_ddatum;
}
