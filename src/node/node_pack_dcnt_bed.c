/* this routine gives the count of double items per node to be packed */

#include "global_header.h"


int node_pack_dcnt_bed(SMODEL *mod)
{
  int nd_ddatum = 0;            /* the number of doubles */

  /* set the size of the node double datum - look at node_packd to count the data */
  if (mod->flag.SW3_FLOW) {
    nd_ddatum += 3;
  }
#ifdef _SEDIMENT
    if(mod->sed != NULL){
    nd_ddatum += 22;
    nd_ddatum += 5*mod->nconti;
    nd_ddatum += mod->nsed*(17+ (mod->nlayers*2));
    nd_ddatum += 12*mod->nlayers;
    }
#endif
  
  /* returns the number of items */
  return nd_ddatum;
}
