/* this routine gives the count of integer items per node to be packed */

#include "global_header.h"

#ifdef _MESSG
int node_pack_icnt(void)
{
  int nd_idatum = 0;            /* the number of datums */

  /* set the size of the node integer datum - look at node_new to count the data */
  nd_idatum = 19;
  /* returns the number of items */
  return (nd_idatum);
}
#else
void node_pack_icnt(void)
{
}
#endif
