/* this routine gives the count of integer items per element to be packed */

#include "global_header.h"

#ifdef _MESSG
int elem1d_pack_icnt(
  void
)
{
  /* nodes and material property */
  return (6 * NDPREDG + 5);
}
#else
int elem1d_pack_icnt(
  void
)
{
}
#endif
