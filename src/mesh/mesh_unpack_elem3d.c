/* unpacks the 3D element buffers */

#include "global_header.h"

#ifdef _MESSG
void mesh_unpack_elem3d(
  NODE_LIST_ITEM ** node_hashtab,	/* tops of the linked lists for the node hash table */
  ELEM3D_LIST_ITEM ** elem3d_hashtab,	/* tops of the linked lists for the element hash table */
  MESSG_BUFFER * elem3d_ibuffer,	/* the buffer for the 3D elements */
  MESSG_BUFFER * elem3d_dbuffer,  /* the buffer for the 3D elements */
  SMODEL *mod
)
{
  int ibuff, dbuff;			/* counter in the buffer */
  int nelem3d_idatum;		/* amount of integer data per 3d element */
  int nelem3d_ddatum = 26; /* x y z coordinates for all 4 nodes */

  /* set the packet size */
  nelem3d_idatum = elem3d_pack_icnt();

  /* unpack the buffer */
  for(ibuff = 0, dbuff = 0; ibuff < elem3d_ibuffer->nitem; ibuff += nelem3d_idatum, dbuff += nelem3d_ddatum)
    elem3d_unpacki(((int *)elem3d_ibuffer->buffer) + ibuff, ((double *)elem3d_dbuffer->buffer) + dbuff, node_hashtab, elem3d_hashtab, mod);
}
#else
void mesh_unpack_elem3d(
  void
)
{
}
#endif
