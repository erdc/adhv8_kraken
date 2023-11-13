/* unpacks the 1D element buffers */

#include "global_header.h"

#ifdef _MESSG
void mesh_unpack_elem1d(
  NODE_LIST_ITEM ** node_hashtab,	/* tops of the linked lists for the node hash table */
  ELEM1D_LIST_ITEM ** elem1d_hashtab,	/* tops of the linked lists for the element hash table */
  MESSG_BUFFER * elem1d_ibuffer,	/* the buffer for the 1D elements */
  MESSG_BUFFER * elem1d_dbuffer, /* the buffer for the 1D elements */
  SMODEL *mod
)
{
  int ibuff, dbuff;			/* counter in the buffer */
  int nelem1d_idatum;		/* amount of integer data per 1d element */
  int nelem1d_ddatum = 11;   /* amount of double data per 1d element x,y,z of nodes*/

  /* set the packet size */
  nelem1d_idatum = elem1d_pack_icnt();

  /* unpack the buffer */
  for(ibuff = 0, dbuff = 0; ibuff < elem1d_ibuffer->nitem; ibuff += nelem1d_idatum, dbuff += nelem1d_ddatum)
    elem1d_unpacki(((int *)elem1d_ibuffer->buffer) + ibuff, ((double *)elem1d_dbuffer->buffer) + dbuff, node_hashtab, elem1d_hashtab, mod);
}
#else
void mesh_unpack_elem1d(
  void
)
{
}
#endif
