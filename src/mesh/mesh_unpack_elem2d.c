/* unpacks the 2D element buffers */
/* JRC: + unpack merged element data structure for mass-conserved computation */

#include "global_header.h"

#ifdef _MESSG
void mesh_unpack_elem2d(
  NODE_LIST_ITEM ** node_hashtab,	/* tops of the linked lists for the node hash table */
  ELEM2D_LIST_ITEM ** elem2d_hashtab,	/* tops of the linked lists for the element hash table */
  MESSG_BUFFER * elem2d_ibuffer,	/* the buffer for the 2D elements */
  MESSG_BUFFER * elem2d_dbuffer, /* the buffer for the 2D elements */
  SMODEL *mod
)
{
  int ibuff, dbuff;			/* counter in the buffer */
  int nelem2d_idatum; /* amount of integer data per 2d element */
  int nelem2d_ddatum = 21;
#ifdef UNREF_CONSV
  int total_int, *movein_elem, nel, i;
  UNREFelem_size;   /*** get size_char and size_int ***/
  CONSVelem_size;
#endif
  
  /* set the packet size */
  nelem2d_idatum = elem2d_pack_icnt();
  
  /* unpack the buffer */

  for (ibuff = 0, dbuff = 0; ibuff < elem2d_ibuffer->nitem; ibuff += nelem2d_idatum, dbuff += nelem2d_ddatum)
    {
      elem2d_unpacki(((int *) elem2d_ibuffer->buffer) + ibuff, ((double *) elem2d_dbuffer->buffer) + dbuff, node_hashtab, elem2d_hashtab, mod);
    }
  
}
#else
void mesh_unpack_elem2d(
  void
)
{
}
#endif
