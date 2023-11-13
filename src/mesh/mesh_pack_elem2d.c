#include "global_header.h"

/*!  
   \brief Pack 2D elements as indicated by the node list 
 */
void mesh_pack_elem2d(int *pack_node,   /* the new partition of the nodes I own */
                      MESSG_BUFFER * ibuffer, /* the integer buffer to be packed */
                      MESSG_BUFFER * dbuffer, /* the double buffer to be packed */
                      SGRID *g
  )
{
  int ie = 0;                   /* loop counter over the elements */
  int ibuff = 0, dbuff = 0;                /* position in buffer */
  int nelem2d_out = 0;          /* the number of 2d elements going out */
  int nelem2d_idatum = 0;       /* amount of integer data per 2d element */
  int nelem2d_ddatum = 21;       /* amount of double data per 2d element x, y, z for each node*/

  /* Checking */
  assert(pack_node != NULL);
  assert(ibuffer != NULL);

  /* Count the number of elements going out */
  for (ie = 0, nelem2d_out = 0; ie < g->nelems2d; ie++)
    {
    //  if (g->elem2d[ie].string != UNSET_INT)
    //    {
          if (pack_node[g->elem2d[ie].nodes[0]] == YES || pack_node[g->elem2d[ie].nodes[1]] == YES || pack_node[g->elem2d[ie].nodes[2]] == YES)
            {
              nelem2d_out++;
            }
      //  }
    }

  /* Allocate the element buffers */
  nelem2d_idatum = elem2d_pack_icnt();
  messg_buffer_alloc(nelem2d_out * nelem2d_idatum, sizeof(int), ibuffer);
  ibuffer->type = MESSG_INT;
  messg_buffer_alloc(nelem2d_out * nelem2d_ddatum, sizeof(double), dbuffer);
  dbuffer->type = MESSG_DOUBLE;

  /* Load the element buffers */
  for (ie = 0, ibuff = 0, dbuff = 0; ie < g->nelems2d; ie++)
    {
    //  if (g->elem2d[ie].string != UNSET_INT)
    //    {
          if (pack_node[g->elem2d[ie].nodes[0]] == YES || pack_node[g->elem2d[ie].nodes[1]] == YES || pack_node[g->elem2d[ie].nodes[2]] == YES)
            {
              
              elem2d_packi(((int *) ibuffer->buffer) + ibuff, ie, g);
              ibuff += nelem2d_idatum;
              elem2d_packd(((double *) dbuffer->buffer) + dbuff, ie, g);
              dbuff += nelem2d_ddatum;
            }
     //   }
    }
  return;
}
