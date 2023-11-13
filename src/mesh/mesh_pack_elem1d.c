#include "global_header.h"

/*! 
   \brief Pack 1D elements as indicated by the node list 
 */
void mesh_pack_elem1d(int *pack_node,   /* the new partition of the nodes I own */
                      MESSG_BUFFER * ibuffer, /* the integer buffer to be packed */
                      MESSG_BUFFER * dbuffer, /* the double buffer to be packed x, y, z for each node*/
                      SGRID *grid
  )
{
  int ie = 0;                   /* loop counter over the elements */
  int ibuff = 0, dbuff = 0;                /* position in buffer */
  int nelem1d_out = 0;          /* the number of 1d elements going out */
  int nelem1d_idatum = 0;       /* amount of integer data per 1d element */
  int nelem1d_ddatum = 11;       /* amount of integer data per 1d element */

  /* Checking */
  assert(pack_node != NULL);
  assert(ibuffer != NULL);

  /* Count the number of elements going out */
  for (ie = 0, nelem1d_out = 0; ie < grid->nelems1d; ie++)
    {
       if (pack_node[grid->elem1d[ie].nodes[0]] == YES || pack_node[grid->elem1d[ie].nodes[1]] == YES)
            {
              nelem1d_out++;
            }
        
    }

  /* Allocate the element buffers */
  nelem1d_idatum = elem1d_pack_icnt();
  messg_buffer_alloc(nelem1d_out * nelem1d_idatum, sizeof(int), ibuffer);
  ibuffer->type = MESSG_INT;
  messg_buffer_alloc(nelem1d_out * nelem1d_ddatum, sizeof(double), dbuffer);
  dbuffer->type = MESSG_DOUBLE;

  /* Load the element buffers */
  for (ie = 0, ibuff = 0, dbuff = 0; ie < grid->nelems1d; ie++)
    {
     
          if (pack_node[grid->elem1d[ie].nodes[0]] == YES || pack_node[grid->elem1d[ie].nodes[1]] == YES)
            {
              elem1d_packi(((int *) ibuffer->buffer) + ibuff, ie, grid);
              ibuff += nelem1d_idatum;
              elem1d_packd(((double *) dbuffer->buffer) + dbuff, ie, grid);
              dbuff += nelem1d_ddatum;
            }
        
    }
  return;
}
