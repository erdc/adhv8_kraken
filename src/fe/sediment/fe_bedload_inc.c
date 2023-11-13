/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Increments the bed load solutions after a Newton iterate.
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Gary Brown, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[inout]  mod (SMODEL *) pointer to the model struct
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_bedload_inc(SMODEL *mod) {
    
  int i=UNSET_INT, inode=UNSET_INT;
  
  for (i = 0; i < mod->grid->nnodes_bed; i++) {
      inode = i;
      if (mod->grid->ndim == 3) {
          inode = mod->grid->nodeID_2d_to_3d_bed[i];
      }
      mod->sed->bedload[mod->ised].c[i] += mod->sol[inode];
  }
}
