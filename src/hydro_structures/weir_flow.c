/* Version x.x */

#include "global_header.h"

void weir_flow(SNODE *node, int nnodes, int nweir, double *ol_head, SVECT2D *ol_vel, STR_VALUE *str_values, SWEIR_C *weir, SFLAGS flag, double gravity, SGRID *grid)
{
  int i;
  for(i = 0; i < nweir; i++) 
   { 
     weir_compound_flux(i, node, nnodes, ol_head, ol_vel, nweir, str_values, weir, gravity, grid);
   }
}
