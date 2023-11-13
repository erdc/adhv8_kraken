/* Version x.x */

#include "global_header.h"

void flap_flow(SNODE *node, int nnodes, int nflap, double *ol_head, SVECT2D *ol_vel, 
	           STR_VALUE *str_values, SFLAP_C *flap, SFLAGS flag, double gravity, SGRID *grid)
  {
  int i;
  for(i = 0; i < nflap; i++) 
   {
     flap_compound_flux(i, node, nnodes, ol_head, ol_vel, nflap, str_values, flap, gravity, grid);
   }
}
