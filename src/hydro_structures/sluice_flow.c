/* Version x.x */

#include "global_header.h"

void sluice_flow(SNODE *node, int nnodes, int nsluice, double *ol_head, SVECT2D *ol_vel, STR_VALUE *str_values, SSERIES *series_head, SSLUICE_C *sluice, SFLAGS flag, double gravity, SGRID *grid)
{
  int i;
  for(i = 0; i < nsluice; i++) 
   { 
     sluice_compound_flux(i, node, nnodes, ol_head, ol_vel, nsluice, str_values, series_head, sluice, gravity, grid);
   }
}
