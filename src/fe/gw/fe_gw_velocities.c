
#include "global_header.h"

/*******
-k(\grad \phi_h -\rho g)
-k\grad  \phi_h + k\rho g = -k\grad \phi_h -k\rho\grad z = -k(\grad \phi + \rho \grad z)
 */

void fe_gw_elem_flux(int nnodes, TENSOR k, double *psi, double *z, double *dens,
		     SVECT3D *grad_phi, SVECT3D darcy_flux) {
  int i;
  SVECT3D grad_pot,gradient_z;
  vel.x = 0.; vel.y=0.; vel.z=0.;
  gradient_z.x=0.; gradient_z.y=0.; gradient_z.z=1.;
  grad_pot.x=0.; grad_pot.y=0.; grad_pot.z=0.;
  double avg_dens = 0.25*(dens[0]+dens[1]+dens[2]+dens[3]);
  for (i=0; i < nnodes; i++) {
    grad_pot.x += grad_phi[i].x*psi[i] + avg_dens*gradient_z.x;
    grad_pot.y += grad_phi[i].y*psi[i] + avg_dens*gradient_z.y;
    grad_pot.z += grad_phi[i].z*psi[i] + avg_dens*gradient_z.z;
  }
  VT_3D_TENS_VECT_PROD(darcy_flux,k,grad_pot);
  VT_3D_VSCALE(darcy_flux,-1.0);
}
