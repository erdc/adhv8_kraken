/*! \file  fr_bedshstr_drag_coef.c This file collections functions responsible for computinf friction coefficient  */
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief This routine calculates a coefficient c_f, such that the shear stress tau is found as follows:
 *  \f$ \tau = 1/2 \rho c_f \|V\|^2 \f$
 *  This coefficient is for depth integrated model velocity.
 *  
 *  @param[in] depth (double) - water depth
 *  @param[in] roughness_height (double) - roughness height
 *  \returns double friction coefficient \f$c_f\f$
 *  
 *  \note The minimum ratio of depth to roughness_height is 1./29.7
 *  \author Brown Jul 24 2006
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double fr_bedshstr_drag_coef(
  /* bed shear stress drag coefficient */
  double depth,			/* local depth, L */
  double roughness_height	/* the roughness height, L */
)
{
  double depth_rh;		/* the depth divided by roughness_height ratio */
  double bta;
  double phi;
  double c_f;			/* the coefficient of friction */

  depth = MAX(depth, 1.E-8);
  roughness_height = MAX(roughness_height, 1.E-8);

  if(roughness_height > 29.7 * depth)
    depth_rh = 1. / 29.7;
  else
    depth_rh = depth / roughness_height;

  bta = 29.7 * depth_rh;
  phi = pow((.4 * bta) / ((bta + 1) * (log(bta + 1) - 1) + 1), 2.);
  c_f = 2. * phi;

  return (c_f);
}
