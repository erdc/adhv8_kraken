/*! \file  fr_urv_drag_coef.c This file calculates a coefficient c_f, such that the shear stress may be computed  */
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief This routine calculates a coefficient c_f, such that the shear stress tau is found as follows:
 *  \f$ \tau = 1/2 \rho c_f vel*vel \f$
 *  \param[in] depth (double) - water depth
 *  \param[in] roughness_height (double) - the effective roughness height of the stems, L
 *  \param[in] stem_diameter (double) - the average stem diameter, L
 *  \param[in] stem_density (double) - the average stem density (stems per unit area), L^-2
 *  \returns c_f (double) - friction coefficient
 *  \note The minimum ratio of depth to grain roughness_height is 1.00
 *  \author Brown Jul 24 2006
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double fr_urv_drag_coef(
  /* unsubmerged rigid vegetation drag coefficient */
  double depth,			/* local depth, L */
  double roughness_height,	/* the roughness height, L */
  double stem_diameter,		/* the average stem diameter, L */
  double stem_density		/* the average stem density (stems per unit area), L^-2 */
)
{
  double depth_rh;		/* the depth divided by roughness_height ratio */
  double cd;			/* drag coefficient for stems */
  double c_f;			/* the coefficient of friction */
  double c1;
  double c2;
  double c3;

  depth = MAX(depth, 1.E-8);
  roughness_height = MAX(roughness_height, 1.E-8);

  if(roughness_height > depth)
    depth_rh = 1.;
  else
    depth_rh = depth / roughness_height;

  cd = 0.4;

  c1 = 0.32 * (1. - stem_density * PI * stem_diameter * stem_diameter / 4.);
  c2 = log(10.94 * depth_rh + 1.);
  c3 = cd * depth * stem_density * stem_diameter;

  c_f = c1 / (c2 * c2) + c3;

  return (c_f);
}
