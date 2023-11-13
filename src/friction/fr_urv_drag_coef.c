/* ADH Version 2.0.0 6/04 */

/* This routine calculates a friction coefficient c_f, such that the shear stress tau is found as follows:
   tau = 1/2 rho c_f vel*vel
   This friciton coefficient expresses the friction induced by rigid, unsubmerged vegetation
   This coefficient is found using depth integrated model velocity.
   The input required is:
   double depth
   double roughness_height
   double stem_diameter
   double stem_density
   The function returns
   double c_f
   The minimum ratio of depth to grain roughness_height is 1.00

   Brown Jul 24 2006
 */

#include "fr_defs.h"

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
