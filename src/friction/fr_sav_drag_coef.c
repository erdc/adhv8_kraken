/* ADH Version 2.0.0 6/04 */

/* This routine calculates a coefficient c_f, such that the shear stress tau is found as follows:
   tau = 1/2 rho c_f vel*vel
   This friction coefficient expresses the friction induced by submerged aquatic vegetation (SAV)
   This coefficient is for depth integrated model velocity.
   The input required is:
   double depth
   double sav_stem_height
   The function returns
   double c_f
   The minimum ratio of depth to sav_stem_height is 1.

   Brown Jul 24 2006
 */

#include "fr_defs.h"

double fr_sav_drag_coef(
  /* submerged aquatic vegetation drag coefficient */
  double depth,			/* local depth, L */
  double roughness_height,	/* the effective roughness height of the stems, L */
  double sav_stem_height /* the undeflected sav stem height, L */  
) 
{
  double depth_sh;		/* the depth divided by sav_stem_height ratio */
  double bta;
  double lda;
  double btaplda;
  double phi;
  double ts;			/*the apparent thickness of the sav */
  double rhs;			/*the apparent roughness of the sav */
  double c_f;			/* the coefficient of friction */

  depth = MAX(depth, 1.E-8);
  sav_stem_height = MAX(sav_stem_height, 1.E-8);

  if(sav_stem_height > depth)
    depth_sh = 1.;
  else
    depth_sh = depth / sav_stem_height;

  ts = .667 * depth / depth_sh;
  rhs = (roughness_height / sav_stem_height) * depth / depth_sh;
  bta = 29.7 * depth / rhs;
  lda = 1 - 29.7 * ts / rhs;
  btaplda = bta + lda;

  phi = pow((.4 * bta) / (btaplda * (log(btaplda) - 1) + 1), 2.);
  c_f = 2. * phi;

  return (c_f);
}
