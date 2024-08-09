/* ADH Version 2.0.0 6/04 */

/* This routine calculates a coefficient c_eda, such that the shear stress tau is found as follows:
   tau = c_eda tau_nbv
   where tau_nbv is the shear stress as calcuated using near bed velocity.
   Hence, this routine yields a coefficient that will convert the shear stress found
   using near bed velocity into an approximate shear stress that would result from using
   depth averaged valocity.
   The input required is:
   double depth
   double roughness_height
   The function returns
   double c_eda
   The minimum ratio of depth to roughness_height is 1./29.7

   Brown Jul 24 2006
 */

#include "fr_defs.h"

double fr_equiv_depavg_shstr_coef(
  /* equivalent depth averaged shear stress coefficient */
  double depth,			/* local depth, L */
  double roughness_height	/* the roughness height, L */
)
{
  double depth_rh;		/* the depth divided by roughness_height ratio */
  double bta;
  double c_eda;			/* the shear stress conversion coefficient */

  depth = MAX(depth, 1.E-8);
  roughness_height = MAX(roughness_height, 1.E-8);

  if(roughness_height > 29.7 * depth)
    depth_rh = 1. / 29.7;
  else
    depth_rh = depth / roughness_height;

  bta = 29.7 * depth_rh;
  c_eda = pow(log(0.368 * bta + 1.) / log(0.05 * bta + 1.), 2.);

  return (c_eda);
}
