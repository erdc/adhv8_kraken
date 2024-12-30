/*! \file  fr_sav_drag_coef.c This file calculates a coefficient c_f, such that the shear stress may be computed  */
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief TThis routine calculates a coefficient c_f, such that the shear stress tau is found as follows:
 *  \f$ \tau = 1/2 \rho c_f vel*vel \f$
 *  \param[in] depth (double) - water depth
 *  \param[in] roughness_height (double) - the effective roughness height of the stems, L
 *  \param[in] sav_stem_height (double) - the undeflected sav stem height, L
 *  \returns c_f (double) - friction coefficient
 *  
 *  \author Brown Jul 24 2006
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
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
