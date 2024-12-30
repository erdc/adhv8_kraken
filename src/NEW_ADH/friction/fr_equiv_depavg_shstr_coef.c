/*! \file  fr_equiv_depavg_shstr_coef.c This file collections functions responsible for computing friction coefficient  */
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief This routine calculates a coefficient \f$c_{eda}\f$, such that the shear stress \f$tau\f$ is found as follows:
 *  \f$ \tau = c_{eda} \tau_{nbv}\f$
 *  where \f$\tau_{nbv}\f$ is the shear stress as calcuated using near bed velocity.
 *  Hence, this routine yields a coefficient that will convert the shear stress found
 *  using near bed velocity into an approximate shear stress that would result from using
 *  depth averaged valocity.
 *  \param[in] depth (double) - water depth 
 *  \param[in] roughness_height (double) - roughness height
 *  \returns c_eda (double) 
 *  
 *  \note The minimum ratio of depth to roughness_height is 1./29.7
 *  \author Brown Jul 24 2006
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
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
