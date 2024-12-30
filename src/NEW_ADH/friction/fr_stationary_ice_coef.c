/*! \file  fr_stationary_ice_coef.c This file calculates a coefficient c_f for stationary ice coefficient of friction, such that the shear stress may be computed  */
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief This routine calculates a coefficient c_f for stationary ice coefficient, such that the shear stress tau is found as follows:
 *  \f$ \tau = 1/2 \rho c_f vel*vel \f$
 *  \param[in] depth (double) - water depth
 *  \param[in] velocity (double) - velocity magnitude
 *  \param[in] rheight_bed (double) - bed roughness height
 *  \param[in] rheight_ice (double) - ice roughness height
 *  \param[in] rho (double) - density
 *  \param[in] which (int) - which value to return (1 = bed shear, 2 = ice shear, 3 = total shear)
 *  \returns c_f (double) - friction coefficient
 *  \author Brown Jul 24 2006
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double fr_stationary_ice_coef(
  double depth,			/* local depth */
  double velocity,		/* velocity magnitude */
  double rheight_bed,		/* bed roughness height */
  double rheight_ice,		/* ice roughnes height */
  double rho,			/* density */
  int which			/* which value to return (1 = bed shear, 2 = ice shear, 3 = total shear) */
)
{
  double alpha_bed, beta_bed, z_mv, z_mv1;
  double alpha_ice, beta_ice, cd_bed, cd_ice;
  double A, alpha_ibr, kappa;
  double coeff_bed, coeff_ice, coeff_total;
  double shear, shear_bed, shear_ice;
  double ustar_bed, ustar_ice, ev_bed, ev_ice, shear_exchange;
  double delta_mv;
  int i;

  kappa = 0.4;
  delta_mv = 0.368;
  z_mv1 = 0.5 * depth;
  z_mv = 0.5 * depth;
  if (rheight_bed > depth) rheight_bed = depth;
  if (rheight_ice > depth) rheight_ice = depth;

  for(i = 0; i < 100; i++)
    {
      beta_bed = 29.7 * (z_mv1 / rheight_bed) + 1.0;
      beta_ice = 29.7 * ((depth - z_mv1) / rheight_ice) + 1.0;
      A = log(beta_bed) / log(beta_ice);
      z_mv1 = depth / (A * A + 1.0);
      if (fabs (z_mv1 - z_mv) < 1.e-3 * depth)
	break;
      z_mv = z_mv1;
    }

   beta_bed = 29.7 * (z_mv / rheight_bed) + 1.0;
   beta_ice = 29.7 * ((depth - z_mv) / rheight_ice) + 1.0;

  cd_bed =
    2. * pow(((kappa * (beta_bed - 1.)) / (beta_bed * (log(beta_bed) - 1.) + 1.)), 2.);
  cd_ice =
    2. * pow(((kappa * (beta_ice - 1.)) / (beta_ice * (log(beta_ice) - 1.) + 1.)), 2.);

  alpha_ibr = pow(((cd_bed / cd_ice) * ((depth - z_mv) / z_mv)), 0.5);
  alpha_bed = 1. / (alpha_ibr * ((depth - z_mv) / depth) + z_mv / depth);
  alpha_ice = 1. / ((depth - z_mv) / depth + z_mv / (alpha_ibr * depth));

  coeff_bed = cd_bed * pow(alpha_bed, 2.);
  coeff_ice = cd_ice * pow(alpha_ice, 2.);
  coeff_total = coeff_bed + coeff_ice;

  shear_bed = .5 * rho * cd_bed * pow((alpha_bed * velocity), 2.);
  shear_ice = .5 * rho * cd_ice * pow((alpha_ice * velocity), 2.);
  shear = shear_bed + shear_ice;

 /* finally, correct for cross-profile exchange of energy */

  ustar_bed = sqrt(shear_bed/rho);
  ustar_ice = sqrt(shear_ice/rho);
  ev_bed = 0.25 * kappa * ustar_bed * z_mv;
  ev_ice = 0.25 * kappa * ustar_ice * (depth - z_mv);
  shear_exchange = velocity * rho * 0.5 * (ev_bed + ev_ice) *
                   (alpha_bed - alpha_ice) /
                   (depth * (1. - delta_mv));

  shear_bed += shear_exchange;
  shear_ice -= shear_exchange;

  if (velocity > 1.e-8) coeff_bed = shear_bed / (.5 * rho *  pow (velocity, 2.));
  if (velocity > 1.e-8) coeff_ice = shear_ice / (.5 * rho *  pow (velocity, 2.));
                   
  switch (which)
    {
    case BED_CO:
      return coeff_bed;
    case ICE_CO:
      return coeff_ice;
    case TOT_CO:
      return coeff_total;
    default:
      return ERROR;
    }
}
