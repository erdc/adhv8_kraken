/** \file
 * \brief fr_rheight_to_mannings
 * 
 * This file converts a value of the equivalent sand roughness height to a value of Manning's n.
   The conversion is calculated for Metric or English units.
   The input required is:
   double ks
   double g
   The function returns
   double n

   Brown Jul 24 2006
 *  */



#include "fr_defs.h"

 /** The distance between \f$(x_1,y_1)\f$ and \f$(x_2,y_2)\f$ is 
  \f$\sqrt{(x_2-x_1)^2+(y_2-y_1)^2}\f$. and \f$ \dwdx \f$
  */
double fr_rheight_to_manningsn(
  /* convert roughness height to Manning's n */
  double ks,			/* the equivalent sand roughness height, L */
  double g			/* the acceleration due to gravity, LT^(-2) */
)
{
  double kn;			/* the Manning's equation constant */
  double n;			/* Manning's n, L^(1/3)T */

  if (ks < 1.E-8)
  {
   n = 0.0;
  }
  else
  {
    ks = MAX (ks, 1.E-8);
    g = MAX (g, 1.E-8);

  if(g > 10.0)
    kn = 1.486;
  else
    kn = 1.0;

    n = kn / (8.25 * sqrt (g) / pow (ks, (1. / 6.)));
  }

  return (n);
}
