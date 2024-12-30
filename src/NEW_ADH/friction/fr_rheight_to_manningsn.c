/*! \file  fr_rheight_to_manningsn.c This file converts a value of the equivalent sand roughness height to a value of Manning's n.  */
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief This file converts a value of the equivalent sand roughness height to a value of Manning's n.
 *  The conversion is calculated for Metric or English units.
 *  \param[in] ks (double) - the equivalent sand roughness height, L
 *  \param[in] g (double) - gravitational constant
 *  \returns n (double) - Manning's n, L^(1/3)T 
 *  
 *  \author Brown Jul 24 2006
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
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
