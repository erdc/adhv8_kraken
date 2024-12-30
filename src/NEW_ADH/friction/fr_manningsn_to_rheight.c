/*! \file  fr_manningsn_to_rheight.c This file collections functions responsible for computing friction coefficient  */
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief This routine converts a value of Manning's n to a value of the equivalent sand roughness height.
 *  \param[in] n (double) - manning's n L^(1/3) T
 *  \param[in] kn (double) - Manning's units conversion constant L^(1/6)
 *  \param[in] g (double) - gravitational constant
 *  \returns ks (double) - equivalent sand roughness height
 *  
 *  \author Brown Jul 24 2006
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double fr_manningsn_to_rheight(
  /* convert Manning's n to roughness height */
  double n,			/* Manning's n, L^(1/3)T */
  double kn,    /* Manning's units constant L^(1/6) */
  double g			/* the acceleration due to gravity, LT^(-2) */
)
{

  double ks;			/* the equivalent sand roughness height */

  if (n < 1.E-8)
  {
   ks = 0.0; 
  }
  else
  {
    ks = pow ((8.25 * sqrt (g) * n / kn), 6.);
  }

  return (ks);
}
