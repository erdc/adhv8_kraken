/* ADH Version 2.0.0 6/04 */

/* This routine converts a value of Manning's n to a value of the equivalent sand roughness height.
   The input required is:
   double n
   double manning's conversion constant
   double g
   The function returns
   double ks

   Brown Jul 24 2006
 */

#include "fr_defs.h"

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
