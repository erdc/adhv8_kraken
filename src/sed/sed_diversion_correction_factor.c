/* the calculates a correction facctor for sediment diversions. */
/* the user defines a string across a diversion channel */
/* the user then defines a diversion invert elevation, and a main channel bed elevation */
/* (the main channel bed elevation must be less than (lower than) the invert elevation) */
/* This routine then calcualtes the fraction of the main channel sediment flux that is */
/* diverted, based on the assumtion that the water is simply "skimmed" off of the main channel */
/* (i.e. there is no convergence or divergence of vertical streamlines as the water is */
/* diverted from the main channel) */

#include "share_extern.h"

double
sed_diversion_correction_factor (
                                  double esrh,
                                  double l_bedelev,
				  double l_depth,
                                  double l_disp,	
                                  double l_top,
                                  double l_bot,
                                  double l_mainch,
				  double rdf
                                  )
{
  double k;			/* the sediment profile exponent */
  double kh;                    /* the dimensionless sediment profile exponent */
  double poly;			/* a polynomial expansion needed for the computation */
  double fact;			/* parameter used to compute factorial */
  double beta;			/* a dimensionless roughness parameter */
  double x1;
  double x2;
  double x3;
  double z;
  double z_top;
  double z_bot;
  double z_mainch;
  double z_nbed;
  double z_wsel;
  double sed_flux[4];
  double dcf;
  int i, j;			/* loop counter */

  kh = 1.0;
  if (rdf < 1.0) rdf = 1.0;
  for (i = 0; i < 100; i++)
    {
      kh = rdf * (1. - exp (-kh));
    }

  if (l_depth < (esrh/29.7 + 1.E-8)) l_depth = esrh/29.7 + 1.E-8;
  z_wsel = l_bedelev + l_depth + l_disp;
  z_mainch = l_mainch;
  z_nbed = l_mainch + esrh/29.7;
  z_bot = l_bot;
  z_top = l_top;
  if (z_bot < z_nbed) z_bot = z_nbed;
  if (z_top < z_nbed) z_top = z_nbed;
  if (z_bot > z_wsel) z_bot = z_wsel;
  if (z_top > z_wsel) z_top = z_wsel;


  for (j = 0; j < 4; j++)
  {
    if (j == 0) z = z_nbed;
    if (j == 1) z = z_bot;
    if (j == 2) z = z_top;
    if (j == 3) z = z_wsel;
    k = kh * (z - z_nbed) / (z_wsel - z_mainch);
    poly = 0.0;
    fact = 1.0;
    for (i = 0; i < 50; i++)
      {
        fact *= (i + 1);
        poly += pow (-k, i + 1) / ((i + 1) * fact);
      }
    beta = 29.7 * (z - z_mainch) / esrh;
    x1 = log (beta);
    x2 = 1. -  exp (-k);
    x3 = poly;

    sed_flux[j] = x1 * x2 + x3;
  }

  dcf = (sed_flux[2] - sed_flux[1])/(sed_flux[3] - sed_flux[0]);
  dcf *= (z_wsel - z_mainch) / (z - z_mainch);

  return (dcf);

}
