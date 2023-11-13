/* This evaluates the 2nd order Mellor-Yamada */
/* turbulent viscosity and vertical diffusion*/
/* Based on " Development of Turbulence Closure*/
/* Model for Geophysical Fluid Problems*/
/* Reviews of Geophysics and Space Physics Vol 28, No 4, 851-825*/
/* Mellor-Yamada*/

/* "Some effects of Suspended Sediment Stratification on an Oceanic Bottom Boundary Layer*/
/* Adams and Weatherly*/
/* J. of Geo physical Research Vol 86, No C5*/

/* Includes adjustments from Henderson-Sellers 1984*/
/* New Formulation for Eddy Diffusion Thermocline Models*/

/* GSavant Nov, 2010 */

#include "global_header.h"


double tur_MY_2 ( double grav,
    double elem_depth,
    double dist_above_bed,
    double dudz,
    double dvdz,
    double drdz,
    double ave_den,
	int supression,
    int want
)
{
	
	/* A1,A2,B1,B2 and C2 are Empirical Coefficients*/
	/* Already determined by Mellor-Yamada in '82 publication*/
	double A1 = 0.92;
	double B1 = 16.6;
	double A2 = 0.74;
	double B2 = 10.1;
	double C1 = 0.08;
	double gamma1;
	double gamma2;
	double Shh;
	double Smm;
	double Ghh;
	double Gmm;
	double mix_l;
	double kappa = 0.4;   /* Von Karman*/
	double grad_rich;
	double tke;
	double tur_evxz;
	double tur_evyz;
	double tur_diff;
	
	gamma1 = (1./3.) - (2*A1/B1);
	gamma2 = (B2/B1)  + (6*A1/B1);
	Shh = (3.* A2) * (gamma1); /* Stabilization param for salt*/
	Smm = (Shh * A1 / A2) * (B1 * (gamma1 - C1)/(B1 * gamma1));/* Stabilization Param for momentum*/

	mix_l = kappa * dist_above_bed*(fabs(1-dist_above_bed/elem_depth));

	Ghh = -( grav * drdz / ave_den);
	Gmm = dudz*dudz + dvdz*dvdz;
	/* Compute Gradient Richardson Number */
        /* HENSLEY MESSING AROUND HERE */

	if (is_double_small(Gmm) == 1) /* cjt */
        /* if (Gmm == 0.0) */
          {
            grad_rich = 0.0;
          }
        else
          {
            grad_rich = Ghh/Gmm;
          }

	tke = pow(B1 * mix_l * mix_l * Smm * Gmm, 0.5);

	/* Eddy Viscosity and turb diffusion */
	/* Added by RCB to make sure unstable stratification doesn't suppress diffusion */
    if(grad_rich<0.)grad_rich=0.;

	tur_evxz = ave_den*Smm * mix_l * tke / (1 + 0.74 * grad_rich);     /* Eqn 37 in Henderson-Sellers paper */
	tur_evyz = ave_den*Smm * mix_l * tke / (1 + 0.74 * grad_rich);
	tur_diff = Shh *  mix_l * tke / (1 + 37. * grad_rich * grad_rich);    /* Eqn 42 in Henderson-Sellers paper */

	if (tur_evxz < 0)
		tur_evxz = 0;
	if (tur_evyz < 0)
		tur_evyz = 0;
	if (tur_diff < 0)
		tur_diff = 0;
	if (want == 0)
		return(tur_evxz);
	else if (want == 1)
		return(tur_evyz);
	else if (want == 2)
		return(tur_diff);
	else
	{
		printf("\nUnknown want coefficient from tur_MY_2, returning default of 0\n");
		tur_evxz = 0;
		return(tur_evxz);
	}

}
