/* This routine computes the ev.xx, ev.yy, ev.xy and ev.yx
   using the scheme proposed by Smagorinski */

/* GSavant Nov, 2010 */

#include "stdio.h"
#include "math.h"

double tur_smag ( double, double, double, double, double, double);

double tur_smag ( 
    double dudx,
    double dudy,
    double dvdx,
    double dvdy,
    double djc,
    double smag_coeff
)
{
	double sxx;
	double syy;
	double sxy;
	double syx;
	double mag_strain;
	double turviscosity;
	double length;

/*	mag_strain = sqrt(0.5 * ((dudy * dudy) + (dvdx * dvdx) + 2 * dudy * dvdx)); */
	mag_strain = sqrt(2.*(dudx*dudx+dvdy*dvdy+(dudy+dvdx)*(dudy+dvdx)));
    length=pow(djc,0.33333333333333);
	turviscosity= (smag_coeff*length)*(smag_coeff*length)*mag_strain;
	/* turviscosity = pow (smag_coeff * pow (djc, 1./3.), 2) * mag_strain; */
	/*printf("THe Turbulent Vis is %lf \t DUDY %lf \t DVDX %lf Smag %lf mag strain %lf  \t\n",turviscosity, dudy, dvdx, smag_coeff, mag_strain);*/
	if (turviscosity < 0)
		turviscosity = 0.;
	return(turviscosity);
	
}
