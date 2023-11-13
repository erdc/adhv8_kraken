/* This routine computes the ev.xx, ev.yy, ev.xy and ev.yx
   using the scheme proposed by Webel and Schatzmann (1984) */

/* GBrown Jun, 2013 */

#include "stdio.h"
#include "math.h"

double tur_ws (double, double, double);

double tur_ws ( 
    double vmag,
    double depth,
    double mix_coeff
)
{
	double turviscosity;
	double drag;

       /* assume depth to roughness ratio of 100 calculate the drag coefficient from Christensen = 0.0033 */
         drag = 0.0066;
         turviscosity = mix_coeff * 0.92 * depth * sqrt(drag) * vmag; 
         return(turviscosity);
	
}
