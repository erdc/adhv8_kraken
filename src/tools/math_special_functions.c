#include "global_header.h"

// erf and erfc adapted from Numerical recipes
// erfc = 1 -erf

double erf (double x) {
    double t,z,ans;
    z=fabs(x);
    t=1.0/(1.0+0.5*z);
    
    ans=1-t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
    
    return x >= 0.0 ? ans : -ans;
}

double erfc (double x) {
    double t,z,ans;
    z=fabs(x);
    t=1.0/(1.0+0.5*z);
    ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
    return x >= 0.0 ? ans : 2.0-ans;
}

/* factorial */
double factorial(int n){
    double fact;
    int c;
 
    fact = 1;
    c = 1;
    for (c = 1; c <= n; c++) {
        fact = fact * c;
    }
    return(fact);
}

/* Bessel Function of the First Kind */
double besselj(double p, double z){
    
    double J;
    int i;
    J = 0;
    for (i = 0; i <= 10; i++) {
        J += pow(-1,i) * pow(z / 2, 2 * i + p) / (factorial(i) * factorial(i+p));
    }
    return (J);
}
