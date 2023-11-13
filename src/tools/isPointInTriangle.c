#include "global_header.h"

//------------------------------------------------------------------------------------
// if stand-alone uncomment below ----------------------------------------------------
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <float.h>
//#define SMALL DBL_EPSILON   /* the smallest double */
//#define YES 1           /* yes */
//#define NO -3           /* no */
//
//int compare_double( double v1, double v2) {
//    if ( fabs(v1 - v2) < SMALL) return YES;
//    return NO;
//}
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

/* A utility function to calculate area of triangle formed by (x1, y1), (x2, y2) and (x3, y3) */
double area(double x1, double y1, double x2, double y2, double x3, double y3) {
    return fabs((x1*(y2-y3) + x2*(y3-y1)+ x3*(y1-y2))/2.0);
}

/* A function to check whether point P(x, y) lies inside the triangle formed by A(x1, y1), B(x2, y2) and C(x3, y3) */
int isPointInTriangle(double x1, double y1, double x2, double y2, double x3, double y3, double x, double y) {
    /* Calculate area of triangle ABC */
    double A = area (x1, y1, x2, y2, x3, y3);
    
    /* Calculate area of triangle PBC */
    double A1 = area (x, y, x2, y2, x3, y3);
    
    /* Calculate area of triangle PAC */
    double A2 = area (x1, y1, x, y, x3, y3);
    
    /* Calculate area of triangle PAB */
    double A3 = area (x1, y1, x2, y2, x, y);
 
    int isIt = compare_double(A, A1 + A2 + A3);

    /* Check if sum of A1, A2 and A3 is same as A */
    return (isIt);
}

///* Driver program to test above function */
//int main() {
//    
//    /* Let us check whether the point P(10, 15) lies inside the triangle formed by A(0, 0), B(20, 0) and C(10, 30) */
//    if (isPointInTriangle(0, 0, 20, 0, 10, 30, 10, 0) == YES)
//        printf ("Inside\n");
//    else
//        printf ("Not Inside\n");
//    
//    return 0;
//}
