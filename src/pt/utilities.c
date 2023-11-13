#include <stdlib.h>
#include <math.h>

#define PI 3.14159265359

//#define EPSILON  0.0000001
//#define MODULUS(p) (sqrt(p.x*p.x + p.y*p.y + p.z*p.z))
//#define TWOPI 6.283185307179586476925287
//#define RTOD 57.2957795

//===========================================================================
//= find a random number over a user-specified range                        =
//===========================================================================
double find_random_in_range(double min, double max) {
    double scale = rand() / (double) RAND_MAX; /* [0, 1.0] */
    return min + scale * ( max - min );      /* [min, max] */
}

//=========================================================================
//= Multiplicative LCG for generating uniform(0.0, 1.0) random numbers    =
//=   - x_n = 7^5*x_(n-1)mod(2^31 - 1)                                    =
//=   - With x seeded to 1 the 10000th x value should be 1043618065       =
//=   - From R. Jain, "The Art of Computer Systems Performance Analysis," =
//=     John Wiley & Sons, 1991. (Page 443, Figure 26.2)                  =
//=========================================================================
double rand_val(int seed)
{
  const long  a =      16807;  // Multiplier
  const long  m = 2147483647;  // Modulus
  const long  q =     127773;  // m div a
  const long  r =       2836;  // m mod a
  static long x;               // Random int value
  long        x_div_q;         // x divided by q
  long        x_mod_q;         // x modulo q
  long        x_new;           // New x value

  // Set the seed if argument is non-zero and then return zero
  if (seed > 0)
  {
    x = seed;
    return(0.0);
  }

  // RNG using integer arithmetic
  x_div_q = x / q;
  x_mod_q = x % q;
  x_new = (a * x_mod_q) - (r * x_div_q);
  if (x_new > 0)
    x = x_new;
  else
    x = x_new + m;

  // Return a random value between 0.0 and 1.0
  return((double) x / m);
}

//===========================================================================
//=  Function to generate normally distributed random variable using the    =
//=  Box-Muller method                                                      =
//=    - Input: mean and standard deviation                                 =
//=    - Output: Returns with normally distributed random variable          =
//===========================================================================
double find_random_from_normal(double mean, double std_dev) {
    
  double   u, r, theta;           // Variables for Box-Muller method
  double   x;                     // Normal(0, 1) rv
  double   norm_rv;               // The adjusted normal rv

  // Generate u
  u = 0.0;
  while (u == 0.0) u = rand_val(0);

  // Compute r
  r = sqrt(-2.0 * log(u));

  // Generate theta
  theta = 0.0;
  while (theta == 0.0) theta = 2.0 * PI * rand_val(0);

  // Generate x value
  x = r * cos(theta);

  // Adjust x value for specified mean and variance
  norm_rv = (x * std_dev) + mean;

  // Return the normally distributed RV value
  return(norm_rv);
}

//*************************************************************************
//*************************************************************************
//****                                                                 ****
//****                  SUBROUTINE RandomFromvonMises                  ****
//****                                                                 ****
//****              developed by R. Andrew Goodwin, Ph.D.              ****
//****                                                                 ****
//*************************************************************************
//*************************************************************************
//I                              mu,
//I                              kappa,
//I                              Seed,
//O                              Angle_vonMises
double RandomFromvonMises(double mu, double kappa, int Seed) {
    
    // Algorithm from:
    // "Statistical Analysis of Circular Data", Author: N. I. Fisher
    // Random # is between -pi and pi
    
    double Angle_vonMises = 0., RRR;
    
    double a = 1.0 + sqrt(1.0 + 4.0 * (kappa * kappa));
    double b = (a - sqrt(2.0 * a)) / (2.0 * kappa);
    double r = (1.0 + (b * b)) / (2.0 * b);
    
    double dRRR1 = 0., dRRR2 = 0., dRRR3 = 0., f=0., z=0., c=0., PosOrNeg = 0.0;
    while (1) {
        dRRR1 = find_random_in_range(0.0, 1.0);
        dRRR2 = find_random_in_range(0.0, 1.0);
        dRRR3 = find_random_in_range(0.0, 1.0);
        
        z = cos(PI * dRRR1);
        f = (1.0 + r * z) / (r + z);
        c = kappa * (r - f);
        if ((c*(2.0-c)-dRRR2) > 0.0) {
            break;
        } else if  ((log(c/dRRR2)+1.0-c) >= 0.0) {
            break;
        }
    }
    
    if (dRRR3 >= 0.5) {
        PosOrNeg =  1.0;
    } else {
        PosOrNeg = -1.0;
    }
    Angle_vonMises = PosOrNeg * acos(f) + mu;
    Angle_vonMises = fmod(Angle_vonMises,(2.0*PI));
    if (Angle_vonMises > PI) {
        Angle_vonMises = Angle_vonMises - (2.0 * PI);
    }
    
    return Angle_vonMises;
}

//// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//// determine if a point lays within a 2D user-specified polygon
//// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//int InsidePolygon2D(POLYGON *poly,SVECT p) {
//    int i;
//    double angle=0;
//    SVECT p1,p2;
//    
//    for (i=0;i<poly->np;i++) {
//        p1.x = poly->p[i].x - p.x;
//        p1.y = poly->p[i].y - p.y;
//        p2.x = poly->p[(i+1)%poly->n].x - p.x;
//        p2.y = poly->p[(i+1)%poly->n].y - p.y;
//        angle += Angle2D(p1.x,p1.y,p2.x,p2.t);
//    }
//    
//    if (ABS(angle) < PI)
//        return(FALSE);
//    else
//        return(TRUE);
//}
//
//// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
////   Return the angle between two vectors on a plane
////   The angle is from vector 1 to vector 2, positive anticlockwise
////   The result is between -pi -> pi
//// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//double Angle2D(double x1, double y1, double x2, double y2) {
//   double dtheta,theta1,theta2;
//
//   theta1 = atan2(y1,x1);
//   theta2 = atan2(y2,x2);
//   dtheta = theta2 - theta1;
//   while (dtheta > PI)
//      dtheta -= TWOPI;
//   while (dtheta < -PI)
//      dtheta += TWOPI;
//
//   return(dtheta);
//}
//
//// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//// determine if a point lays within a 3D user-specified polygon
//// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//int InsidePolygon3D(POLYGON *poly, SVECT p) {
//    
//    double angle_sum = CalcAngleSum(p,poly->p,poly->np);
//    
//    if (fabs(angle_sum - TWOPI) < 1e-10) {
//        return(TRUE);
//    }
//    return(FALSE);
//}
//
//// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//// The following code snippet returns the angle sum between the test point
//// q and all the vertex pairs. Note that the angle sum is returned in radians.
//// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//double CalcAngleSum(SVECT q, SVECT *p, int n) {
//   int i;
//   double m1,m2;
//   double anglesum=0,costheta;
//   XYZ p1,p2;
//
//   for (i=0; i<n; i++) {
//
//      p1.x = p[i].x - q.x;
//      p1.y = p[i].y - q.y;
//      p1.z = p[i].z - q.z;
//      p2.x = p[(i+1)%n].x - q.x;
//      p2.y = p[(i+1)%n].y - q.y;
//      p2.z = p[(i+1)%n].z - q.z;
//
//      m1 = MODULUS(p1);
//      m2 = MODULUS(p2);
//      if (m1*m2 <= EPSILON)
//         return(TWOPI); /* We are on a node, consider this inside */
//      else
//         costheta = (p1.x*p2.x + p1.y*p2.y + p1.z*p2.z) / (m1*m2);
//
//      anglesum += acos(costheta);
//   }
//   return(anglesum);
//}
