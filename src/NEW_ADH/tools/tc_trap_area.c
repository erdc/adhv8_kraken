/* calculates the area of a trapezoid */

#include "adh.h"

double tc_trap_area(double lvalue,  /* the left value */
                    double rvalue,  /* the right value */
                    double ltime,   /* the left time */
                    double rtime    /* the right time */
  )
{
  double delta;                 /* the length of the interval */

  /* computes the base length */
  delta = rtime - ltime;

  /* if the base is too small, then return zero */
  if (delta < SMALL)
    return (0.0);
  /* otherwise, return (lvalue+rvalue)*0.5*delta */
  else
    return ((lvalue + rvalue) * 0.5 * delta);
}
