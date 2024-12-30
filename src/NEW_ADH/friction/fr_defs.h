/* ADH Version 2.0.0 6/04 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* constants */
#ifndef PI
#define PI 3.14159265358979
#endif
#define TOL 1e-18
#define BED_CO 1
#define ICE_CO 2
#define TOT_CO 3
#define ERROR -1

/* function prototypes */
double fr_bedshstr_drag_coef(
  double,
  double
);
double fr_equiv_depavg_shstr_coef(
  double,
  double
);
double fr_manningsn_to_rheight(
  double,
  double,
  double
);
double fr_rheight_to_manningsn(
  double,
  double
);
double fr_sav_drag_coef(
  double,
  double,
  double
);
double fr_sav_sed_drag_coef(
  double,
  double
);
double fr_stationary_ice_coef(
  double,
  double,
  double,
  double,
  double,
  int
);
double fr_urv_drag_coef(
  double,
  double,
  double,
  double
);
double fr_urv_sed_drag_coef(
  double,
  double,
  double,
  double
);
