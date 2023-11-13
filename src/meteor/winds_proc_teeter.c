/* Teeter wind velocity to stress transform (cjt) */
#include "global_header.h"

/* Calculate wind shear components using Teeter's formulation when speeds are input */
void winds_proc_teeter(double gravity, double h, double x_wspeed, double y_wspeed, double *x_wshear, double *y_wshear) {
  double wmag, x_ws, y_ws, hlim, CD, wcfrac;

  wmag = sqrt(pow(x_wspeed, 2.0) + pow(y_wspeed, 2.0) + 1e-15);
  x_ws = x_wspeed;
  y_ws = y_wspeed;
  if(gravity > 10.) {
        h    *= 0.3048;
        wmag *= 0.3048;
        x_ws *= 0.3048;
        y_ws *= 0.3048;
  }
  hlim = MAX(h, .3048);
  CD = pow((0.4 / (16.11 - 0.5 * log(hlim) - 2.48 * log(wmag))), 2.0);
  wcfrac = 0.97;
  if(wmag > 5.0) wcfrac = 2.169 / sqrt(wmag);
  if(hlim > 2.001) wcfrac *= 1.0 / exp(-0.6 * (MIN(hlim, 10.0) - 2.0));
  wcfrac = MIN(wcfrac, 1.0);
  x_ws *= wcfrac * 1.225 * CD * wmag;
  y_ws *= wcfrac * 1.225 * CD * wmag;
  if(gravity > 10.) {
        *x_wshear = x_ws * 0.020868 ;
        *y_wshear = y_ws * 0.020868 ;
  }
  else {
        *x_wshear = x_ws ;
        *y_wshear = y_ws ;
  }
}
