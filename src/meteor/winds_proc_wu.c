/* Wu wind velocity to stress transform (cjt) */
/* density of air = 1.225 */
#include "global_header.h"

void winds_proc_wu(double gravity, double x_wspeed, double y_wspeed, double *x_wshear, double *y_wshear) {
  double wmag, x_ws, y_ws, CD;

  wmag = sqrt(pow(x_wspeed, 2.0) + pow(y_wspeed, 2.0) + 1e-10);
  x_ws = x_wspeed;
  y_ws = y_wspeed;
  if(gravity > 10.) {
        wmag *= 0.3048;
        x_ws *= 0.3048;
        y_ws *= 0.3048;
  }
  CD = (0.8 + 0.065 * wmag) * .001;
  if(gravity > 10.) {
        *x_wshear = 1.225 * CD * x_ws * wmag * 0.020885;
        *y_wshear = 1.225 * CD * y_ws * wmag * 0.020885;
  }
  else {
        *x_wshear = 1.225 * CD * x_ws * wmag;
        *y_wshear = 1.225 * CD * y_ws * wmag;
  }
}
