/* calculates a value of a time series at a given time */

#include "adh.h"

double tc_eval_series(
        SSERIES series, /* the series */
        int interval,   /* index to the left interval point containing the time */
        double time,    /* the time */
        int index       /* which value */
  ) 
{
  int ileft, iright;       /* the left/right indices of the endpoints of the time interval */
  double lvalue, rvalue;   /* values of the endpoints of the interval containing the time */
  double ltime, rtime;     /* times of the endpoints of the  interval containing the time */
  double theta;            /* parameter telling how far the time is into the time interval */
  double delta;            /* the length of the interval */
  double value;            /* the value of the time series */

  /* checks to see if this is the last entry */
  if (interval == series.size - 1)
    return (series.entry[series.size - 1].value[index]);

  /* finds the interval containing the time */
  ileft = interval;
  iright = interval + 1;

  /* sets the values and the times */
  ltime = series.entry[ileft].time;
  rtime = series.entry[iright].time;
  lvalue = series.entry[ileft].value[index];
  rvalue = series.entry[iright].value[index];

  /* finds how far the time is into the time interval */
  delta = rtime - ltime;
  theta = (rtime - time) / delta;

  /* computes the value */
  value = theta * lvalue + (1.0 - theta) * rvalue;

  /* returns the value */
  return (value);
}
