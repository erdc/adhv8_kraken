/* Writes the Timeunits card in the output data files for SMS...9.2 or higher */
/* March 4, 2010, JNT */

#include "adh.h"
void tc_timeunits(FILE * fp_out, int outfact)
{
  double conv1, conv2, conv3, conv4;    /* conversion factors for the units from seconds */

  /*  if(time <= tinitial){ */
  conv1 = pow(60.0, -1.0);
  conv2 = pow(3600.0, -1.0);
  conv3 = pow(86400.0, -1.0);
  conv4 = pow(604800.0, -1.0);
  if (outfact == 1.) {
    fprintf(fp_out, "TIMEUNITS SECONDS\n");
  }

  else if (outfact == conv1) {
    fprintf(fp_out, "TIMEUNITS MINUTES\n");
  }

  else if (outfact == conv2) {
    fprintf(fp_out, "TIMEUNITS HOURS\n");
  }

  else if (outfact == conv3) {
    fprintf(fp_out, "TIMEUNITS DAYS\n");
  }

  else if (outfact == conv4) {
    fprintf(fp_out, "TIMEUNITS WEEKS\n");
  }

  else
    tl_error("Timeunits do not match available types.\n");

}
