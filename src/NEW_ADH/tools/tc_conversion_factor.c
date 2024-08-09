/* returns the conversion factor to or from seconds - dss */

#include "adh.h"

double tc_conversion_factor(
  int units,			/* units to convert */
  int direct			/* TO or FROM seconds */
)
{
  if(direct == TO)
    {
      switch (units)
	{
	case SECONDS:
	  return 1.;
	case MINUTES:
	  return 60.;
	case HOURS:
	  return 3600.;
	case DAYS:
	  return 86400.;
	case WEEKS:
	  return 604800.;
	case MONTHS:
	  return 2629744.;
	case YEARS:
	  return 31556926.;
	default:
	  return -1.;
	}
    }
  else if(direct == FROM)
    {
      switch (units)
	{
	case SECONDS:
	  return 1;
	case MINUTES:
	  return pow(60.0, -1.0);
	case HOURS:
	  return pow(3600.0, -1.0);
	case DAYS:
	  return pow(86400.0, -1.0);
	case WEEKS:
	  return pow(604800.0, -1.0);
	case MONTHS:
	  return pow(2629744.0, -1.0);
	case YEARS:
	  return pow(31556926.0, -1.0);
	default:
	  return -1;
	}
    }
  else
    return -1;
}
