#include "adh.h"
static double IS_NEAR_TOL = 1e-9;
bool is_near(double x, double y){
	bool isNear = false;
	if(fabs(x-y)<IS_NEAR_TOL) isNear=true;
	return isNear;
}
