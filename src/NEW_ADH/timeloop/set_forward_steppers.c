#include  "adh.h"

void set_forward_steppers(int (*forward_stepper[N_TIME_STEPPERS]) (SMODEL_SUPER*)){
	forward_stepper[0] = fe_newton;
}