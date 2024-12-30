#include "adh.h"
void set_models(int (*fe_resid[N_RESID_ROUTINES])(SMODEL_SUPER *, double *, int, double, int, int, int, int)){
	fe_resid[SW2] = fe_sw2_body_resid;
	fe_resid[POISSON] = poisson_residual;

}
