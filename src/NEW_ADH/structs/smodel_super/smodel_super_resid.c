/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smodel_super.c This file collects methods of the SUPER_MODEL structure */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
int smodel_super_resid(SMODEL_SUPER* sm, double *rhs, int ie, double perturbation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG, int (*fe_resid)(SMODEL_SUPER *, double *, int, double, int, int, int, int)){
	//a simple function wrapper
	return fe_resid(sm,rhs,ie, perturbation, perturb_node, perturb_var, perturb_sign, DEBUG);
}
