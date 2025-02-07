/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  set_models.c This file sets the global arrays of function pointers
 *  to all residual and init routines as created in include/model_codes.h             */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
void set_models( int (*fe_resid[N_RESID_ROUTINES])(SMODEL_SUPER *, double *, int, double, int, int, int, int), int (*fe_init[N_INIT_ROUTINES])(SMODEL_SUPER *)){
	
	//assign all values to residual vector of function pointers
	fe_resid[SW2_BC_DISCHARGE] = fe_sw2_bc_discharge;
	fe_resid[SW2_BC_ELE] = fe_sw2_bc_ele;
	fe_resid[SW2_BC_FLAPD] = fe_sw2_bc_flapd;
	fe_resid[SW2_BC_FLAPU] = fe_sw2_bc_flapu;
	fe_resid[SW2_BC_FLUX] = fe_sw2_bc_flux;
	fe_resid[SW2_BC_H] = fe_sw2_bc_h;
	fe_resid[SW2_BC_HYBRID] = fe_sw2_bc_hybrid;
	fe_resid[SW2_BC_OUTFLOW] = fe_sw2_bc_outflow;
	fe_resid[SW2_BC_SLUICED] = fe_sw2_bc_sluiced;
	fe_resid[SW2_BC_SLUICEU] = fe_sw2_bc_sluiceu;
	fe_resid[SW2_BC_TAILWATER] = fe_sw2_bc_tailwater;
	fe_resid[SW2_BC_VEL] = fe_sw2_bc_vel;
	fe_resid[SW2_BC_VEL_ELE] = fe_sw2_bc_vel_ele;
	fe_resid[SW2_BC_WEIRD] = fe_sw2_bc_weird;
	fe_resid[SW2_BC_WEIRU] = fe_sw2_bc_weiru;
	//assign all values to residual vector of function pointers
	fe_resid[SW2] = fe_sw2_body_resid;
	fe_resid[POISSON] = poisson_residual;
	fe_resid[HEAT] = heat_residual;


	//assign all values to init vector of function pointers
	fe_init[SW2_INIT] = fe_sw2_init;

}
