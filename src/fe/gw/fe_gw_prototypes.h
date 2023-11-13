#ifndef H_FE_GW_PROTOTYPES_
#define H_FE_GW_PROTOTYPES_

#ifdef _ADH_GROUNDWATER
/***********************************************************/
/***********************************************************/
/***********************************************************/
int  fe_gw(SMODEL *mod);
void fe_gw_inc(SSUPER_MODEL *sm, int imod);
void fe_gw_init(SSUPER_MODEL *sm, int imod);
void fe_gw_update(SSUPER_MODEL *sm, int imod);
void fe_gw_load(SSUPER_MODEL *sm, int imod);
void fe_gw_resid(SSUPER_MODEL *sm, int imod);
int fe_gw_solve(SSUPER_MODEL *sm, int isubModel);
void fe_gw_body_resid(SMODEL *mod, double *elem_rhs, int ie, double perturbation,
		      int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
void fe_gw_boundary_resid(SMODEL *mod, double *elem_rhs, int ie, double perturbation,
			  int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
/***********************************************************/
/***********************************************************/
/***********************************************************/
#endif

#endif
