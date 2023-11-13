#ifndef H_FE_DIFFUSIVE_PROTOTYPES_
#define H_FE_DIFFUSIVE_PROTOTYPES_

/***********************************************************/
/***********************************************************/
/***********************************************************/

int fe_diffusive(SMODEL *mod);
int fe_diffusive_solve(SSUPER_MODEL *sm, int imod);
void fe_diffusive_inc(SSUPER_MODEL *sm, int imod);
void fe_diffusive_init(SSUPER_MODEL *sm, int imod);
void fe_diffusive_load(SSUPER_MODEL *sm, int imod);
void fe_diffusive_resid(SSUPER_MODEL *sm, int imod);
void fe_diffusive_update(SSUPER_MODEL *sm, int imod);
void fe_diffusive_boundary_resid(SMODEL *mod, double *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
void fe_diffusive_body_resid(SMODEL *mod, double *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
void fe_diffusive_conv(SVECT *elem_nds, double *depth, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *wd_vel, double *f, double *wd_head, double djac, double *vars, DOF_3 *elem_rhs);
SVECT2D getDiffusiveWaveVelocities(int nnodes, double *depth, double avg_depth, double *z, SVECT2D *grad_phi, double roughness, double manning_units_constant, SVECT2D *vel);
double fe_diffusive_wave_get_tau(double djac, SVECT2D *elem_vel);


/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
