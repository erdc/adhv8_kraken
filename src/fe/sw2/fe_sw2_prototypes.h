#ifndef H_FE_SW2_PROTOTYPES_
#define H_FE_SW2_PROTOTYPES_

/***********************************************************/
/***********************************************************/
/***********************************************************/

int fe_sw2(SMODEL *mod);
void fe_sw2_boundary_resid(SMODEL *mod, DOF_3 *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
void fe_sw2_body_resid(SMODEL *mod, DOF_3 *elem_rhs, int ie, double pertubation, int perturb_node, int perturb_var, int perturb_sign, int DEBUG);
void fe_sw2_get_EEVF(int eev_mode, double eev_coef, double g, double drying_lower_limit, double h_fric, double roughness, double area, SVECT2D grad_u, SVECT2D grad_v, double avg_vmag, int wd_flag, double *ev_st, double *ev_tr);
double fe_sw2_get_roughness(SMODEL *mod, double depth, double velocity, int id, int which);
void fe_sw2_inc(SSUPER_MODEL *sm, int imod);
void fe_sw2_init(SSUPER_MODEL *sm, int imod);
void fe_sw2_load(SSUPER_MODEL *sm, int imod);
void fe_sw2_resid(SSUPER_MODEL *sm, int imod);
void fe_sw2_update(SSUPER_MODEL *sm, int imod);
void fe_sw2_wdflag(SSW_2D *sw2d, SGRID *grid);
int  fe_sw2_solve(SSUPER_MODEL *sm, int imod);

//fe_sw2_wet_dry_integrations.c
void fe_sw2_wd_average_tri(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double *vars, DOF_3 *rhs);
void fe_sw2_wd_average(SVECT *elem_nds, double *depth, SVECT2D *v_wd, double *f_wd, double djac, double *f_wd_avg, SVECT2D *v_wd_avg);
void fe_sw2_wd_integrate_triangle_f(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, DOF_3 *rhs);
void fe_sw2_wd_continuity_temporal_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, DOF_3 *rhs);
void fe_sw2_wd_convection_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, DOF_3 *rhs);
void fe_sw2_wd_gls_convection_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, DOF_3 *rhs);
void fe_sw2_wd_pressure_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, DOF_3 *rhs);
void fe_sw2_wd_bodyForce_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, DOF_3 *rhs);
void fe_sw2_wd_boundaryPressure_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, DOF_3 *rhs);
void fe_sw2_wd_densityPressure_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, DOF_3 *rhs);
void fe_sw2_wd_densityBodyForce_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, DOF_3 *rhs);
void fe_sw2_wd_densityBoundaryPressure_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, DOF_3 *rhs);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
