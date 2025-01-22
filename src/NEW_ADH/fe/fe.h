#ifndef _H_FE_
#define _H_FE_
//1d segments
void integrate_line_phi_h_h(double djac, double c, double h1, double h2, double *integral);
void integrate_line_phi_h_h_g(double djac, double c, double h1, double h2, double g1, double g2, double *integral); 
double integrate_line_h_f(double djac, double c, double *h, double *f);
//triangles
double integrate_triangle_f(double djac, double c, double *f);
double integrate_triangle_f_f(double djac, double c, double *f);
void integrate_triangle_phi_f(double djac, double c, double *f, double *integral);
void integrate_triangle_phi_f_lump(double djac, double c, double *f, double *integral);
void integrate_triangle_gradPhi_dot_vbar(SVECT2D *grad_shp, double djac, double c,  SVECT2D vbar, double *integral);
void grad2d_phi_f(SVECT2D *grad_phi, double *f, SVECT2D *grad, int ndof);
void grad2d_phi_dot_v(SVECT2D *grad_phi, SVECT2D *v, SVECT2D *grad_x, SVECT2D *grad_y, int ndof);
double integrate_quadrilateral_area(SVECT *nd, double c);
double integrate_quadrilateral_f(SVECT *nd, double c, double *f);
SVECT2D integrate_quadrilateral_gradF(SVECT *nd, double c, double *f);
void integrate_triangle_gradPhi_dot_f_v(SVECT2D *grad_shp, double djac, double c, double *f, SVECT2D *v, double *integral);
void integrate_triangle_dphi_f_g_h(SVECT2D *grad_shp, double djac, double c, double *f, double *g, double *h, double *integral_x, double *integral_y);
void integrate_triangle_dphi_f_f_h(SVECT2D *grad_shp, double djac, double c, double *f, double *h, double *integral_x, double *integral_y);
void integrate_triangle_phi(double djac, double c, double *integral);
void integrate_triangle_dphi_f_f(SVECT2D *grad_shp, double djac, double c, double *f, double *integral_x, double *integral_y);
void integrate_triangle_phi_h_df(double djac, double c, double *h, SVECT2D df, double *integral_x, double *integral_y);
void integrate_triangle_phi_h_g_df(double djac, double c, double *h, double *g, SVECT2D df, double *integral_x, double *integral_y);
void integrate_triangle_phi_f_g(double djac, double c, double *f, double *g, double *integral);
//quads
void extractNodesQuad(SVECT *v, double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4);
void extractVect2DQuad(SVECT2D *v, double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4) ;
void extractFunctionQuad(double *f, double *f1, double *f2, double *f3, double *f4);
void integrate_quadrilateral_gradPhi_dot_f_v(SVECT *nd, double c, double *f, SVECT2D *v, double *integral);
void integrate_quadrilateral_phi_f(SVECT *nd, double c, double *f, double *integral);
void integrate_quadrilateral_phi_f_lump(SVECT *nd, double c, double *f, double *integral);
void integrate_quadrilateral_gradPhi_dot_Dv(SVECT *nd, SQUAD *quad, double c, STENSOR2DAI T, SVECT2D *v, double *integral_x, double *integral_y);
void integrate_quadrilateral_phi_f_g(SVECT *nd, double c, double *f, double *g, double *integral);
void integrate_quadrilateral_gradPhi_dot_f_g_v(SVECT *nd, double c, double *f, double *g, SVECT2D *v, double *integral);
void integrate_quadrilateral_phi(SVECT *nd, double c, double *integral);
void integrate_quadrilateral_dphi_f_f(SVECT *nd, double c, double *f, double *integral_x, double *integral_y);
void integrate_quadrilateral_dphi_f_h_h(SVECT *nd, double c, double *f, double *h, double *integral_x, double *integral_y);
void integrate_quadrilateral_phi_h_g_df(SVECT *nd, double c, double *df, double *h, double *g, double *integral_x, double *integral_y);
void integrate_quadrilateral_gradPhi_dot_vbar(SVECT *nd, double c, SVECT2D vbar, double *integral);
void integrate_quadrilateral_phi_h_df(SVECT *nd, double c, double *df, double *h, double *integral_x, double *integral_y);
double get_quadrilateral_linear_djac2d(double xhat, double yhat, SVECT *nd);
double get_quadrilateral_linear_djac_gradPhi(double xhat, double yhat, SVECT *nd, SVECT *grad_shp);
//other
double integrate_quadZ_f(SVECT *nd, double c, double *f);
//wd integrations
double fe_sw2_wet_dry_factor(SVECT *x, double *h, double djac);
double fe_sw2_wet_dry_wrapper(double *elem_rhs, SVECT *x, double *h, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, int REDISTRIBUTE_FLAG, int DEBUG, double *, void (*fe_sw_func) (SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *elem_rhs));
double fe_sw2_wet_dry_wrapper_1d(double *elem_rhs, SVECT *x, double *h, SVECT2D *nrml, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, int redistribute_flag, int DEBUG, double *vars, void (*fe_sw_func) (SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *elem_rhs));
void fe_sw2_wd_convection_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *elem_rhs);
void fe_sw2_wd_average(SVECT *elem_nds, double *depth, SVECT2D *v_wd, double *f_wd, double djac, double *f_wd_avg, SVECT2D *v_wd_avg);
void fe_sw2_wd_continuity_temporal_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *rhs);
void fe_sw2_wd_boundaryPressure_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *elem_rhs);
void fe_sw2_wd_densityPressure_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *rhs);
void fe_sw2_wd_densityBodyForce_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *rhs);
void fe_sw2_wd_densityBoundaryPressure_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *rhs);
void fe_sw2_wd_integrate_triangle_f(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *rhs);
void fe_sw2_wd_gls_convection_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *rhs);
void fe_sw2_wd_bodyForce_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *rhs);
void fe_sw2_wd_pressure_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *rhs);
#endif
