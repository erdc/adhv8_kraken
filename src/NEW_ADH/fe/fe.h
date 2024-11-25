#ifndef _H_FE_
#define _H_FE_

double integrate_triangle_f(double djac, double c, double *f);
void integrate_triangle_phi_f(double djac, double c, double *f, double *integral);
void integrate_triangle_gradPhi_dot_vbar(SVECT2D *grad_shp, double djac, double c,  SVECT2D vbar, double *integral);

#endif
