/*! \file  fe_integrals_2d.c This file collections functions responsible for
 *          the 2D FE integrations             */
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \big( c \, f(\widehat{x},\widehat{y})\big) d\widehat{y}\,d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline double integrate_triangle_f(double djac, double c, double *f) {
    return (c * djac * (f[0] + f[1] + f[2]) / 3.0);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \big( c \, \phi_i(\widehat{x},\widehat{y}) \, f(\widehat{x},\widehat{y})\big) d\widehat{y}\,d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triangle_phi_f(double djac, double c, double *f, double *integral) {
    double f_sum = f[0] + f[1] + f[2];
    double t1 = djac * c * one_12;
    integral[0] += t1 * (f_sum + f[0]);
    integral[1] += t1 * (f_sum + f[1]);
    integral[2] += t1 * (f_sum + f[2]);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \big( c\, \mathbf{\nabla} \phi_i(\widehat{x},\widehat{y}) \cdot \mathbf{\overline{v}}  \big)  d\widehat{y}\,d\widehat{x}\f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triangle_gradPhi_dot_vbar(SVECT2D *grad_shp, double djac, double c,  SVECT2D vbar, double *integral) {
    double t1 = c * djac;
    integral[0] += t1 * (vbar.x * grad_shp[0].x + vbar.y * grad_shp[0].y);
    integral[1] += t1 * (vbar.x * grad_shp[1].x + vbar.y * grad_shp[1].y);
    integral[2] += t1 * (vbar.x * grad_shp[2].x + vbar.y * grad_shp[2].y);
}
