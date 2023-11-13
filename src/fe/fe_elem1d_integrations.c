/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  \brief fe_elem1d_integrations.c This file collects functions that integrate over 
            1D elements. NOTE :: These routines *ADD* to the integral, so zero before call! */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following line segment integration:
 *  \f$ \int_{0}^{1} c \,  f(\widehat{x})  \,  g(\widehat{x})  d\widehat{x}\f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline double integrate_line_h_f(double djac, double c, double *h, double *f) {
    double t1 = djac * one_6 * c;
    return (t1 * ((2*h[0] + h[1])*f[0] + (h[0] + 2*h[1])*f[1]));
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following line segment integration:
 *  \f$ \int_{0}^{1} c \, \phi_i(\widehat{x}) \, f(\widehat{x}) \, (\mathbf{v(\widehat{x},\widehat{y})} \cdot \mathbf{n}) d\widehat{x}\f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_line_phi_h_v_dot_n(double djac, double c, double *h, SVECT2D *v, SVECT2D nrml, double *integral) {
    double t1 = djac * one_12 * c;
    integral[0] += t1 * ( nrml.x * ((3*h[0] + h[1])*v[0].x + (h[0] +   h[1])*v[1].x) + nrml.y * ((3*h[0] + h[1])*v[0].y + (h[0] +   h[1])*v[1].y) );
    integral[1] += t1 * ( nrml.x * ((  h[0] + h[1])*v[0].x + (h[0] + 3*h[1])*v[1].x) + nrml.y * ((  h[0] + h[1])*v[0].y + (h[0] + 3*h[1])*v[1].y) );
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following line segment integration:
 *  \f$ \int_{0}^{1} c \, \phi_i(\widehat{x}) \, f(\widehat{x})^2 d\widehat{x}\f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_line_phi_h_h(double djac, double c, double h1, double h2, double *integral) {
    double t1 = djac * one_12 * c;
    integral[0] += t1 * ( 3*h1*h1 + 2*h1*h2 +   h2*h2 );
    integral[1] += t1 * (   h1*h1 + 2*h1*h2 + 3*h2*h2 );
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following line segment integration:
 *  \f$ \int_{0}^{1} c \, \phi_i(\widehat{x}) \, f(\widehat{x})^2 \, g(\widehat{x})\,  d\widehat{x}\f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_line_phi_h_h_g(double djac, double c, double h1, double h2, double g1, double g2, double *integral) {
    double t1 = djac * one_60 * c;
    integral[0] += t1 * ( 3*(4*g1 + g2)*h1*h1 + 2*(3*g1 + 2*g2)*h1*h2 + (2*g1 + 3*g2)*h2*h2 );
    integral[1] += t1 * ( (3*g1 + 2*g2)*h1*h1 + 2*(2*g1 + 3*g2)*h1*h2 + 3*(g1 + 4*g2)*h2*h2 );
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following line segment integration:
 *  \f$ \int_{0}^{1} c \, \phi_i(\widehat{x}) d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_line_phi(double djac, double c, double *integral) {
    double t1 = one_2 * djac * c;
    integral[0] += t1;
    integral[1] += t1;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following line segment integration:
 *  \f$ \int_{0}^{1} c \, \phi_i(\widehat{x}) * f(\widehat{x}) d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_line_phi_f(double djac, double c, double *f, double *integral) {
    double t1 = one_6 * djac * c;
    integral[0] += t1 * ( 2*f[0] +   f[1] );
    integral[1] += t1 * (   f[0] + 2*f[1] );
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following line segment integration:
 *  \f$ \int_{0}^{1} c \, \phi_i(\widehat{x}) * f(\widehat{x}) * g(\widehat{x}) \, d\widehat{x}\f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_line_phi_f_g(double djac, double c, double *f, double *g, double *integral) {
    double t1 = one_12 * djac * c;
    integral[0] += t1 * ( (3*f[0] + f[1])*g[0] + (f[0] +   f[1])*g[1] );
    integral[1] += t1 * (   (f[0] + f[1])*g[0] + (f[0] + 3*f[1])*g[1] );
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following line segment integration:
 *  \f$ \int_{0}^{1} c \, \phi_i(\widehat{x}) \, f(\widehat{x}) \, g(\widehat{x}) \, h(\widehat{x}) \,d\widehat{x}\f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_line_phi_f_g_h(double djac, double c, double f1, double f2, double g1, double g2, double h1, double h2, double *integral) {
    double t1 = one_60 * djac * c;
    integral[0] += t1 * ( (3*(4*f1 + f2)*g1 + (3*f1 + 2*f2)*g2)*h1 + ((3*f1 + 2*f2)*g1 + (2*f1 + 3*f2)*g2)*h2 );
    integral[1] += t1 * ( ((3*f1 + 2*f2)*g1 + (2*f1 + 3*f2)*g2)*h1 + ((2*f1 + 3*f2)*g1 + 3*(f1 + 4*f2)*g2)*h2 );
}
