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