/*! \file  fe_integrals_2d.c This file collections functions responsible for
 *          the 2D FE integrations             */
#include "adh.h"
inline void extractNodesQuad(SVECT *v, double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4) {
    *x1 = v[0].x; *x2 = v[1].x; *x3 = v[2].x; *x4 = v[3].x;
    *y1 = v[0].y; *y2 = v[1].y; *y3 = v[2].y; *y4 = v[3].y;
}

inline void extractVect2DQuad(SVECT2D *v, double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4) {
    *x1 = v[0].x; *x2 = v[1].x; *x3 = v[2].x; *x4 = v[3].x;
    *y1 = v[0].y; *y2 = v[1].y; *y3 = v[2].y; *y4 = v[3].y;
}

inline void extractFunctionQuad(double *f, double *f1, double *f2, double *f3, double *f4) {
    *f1 = f[0]; *f2 = f[1]; *f3 = f[2]; *f4 = f[3];
}
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
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \big( c \, f(\widehat{x},\widehat{y})^2 \big) d\widehat{y}\,d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline double integrate_triangle_f_f(double djac, double c, double *f) {
    return (c * one_6 * djac * (f[0]*f[0] + f[0]*f[1] + f[1]*f[1] + f[0]*f[2] + f[1]*f[2] + f[2]*f[2]));
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
/*! \brief  Peforms following *lumped triangular integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \big( c \, \phi_i(\widehat{x},\widehat{y}) \, f_i  \big)  d\widehat{y}\,d\widehat{x}\f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triangle_phi_f_lump(double djac, double c, double *f, double *integral) {
    double t1 = djac * c * one_3;
    integral[0] += t1 * f[0];
    integral[1] += t1 * f[1];
    integral[2] += t1 * f[2];
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
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular integration:
 *  \f{eqnarray*}{
 *   integral_x &=&  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \bigg(c \, \phi_i(\widehat{x},\widehat{y}) \, \frac{\partial f(\widehat{x},\widehat{y})}{x} \, h(\widehat{x},\widehat{y})\bigg) d\widehat{y}\,d\widehat{x} \\
 *   integral_y &=&  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \bigg(c \, \phi_i(\widehat{x},\widehat{y}) \, \frac{\partial f(\widehat{x},\widehat{y})}{y} \, h(\widehat{x},\widehat{y})\bigg) d\widehat{y}\,d\widehat{x}
 *  \f}
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triangle_phi_h_g_df(double djac, double c, double *h, double *g, SVECT2D df, double *integral_x, double *integral_y) {
    

    double sum_g = g[0] + g[1] + g[2];
    
    double con, sum_h = h[0] + h[1] + h[2];
    if (integral_x != NULL) {
        con = djac * one_60 * c * df.x;
        integral_x[0] += con * ( (2*(sum_g + 2*g[0])*h[0] + (2*sum_g - g[2])  *h[1] + (2*sum_g - g[1])  *h[2]) );
        integral_x[1] += con * ( ((2*sum_g - g[2])  *h[0] + 2*(sum_g + 2*g[1])*h[1] + (2*sum_g - g[0])  *h[2]) );
        integral_x[2] += con * ( ((2*sum_g - g[1])  *h[0] + (2*sum_g - g[0])  *h[1] + 2*(sum_g + 2*g[2])*h[2]) );
    }
    if (integral_y != NULL) {
        con = djac * one_60 * c * df.y;
        integral_y[0] += con * ( (2*(sum_g + 2*g[0])*h[0] + (2*sum_g - g[2])  *h[1] + (2*sum_g - g[1])  *h[2]) );
        integral_y[1] += con * ( ((2*sum_g - g[2])  *h[0] + 2*(sum_g + 2*g[1])*h[1] + (2*sum_g - g[0])  *h[2]) );
        integral_y[2] += con * ( ((2*sum_g - g[1])  *h[0] + (2*sum_g - g[0])  *h[1] + 2*(sum_g + 2*g[2])*h[2]) );
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates 2D grad(phi) times a function
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] grad a 2D vector array of values for grad(phi) * f
 *  @param[in] grad_phi an array of cartesian space shape function gradients on the element
 *  @param[in] f a double array
 *  @param[in] ndof the number of nodes on the element
 *
 * \f{eqnarray*}{
 *  \bbnabla f = \bbnabla(\sum\limits_{i}{\phidd{i} * f_i}) = \sum\limits_{i}{\bbnabla\phidd{i} * f_i}
 * \f}
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void grad2d_phi_f(SVECT2D *grad_phi, double *f, SVECT2D *grad, int ndof) {
    svect2d_init(grad);
    int idof=0;
    for (idof=0; idof<ndof; idof++) {
        grad->x += f[idof] * grad_phi[idof].x;  grad->y += f[idof] * grad_phi[idof].y;
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates grad(phi) times a 2D vector
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[out] grad_x a 2D vector array of values for grad(phi) * u
 *  @param[out] grad_y a 2D vector array of values for grad(phi) * v
 *  @param[out] grad_z a 2D vector array of values for grad(phi) * w
 *  @param[in] grad_phi an array of cartesian space shape function gradients on the element
 *  @param[in] v a 2D vector array
 *  @param[in] ndof the number of nodes on the element
 *
 * \f{eqnarray*}{
 *  \bbnabla u = \bbnabla(\sum\limits_{i}{\phidd{i} * u_i}) = \sum\limits_{i}{\bbnabla\phidd{i} * u_i}
 *  \bbnabla v = \bbnabla(\sum\limits_{i}{\phidd{i} * v_i}) = \sum\limits_{i}{\bbnabla\phidd{i} * v_i}
 *  \bbnabla w = \bbnabla(\sum\limits_{i}{\phidd{i} * w_i}) = \sum\limits_{i}{\bbnabla\phidd{i} * w_i}
 * \f}
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void grad2d_phi_dot_v(SVECT2D *grad_phi, SVECT2D *v, SVECT2D *grad_x, SVECT2D *grad_y, int ndof) {
    svect2d_init(grad_x);  svect2d_init(grad_y);
    int idof=0;
    for (idof=0; idof<ndof; idof++) {
        grad_x->x += v[idof].x * grad_phi[idof].x;  grad_x->y += v[idof].x * grad_phi[idof].y;
        grad_y->x += v[idof].y * grad_phi[idof].x;  grad_y->y += v[idof].y * grad_phi[idof].y;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular integration:
 *  \f$ \int_{-1}^{1} \int_{-1}^{1} c \, d\widehat{y}\,d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline double integrate_quadrilateral_area(SVECT *nd, double c) {
    return (c * one_2 * (-(nd[1].x - nd[3].x)*nd[0].y + (nd[0].x - nd[2].x)*nd[1].y +
                          (nd[1].x - nd[3].x)*nd[2].y - (nd[0].x - nd[2].x)*nd[3].y) );
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following quadrilateral integration:
 *  \f$ \int_{-1}^{1} \int_{-1}^{1} \big( c \, f(\widehat{x},\widehat{y})\big) d\widehat{y}\,d\widehat{x}  \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline double integrate_quadrilateral_f(SVECT *nd, double c, double *f) {
    
    double f1, f2, f3, f4, x1, x2, x3, x4, y1, y2, y3, y4, t1, t2, t3, t4, t5, t6;
    
    extractFunctionQuad(f,&f1,&f2,&f3,&f4);
    extractNodesQuad(nd,&x1,&x2,&x3,&x4,&y1,&y2,&y3,&y4);
    
    t1 = 2*f1 + 2*f2 +   f3 +   f4, t5 = f2 - f4;
    t2 = 2*f1 +   f2 +   f3 + 2*f4, t6 = f1 - f3;
    t3 =   f1 + 2*f2 + 2*f3 +   f4;
    t4 =   f1 +   f2 + 2*f3 + 2*f4;
    
    return ( one_12 * c * (-(t1*x2 - t5*x3 - t2*x4)*y1
                           +(t1*x1 - t3*x3 - t6*x4)*y2
                           -(t5*x1 - t3*x2 + t4*x4)*y3
                           -(t2*x1 - t6*x2 - t4*x3)*y4) );
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following quadrilateral integration:
 *  \f{eqnarray*}{
 *   f.x &=& \int_{-1}^{1} \int_{-1}^{1} \big( c \, \frac{\partial f(\widehat{x},\widehat{y})}{\partial x} \big) d\widehat{y}\,d\widehat{x} \\
 *   f.y &=& \int_{-1}^{1} \int_{-1}^{1} \big( c \, \frac{\partial f(\widehat{x},\widehat{y})}{\partial y} \big) d\widehat{y}\,d\widehat{x}
 *  \f}
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline SVECT2D integrate_quadrilateral_gradF(SVECT *nd, double c, double *f) {
    SVECT2D result;
    result.x =  -(f[1] - f[3])*nd[0].y + (f[0] - 2*f[2] + f[3])*nd[1].y + 2*(f[1] - f[3])*nd[2].y - (f[0] + f[1] - 2*f[2])*nd[3].y;
    result.y =   (f[1] - f[3])*nd[0].x - (f[0] - 2*f[2] + f[3])*nd[1].x - 2*(f[1] - f[3])*nd[2].x + (f[0] + f[1] - 2*f[2])*nd[3].x;
    result.x *= c * one_24;
    result.y *= c * one_24;
    return result;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \big( c \, \mathbf{\nabla}\phi_i(\widehat{x},\widehat{y}) \cdot (f(\widehat{x},\widehat{y}) \, \mathbf{v(\widehat{x},\widehat{y})})  \big)  d\widehat{y}\,d\widehat{x}\f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triangle_gradPhi_dot_f_v(SVECT2D *grad_shp, double djac, double c, double *f, SVECT2D *v, double *integral) {
    
    double p1, p2, p3, tx, ty, con;
    
    p1 = (2*f[0] +   f[1] +   f[2]);
    p2 = (  f[0] + 2*f[1] +   f[2]);
    p3 = (  f[0] +   f[1] + 2*f[2]);
    
    tx = (p1*v[0].x + p2*v[1].x + p3*v[2].x);
    ty = (p1*v[0].y + p2*v[1].y + p3*v[2].y);
    
    con = djac * one_12 * c;
    integral[0] += con * (grad_shp[0].x * tx + grad_shp[0].y * ty);
    integral[1] += con * (grad_shp[1].x * tx + grad_shp[1].y * ty);
    integral[2] += con * (grad_shp[2].x * tx + grad_shp[2].y * ty);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular integration:
 *  \f{eqnarray*}{
 *   integral_x &=&  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \bigg(c \, \frac{\partial \phi_i(\widehat{x},\widehat{y})}{x} \, f(\widehat{x},\widehat{y}) \, g(\widehat{x},\widehat{y}) \, h(\widehat{x},\widehat{y})\bigg) d\widehat{y}\,d\widehat{x} \\
 *   integral_y &=&  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \bigg(c \, \frac{\partial \phi_i(\widehat{x},\widehat{y})}{y} \, f(\widehat{x},\widehat{y}) \, g(\widehat{x},\widehat{y}) \, h(\widehat{x},\widehat{y})\bigg) d\widehat{y}\,d\widehat{x}
 *  \f}
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triangle_dphi_f_g_h(SVECT2D *grad_shp, double djac, double c, double *f, double *g, double *h, double *integral_x, double *integral_y) {
        
    double con, ff1, ff2, ff3, ff4, ff5, ff6, gg1, gg2, gg3;
    
    ff1 = (3*f[0] +   f[1] + f[2]); ff2 = (2*f[0] + 2*f[1] +   f[2]); ff3 = (2*f[0] + f[1] + 2*f[2]);
    ff4 = (  f[0] + 3*f[1] + f[2]); ff5 = (  f[0] + 2*f[1] + 2*f[2]); ff6 = (  f[0] + f[1] + 3*f[2]);
    gg1 = (2*ff1*g[0] +   ff2*g[1] +   ff3*g[2]);
    gg2 = (  ff2*g[0] + 2*ff4*g[1] +   ff5*g[2]);
    gg3 = (  ff3*g[0] +   ff5*g[1] + 2*ff6*g[2]);
    
    con = djac * one_60 * c * (gg1*h[0] + gg2*h[1] + gg3*h[2]);
    integral_x[0] += con * grad_shp[0].x;
    integral_x[1] += con * grad_shp[1].x;
    integral_x[2] += con * grad_shp[2].x;

    integral_y[0] += con * grad_shp[0].y;
    integral_y[1] += con * grad_shp[1].y;
    integral_y[2] += con * grad_shp[2].y;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following EITHER triangular integration:
 *  \f{eqnarray*}{
 *   integral_x &=&  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \bigg(c \, \frac{\partial \phi_i(\widehat{x},\widehat{y})}{x} \, f(\widehat{x},\widehat{y})^2 \, h(\widehat{x},\widehat{y})\bigg) d\widehat{y}\,d\widehat{x} \\
 *   integral_y &=&  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \bigg(c \, \frac{\partial \phi_i(\widehat{x},\widehat{y})}{y} \, f(\widehat{x},\widehat{y})^2 \, h(\widehat{x},\widehat{y})\bigg) d\widehat{y}\,d\widehat{x}
 *  \f}
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: null pointers decide which derivative direction for output
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triangle_dphi_f_f_h(SVECT2D *grad_shp, double djac, double c, double *f, double *h, double *integral_x, double *integral_y) {
    double p1, p2, p3, con;
    
    double f1_2 = f[0] * f[0], f2_2 = f[1] * f[1], f3_2 = f[2] * f[2];
    
    p1 = (3*f1_2 + 2*f[0]*f[1] +   f2_2 + (2*f[0] +   f[1])*f[2] +   f3_2);
    p2 = (  f1_2 + 2*f[0]*f[1] + 3*f2_2 + (  f[0] + 2*f[1])*f[2] +   f3_2);
    p3 = (  f1_2 +   f[0]*f[1] +   f2_2 + 2*(f[0] +   f[1])*f[2] + 3*f3_2);
    
    con = djac * one_30 * c * (p1*h[0] + p2*h[1] + p3*h[2]);
    if (integral_x != NULL) {
        integral_x[0] += con * grad_shp[0].x;
        integral_x[1] += con * grad_shp[1].x;
        integral_x[2] += con * grad_shp[2].x;
    }
    if (integral_y != NULL) {
        integral_y[0] += con * grad_shp[0].y;
        integral_y[1] += con * grad_shp[1].y;
        integral_y[2] += con * grad_shp[2].y;
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \big( c \, \phi_i(\widehat{x},\widehat{y})\big) d\widehat{y}\,d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triangle_phi(double djac, double c, double *integral) {
    double t1 = c * djac * one_3;
    integral[0] += t1;
    integral[1] += t1;
    integral[2] += t1;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular integration:
 *  \f{eqnarray*}{
 *   integral_x &=&  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \bigg(c \, \frac{\partial \phi_i(\widehat{x},\widehat{y})}{x} \, f(\widehat{x},\widehat{y})^2\bigg) d\widehat{y}\,d\widehat{x} \\
 *   integral_y &=&  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \bigg(c \, \frac{\partial \phi_i(\widehat{x},\widehat{y})}{y} \, f(\widehat{x},\widehat{y})^2\bigg) d\widehat{y}\,d\widehat{x}
 *  \f}
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triangle_dphi_f_f(SVECT2D *grad_shp, double djac, double c, double *f, double *integral_x, double *integral_y) {
    double int_f_f = integrate_triangle_f_f(djac, c, f);
    
    int i;
    if (integral_x != NULL) {
        for (i=0; i<NDONTRI; i++) {
            integral_x[i] += grad_shp[i].x * int_f_f;
        }
    }
    
    if (integral_y != NULL) {
        for (i=0; i<NDONTRI; i++) {
            integral_y[i] += grad_shp[i].y * int_f_f;
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular integration:
 *  \f{eqnarray*}{
 *   integral_x &=&  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \bigg(c \, \phi_i(\widehat{x},\widehat{y}) \, \frac{\partial f(\widehat{x},\widehat{y})}{x} \, h(\widehat{x},\widehat{y})\bigg) d\widehat{y}\,d\widehat{x} \\
 *   integral_y &=&  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \bigg(c \, \phi_i(\widehat{x},\widehat{y}) \, \frac{\partial f(\widehat{x},\widehat{y})}{y} \, h(\widehat{x},\widehat{y})\bigg) d\widehat{y}\,d\widehat{x}
 *  \f}
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triangle_phi_h_df(double djac, double c, double *h, SVECT2D df, double *integral_x, double *integral_y) {
    
    double con, sum_h = h[0] + h[1] + h[2];
    if (integral_x != NULL) {
        con = djac * one_12 * c * df.x;
        integral_x[0] += con * (sum_h + h[0]);
        integral_x[1] += con * (sum_h + h[1]);
        integral_x[2] += con * (sum_h + h[2]);
    }
    if (integral_y != NULL) {
        con = djac * one_12 * c * df.y;
        integral_y[0] += con * (sum_h + h[0]);
        integral_y[1] += con * (sum_h + h[1]);
        integral_y[2] += con * (sum_h + h[2]);
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \big( c \, \phi_i(\widehat{x},\widehat{y}) \, f(\widehat{x},\widehat{y}) \, g(\widehat{x},\widehat{y})\big) d\widehat{y}\,d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triangle_phi_f_g(double djac, double c, double *f, double *g, double *integral) {
    double g_sum = 2 * (g[0] + g[1] + g[2]);
    double t1 = c * djac * one_60;
    integral[0] += t1 * (f[0] * (g_sum + 4*g[0]) + f[1] * (g_sum - g[2]) + f[2] * (g_sum - g[1]));
    integral[1] += t1 * (f[1] * (g_sum + 4*g[1]) + f[2] * (g_sum - g[0]) + f[0] * (g_sum - g[2]));
    integral[2] += t1 * (f[2] * (g_sum + 4*g[2]) + f[1] * (g_sum - g[0]) + f[0] * (g_sum - g[1]));
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following quadrilateral integration:
 *  \f$ \int_{-1}^{1} \int_{-1}^{1} \big( c \, \mathbf{\nabla}\phi_i(\widehat{x},\widehat{y}) \cdot (g(\widehat{x},\widehat{y}) \, \mathbf{v(\widehat{x},\widehat{y})})  \big) d\widehat{y}\,d\widehat{x}   \f$
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: 2D convection type integrals
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_quadrilateral_gradPhi_dot_f_v(SVECT *nd, double c, double *f, SVECT2D *v, double *integral) {
    
    double con;
    double dy1, dy2, ty1, ty2, ty3, ty4, ty5, ty6, ty7, ty8, ty9, ty10, ty11, ty12;
    double dx1, dx2, tx1, tx2, tx3, tx4, tx5, tx6, tx7, tx8, tx9, tx10, tx11, tx12;
    double px1, px2, px3, px4, g1, g2, g3, g4, g5, g6;
    double py1, py2, py3, py4, gx1, gx2, gx3, gx4, gx5, gx6;
    double x1, x2, x3, x4, y1, y2, y3, y4, u1, u2, u3, u4, v1, v2, v3, v4, f1, f2, f3, f4;

    extractVect2DQuad(v,&u1,&u2,&u3,&u4,&v1,&v2,&v3,&v4);
    extractNodesQuad(nd,&x1,&x2,&x3,&x4,&y1,&y2,&y3,&y4);
    extractFunctionQuad(f,&f1,&f2,&f3,&f4);
    
    dy1 = v2 - v4; dx1 = u2 - u4;
    dy2 = v1 - v3; dx2 = u1 - u3;
    py1 = v1 + 4*v2 + v3; px1 = u1 + 4*u2 + u3;
    py2 = v1 + v3 + 4*v4; px2 = u1 + u3 + 4*u4;
    py3 = 4*v1 + v2 + v4; px3 = 4*u1 + u2 + u4;
    py4 = v2 + 4*v3 + v4; px4 = u2 + 4*u3 + u4;
    ty1 = 2*v1 + 2*v2 + v3 + v4;  tx1 = 2*u1 + 2*u2 + u3 + u4;
    ty2 = 2*v1 + v2 + v3 + 2*v4;  tx2 = 2*u1 + u2 + u3 + 2*u4;
    ty3 = v1 + 2*v2 + 2*v3 + v4;  tx3 = u1 + 2*u2 + 2*u3 + u4;
    ty4 = v1 + v2 + 2*v3 + 2*v4;  tx4 = u1 + u2 + 2*u3 + 2*u4;
    ty5 = 6*v1 + 3*v2 + v3 + 2*v4; tx5 = 6*u1 + 3*u2 + u3 + 2*u4;
    ty6 = 3*v1 + 6*v2 + 2*v3 + v4; tx6 = 3*u1 + 6*u2 + 2*u3 + u4;
    ty7 = 2*v1 + 6*v2 + 3*v3 + v4; tx7 = 2*u1 + 6*u2 + 3*u3 + u4;
    ty8 = v1 + 3*v2 + 6*v3 + 2*v4; tx8 = u1 + 3*u2 + 6*u3 + 2*u4;
    ty9 = 6*v1 + 2*v2 + v3 + 3*v4; tx9 = 6*u1 + 2*u2 + u3 + 3*u4;
    ty10 = 3*v1 + v2 + 2*v3 + 6*v4; tx10 = 3*u1 + u2 + 2*u3 + 6*u4;
    ty11 = v1 + 2*v2 + 6*v3 + 3*v4; tx11 = u1 + 2*u2 + 6*u3 + 3*u4;
    ty12 = 2*v1 + v2 + 3*v3 + 6*v4; tx12 = 2*u1 + u2 + 3*u3 + 6*u4;
    
    g1 = ((ty5)*f1 + (ty6)*f2 + (ty3)*f3 + (ty2)*f4);   gx1 = ((tx5)*f1 + (tx6)*f2 + (tx3)*f3 + (tx2)*f4);
    g2 = ((dy1)*f1 + (py1)*f2 + (dy1)*f3 - (py2)*f4);   gx2 = ((dx1)*f1 + (px1)*f2 + (dx1)*f3 - (px2)*f4);
    g3 = ((ty9)*f1 + (ty1)*f2 + (ty4)*f3 + (ty10)*f4);  gx3 = ((tx9)*f1 + (tx1)*f2 + (tx4)*f3 + (tx10)*f4);
    g4 = ((ty1)*f1 + (ty7)*f2 + (ty8)*f3 + (ty4)*f4);   gx4 = ((tx1)*f1 + (tx7)*f2 + (tx8)*f3 + (tx4)*f4);
    g5 = ((py3)*f1 + (dy2)*f2 - (py4)*f3 + (dy2)*f4);   gx5 = ((px3)*f1 + (dx2)*f2 - (px4)*f3 + (dx2)*f4);
    g6 = ((ty2)*f1 + (ty3)*f2 + (ty11)*f3 + (ty12)*f4); gx6 = ((tx2)*f1 + (tx3)*f2 + (tx11)*f3 + (tx12)*f4);
    
    con  = one_72 * c;
    integral[0] += con * (-g1*x2 + g2*x3 + g3*x4 + gx1*y2 - gx2*y3 - gx3*y4);
    integral[1] += con * ( g1*x1 - g4*x3 - g5*x4 - gx1*y1 + gx4*y3 + gx5*y4);
    integral[2] += con * (-g2*x1 + g4*x2 - g6*x4 + gx2*y1 - gx4*y2 + gx6*y4);
    integral[3] += con * (-g3*x1 + g5*x2 + g6*x3 + gx3*y1 - gx5*y2 - gx6*y3);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following quadrilateral integration:
 *  \f$ \int_{-1}^{1} \int_{-1}^{1} \big( c \, \mathbf{\nabla}\phi_i(\widehat{x},\widehat{y}) \cdot (f(\widehat{x},\widehat{y}) \, g(\widehat{x},\widehat{y}) \, \mathbf{v(\widehat{x},\widehat{y})})  \big) d\widehat{y}\,d\widehat{x}   \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_quadrilateral_gradPhi_dot_f_g_v(SVECT *nd, double c, double *f, double *g, SVECT2D *v, double *integral) {
    double con;
    double dy1, dy2, py1, py2, py3, py4, py5, py6, py7, py8;
    double ft1, ft2, ft3, ft4, ft5, ft6, ft7, ft8, ft9, ft10, ft11, ft12, ft13, ft14;
    double ft15, ft16, ft17, ft18, ft19, ft20, ft21, ft22, ft23, ft24, ft25, ft26, ft27, ft28;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
    double hty1, hty2, hty3, hty4, hty5, hty6;
    double htx1, htx2, htx3, htx4, htx5, htx6;
    double x1, x2, x3, x4, y1, y2, y3, y4, u1, u2, u3, u4, v1, v2, v3, v4, f1, f2, f3, f4, g1, g2, g3, g4;

    extractVect2DQuad(v,&u1,&u2,&u3,&u4,&v1,&v2,&v3,&v4);
    extractNodesQuad(nd,&x1,&x2,&x3,&x4,&y1,&y2,&y3,&y4);
    extractFunctionQuad(f,&f1,&f2,&f3,&f4);
    extractFunctionQuad(g,&g1,&g2,&g3,&g4);
    
    // operation count = 34
    dy1 = (f2 - f4);
    dy2 = (f1 - f3);
    py1 = (3*f1 + 6*f2 + f3);
    py2 = (3*f1 + f3 + 6*f4);
    py3 = (f1 + 6*f2 + 3*f3);
    py4 = (f1 + 3*f3 + 6*f4);
    py5 = (6*f1 + 3*f2 + f4);
    py6 = (6*f1 + f2 + 3*f4);
    py7 = (3*f2 + 6*f3 + f4);
    py8 = (f2 + 6*f3 + 3*f4);
    
    // operation count = 28 * 6 = 168
    ft1 =  (12*f1 + 4*f2 + f3 + 3*f4);
    ft2 =  (4*f1 + 4*f2 + f3 + f4);
    ft3 =  (3*f1 + 3*f2 + 2*f3 + 2*f4);
    ft4 =  (9*f1 + 3*f2 + 2*f3 + 6*f4);
    ft5 =  (4*f1 + 12*f2 + 3*f3 + f4);
    ft6 =  (3*f1 + 9*f2 + 6*f3 + 2*f4);
    ft7 =  (2*f1 + 6*f2 + 9*f3 + 3*f4);
    ft8 =  (2*f1 + 2*f2 + 3*f3 + 3*f4);
    ft9 =  (6*f1 + 2*f2 + 3*f3 + 9*f4);
    ft10 = (6*f1 + 27*f2 + 6*f3 + f4);
    ft11 = (6*f1 + f2 + 6*f3 + 27*f4);
    ft12 = (12*f1 + 3*f2 + f3 + 4*f4);
    ft13 = (9*f1 + 6*f2 + 2*f3 + 3*f4);
    ft14 = (3*f1 + 2*f2 + 2*f3 + 3*f4);
    ft15 = (4*f1 + f2 + f3 + 4*f4);
    ft16 = (6*f1 + 9*f2 + 3*f3 + 2*f4);
    ft17 = (2*f1 + 3*f2 + 3*f3 + 2*f4);
    ft18 = (2*f1 + 3*f2 + 9*f3 + 6*f4);
    ft19 = (3*f1 + 2*f2 + 6*f3 + 9*f4);
    ft20 = (4*f1 + f2 + 3*f3 + 12*f4);
    ft21 = (f1 + f2 + 4*f3 + 4*f4);
    ft22 = (3*f1 + f2 + 4*f3 + 12*f4);
    ft23 = (f1 + 3*f2 + 12*f3 + 4*f4);
    ft24 = (f1 + 4*f2 + 4*f3 + f4);
    ft25 = (f1 + 6*f2 + 27*f3 + 6*f4);
    ft26 = (3*f1 + 12*f2 + 4*f3 + f4);
    ft27 = (f1 + 4*f2 + 12*f3 + 3*f4);
    ft28 = (27*f1 + 6*f2 + f3 + 6*f4);
    
    // operation count = 20 * 9 = 180
    t1 =  (3*ft2*g1 + 3*ft5*g2 + ft6*g3 + ft3*g4);
    t2 =  (ft3*g1 + ft6*g2 + ft7*g3 + ft8*g4);
    t3 =  (ft4*g1 + ft3*g2 + ft8*g3 + ft9*g4);
    t4 =  (3*ft1*g1 + 3*ft2*g2 + ft3*g3 + ft4*g4);
    t5 =  (3*ft12*g1 + ft13*g2 + ft14*g3 + 3*ft15*g4);
    t6 =  (ft14*g1 + ft17*g2 + ft18*g3 + ft19*g4);
    t7 =  (ft13*g1 + ft16*g2 + ft17*g3 + ft14*g4);
    t8 =  (ft16*g1 + 3*ft26*g2 + 3*ft24*g3 + ft17*g4);
    t9 =  (ft17*g1 + 3*ft24*g2 + 3*ft27*g3 + ft18*g4);
    t10 = (3*ft15*g1 + ft14*g2 + ft19*g3 + 3*ft20*g4);
    t11 = (ft8*g1 + ft7*g2 + 3*ft23*g3 + 3*ft21*g4);
    t12 = (ft9*g1 + ft8*g2 + 3*ft21*g3 + 3*ft22*g4);
    t13 = (ft28*g1 + py5*g2 + dy2*g3 + py6*g4);
    t14 = (3*dy1*g1 + py1*g2 + dy1*g3 - py2*g4);
    t15 = (dy1*g1 + py3*g2 + 3*dy1*g3 - py4*g4);
    t16 = (py1*g1 + ft10*g2 + py3*g3 + dy1*g4);
    t17 = (py2*g1 - dy1*g2 + py4*g3 + ft11*g4);
    t18 = (py5*g1 + 3*dy2*g2 - py7*g3 + dy2*g4);
    t19 = (dy2*g1 - py7*g2 - ft25*g3 - py8*g4);
    t20 = (py6*g1 + dy2*g2 - py8*g3 + 3*dy2*g4);
    
    // operation count = 6 * 2 * 7 = 84
    hty1 = ( t4*v1 +  t1*v2 +  t2*v3 +  t3*v4);  htx1 = ( t4*u1 +  t1*u2 +  t2*u3 +  t3*u4);
    hty2 = (t14*v1 + t16*v2 + t15*v3 - t17*v4);  htx2 = (t14*u1 + t16*u2 + t15*u3 - t17*u4);
    hty3 = ( t5*v1 +  t7*v2 +  t6*v3 + t10*v4);  htx3 = ( t5*u1 +  t7*u2 +  t6*u3 + t10*u4);
    hty4 = ( t7*v1 +  t8*v2 +  t9*v3 +  t6*v4);  htx4 = ( t7*u1 +  t8*u2 +  t9*u3 +  t6*u4);
    hty5 = (t13*v1 + t18*v2 + t19*v3 + t20*v4);  htx5 = (t13*u1 + t18*u2 + t19*u3 + t20*u4);
    hty6 = ( t3*v1 +  t2*v2 + t11*v3 + t12*v4);  htx6 = ( t3*u1 +  t2*u2 + t11*u3 + t12*u4);
    
    // operation count = 4 * 12 + 2 = 50
    con = (1./720.) * c;
    integral[0] += con * (-hty1*x2 + hty2*x3 + hty3*x4 + htx1*y2 - htx2*y3 - htx3*y4);
    integral[1] += con * ( hty1*x1 - hty4*x3 - hty5*x4 - htx1*y1 + htx4*y3 + htx5*y4);
    integral[2] += con * (-hty2*x1 + hty4*x2 - hty6*x4 + htx2*y1 - htx4*y2 + htx6*y4);
    integral[3] += con * (-hty3*x1 + hty5*x2 + hty6*x3 + htx3*y1 - htx5*y2 - htx6*y3);

}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following quadrilateral integration:
 *  \f$ \int_{-1}^{1} \int_{-1}^{1} \big( c \, \phi_i(\widehat{x},\widehat{y}) \, f(\widehat{x},\widehat{y})\big) d\widehat{y}\,d\widehat{x}   \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_quadrilateral_phi_f(SVECT *nd, double c, double *f, double *integral) {
    
    double f1, f2, f3, f4, x1, x2, x3, x4, y1, y2, y3, y4;

    extractFunctionQuad(f,&f1,&f2,&f3,&f4);
    extractNodesQuad(nd,&x1,&x2,&x3,&x4,&y1,&y2,&y3,&y4);
    
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, d1, d2, c1, c2, c3, c4, c5, con;
    
    t1 = 6*f1 + 3*f2 +   f3 + 2*f4; t2 =  6*f1 + 2*f2 +   f3 + 3*f4;
    t3 = 2*f1 + 2*f2 +   f3 +   f4; t4 =  2*f1 +   f2 +   f3 + 2*f4;
    t5 = 2*f1 + 6*f2 + 3*f3 +   f4; t6 =    f1 + 2*f2 + 2*f3 +   f4;
    t7 =   f1 +   f2 + 2*f3 + 2*f4; t8 =    f1 + 3*f2 + 6*f3 + 2*f4;
    t9 = 3*f1 +   f2 + 2*f3 + 6*f4; t10 = 2*f1 +   f2 + 3*f3 + 6*f4;
    t11 =  f1 + 2*f2 + 6*f3 + 3*f4;
    
    d1 = f2 - f4;
    d2 = f1 - f3;
    c1 = 4*f1 +   f2 +   f4;
    c2 =   f1 + 4*f2 +   f3;
    c3 =   f1 +   f3 + 4*f4;
    c4 =   f2 + 4*f3 +   f4;
    
    con = one_72 * c;
    integral[0] += con * (-(t1*x2 - d1*x3 - t2*x4)*y1 + (t1*x1 - t3*x3 - c1*x4)*y2 - (d1*x1 - t3*x2 +  t4*x4)*y3 - (t2*x1 - c1*x2 -  t4*x3)*y4);
    integral[1] += con * (-(t5*x2 - c2*x3 - t3*x4)*y1 + (t5*x1 - t5*x3 - d2*x4)*y2 - (c2*x1 - t5*x2 +  t6*x4)*y3 - (t3*x1 - d2*x2 -  t6*x3)*y4);
    integral[2] += con * (-(t6*x2 - d1*x3 - t7*x4)*y1 + (t6*x1 - t8*x3 + c4*x4)*y2 - (d1*x1 - t8*x2 + t11*x4)*y3 - (t7*x1 + c4*x2 - t11*x3)*y4);
    integral[3] += con * (-(t4*x2 + c3*x3 - t9*x4)*y1 + (t4*x1 - t7*x3 - d2*x4)*y2 + (c3*x1 + t7*x2 - t10*x4)*y3 - (t9*x1 - d2*x2 - t10*x3)*y4);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following quadrilateral integration:
 *  \f$ \int_{-1}^{1} \int_{-1}^{1} \big( c \, \phi_i(\widehat{x},\widehat{y})\big) d\widehat{y}\,d\widehat{x}  \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_quadrilateral_phi(SVECT *nd, double c, double *integral) {
    
   double x1, x2, x3, x4, y1, y2, y3, y4;

    extractNodesQuad(nd,&x1,&x2,&x3,&x4,&y1,&y2,&y3,&y4);
    
    double con = one_12 * c;
    integral[0] += con * (-2*(x2 - x4)*y1 + (2*x1 - x3 - x4)*y2 + (x2 - x4)*y3 - (2*x1 - x2 - x3)*y4);
    integral[1] += con * (-(2*x2 - x3 - x4)*y1 + 2*(x1 - x3)*y2 - (x1 - 2*x2 + x4)*y3 - (x1 - x3)*y4);
    integral[2] += con * (-(x2 - x4)*y1 + (x1 - 2*x3 + x4)*y2 + 2*(x2 - x4)*y3 - (x1 + x2 - 2*x3)*y4);
    integral[3] += con * (-(x2 + x3 - 2*x4)*y1 + (x1 - x3)*y2 + (x1 + x2 - 2*x4)*y3 - 2*(x1 - x3)*y4);
    
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following *lumped quadrilateral integration:
 *  \f$ \int_{-1}^{1} \int_{-1}^{1} \big( c \, \phi_i(\widehat{x},\widehat{y}) \, f_i  \big)  d\widehat{y}\,d\widehat{x}\f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_quadrilateral_phi_f_lump(SVECT *nd, double c, double *f, double *integral) {
    integrate_quadrilateral_phi(nd, c, integral);
    integral[0] += integral[0] * f[0];
    integral[1] += integral[1] * f[1];
    integral[2] += integral[2] * f[2];
    integral[3] += integral[3] * f[3];
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// returns the 2D projected quadrilateral Jacobian
inline double get_quadrilateral_linear_djac2d(double xhat, double yhat, SVECT *nd) {
    double x12, x23, x34, x14, x24, x13, tx, ty, tc, term;
    x12 = (nd[0].x - nd[1].x);
    x23 = (nd[1].x - nd[2].x);
    x34 = (nd[2].x - nd[3].x);
    x14 = (nd[0].x - nd[3].x);
    x24 = (nd[1].x - nd[3].x);
    x13 = (nd[0].x - nd[2].x);
    tx = ( x34*nd[0].y - x34*nd[1].y - x12*nd[2].y + x12*nd[3].y);
    ty = ( x23*nd[0].y - x14*nd[1].y + x14*nd[2].y - x23*nd[3].y);
    tc = (-x24*nd[0].y + x13*nd[1].y + x24*nd[2].y - x13*nd[3].y);
    term = tx*xhat + ty*yhat + tc;
    return (one_8 * term);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// returns the 2D projected quadrilateral Jacobian and cartesian space shape function gradients
inline double get_quadrilateral_linear_djac_gradPhi(double xhat, double yhat, SVECT *nd, SVECT *grad_shp) {
    
    double djac2d = get_quadrilateral_linear_djac2d(xhat, yhat, nd);
    double con = (1./(djac2d * 8.));
    grad_shp[ 0 ].x = con * (       -xhat*nd[2].y + (xhat - 1)*nd[3].y - (nd[1].y - nd[2].y)*yhat + nd[1].y );
    grad_shp[ 1 ].x = con * (  (xhat + 1)*nd[2].y -       xhat*nd[3].y + (nd[0].y - nd[3].y)*yhat - nd[0].y );
    grad_shp[ 2 ].x = con * (        xhat*nd[0].y - (xhat + 1)*nd[1].y - (nd[0].y - nd[3].y)*yhat + nd[3].y );
    grad_shp[ 3 ].x = con * ( -(xhat - 1)*nd[0].y +       xhat*nd[1].y + (nd[1].y - nd[2].y)*yhat - nd[2].y );
    
    grad_shp[ 0 ].y = con * (  (nd[2].x - nd[3].x)*xhat + (nd[1].x - nd[2].x)*yhat - nd[1].x + nd[3].x );
    grad_shp[ 1 ].y = con * ( -(nd[2].x - nd[3].x)*xhat - (nd[0].x - nd[3].x)*yhat + nd[0].x - nd[2].x );
    grad_shp[ 2 ].y = con * ( -(nd[0].x - nd[1].x)*xhat + (nd[0].x - nd[3].x)*yhat + nd[1].x - nd[3].x );
    grad_shp[ 3 ].y = con * (  (nd[0].x - nd[1].x)*xhat - (nd[1].x - nd[2].x)*yhat - nd[0].x + nd[2].x );
    
    return djac2d;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following quadrilateral integration:
 *  \f{eqnarray*}{
 *   integral_x &=& \int_{-1}^{1} \int_{-1}^{1} \bigg(c \, \mathbf{\nabla}\phi_i(\widehat{x},\widehat{y}) \cdot (\mathbf{D(\mathbf{\nabla}\mathbf{u}(\widehat{x},\widehat{y}))}) \bigg) d\widehat{y}\,d\widehat{x} \\
 *   integral_y &=& \int_{-1}^{1} \int_{-1}^{1} \bigg(c \, \mathbf{\nabla}\phi_i(\widehat{x},\widehat{y}) \cdot (\mathbf{D(\mathbf{\nabla}\mathbf{v}(\widehat{x},\widehat{y}))}) \bigg) d\widehat{y}\,d\widehat{x}
 *  \f}
 *  \n where \n
 *  \f$ D_x = \sigma_{xx} \, \frac{\partial u}{\partial x} + \sigma_{xy} \, \frac{\partial u}{\partial y} \f$ \n
 *  \f$ D_y = \sigma_{xy} \, \frac{\partial v}{\partial x} + \sigma_{yy} \, \frac{\partial v}{\partial y} \f$
 *  \author  Corey Trahan, Ph.D.
*  \note CJT \:: Used for vector (momentum) A/D equations w/ potential anisotropic diffusion
 *  \note CJT \:: Use quadrature here, since dphi * dphi integrals will not solve analytically
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_quadrilateral_gradPhi_dot_Dv(SVECT *nd, SQUAD *quad, double c, STENSOR2DAI T, SVECT2D *v, double *integral_x, double *integral_y) {
    
    int i, iqp, quad_order = 3; // quadrature order
    SVECT2D qp_grad_u, qp_grad_v;
    SQUAD_PT *qp = NULL;
    
    double resultX[NDONQUAD] = {0., 0., 0., 0.};
    double resultY[NDONQUAD] = {0., 0., 0., 0.};
    
    for (iqp=0; iqp<quad[quad_order].n; iqp++) {
        qp = &(quad[quad_order].pt[iqp]);
        
        // evaluate quadrilaterial djac and cartesian shape function gradients at quadrature points
        qp->djac = get_quadrilateral_linear_djac_gradPhi(qp->xhat, qp->yhat, nd, qp->grad_shp);
        
        // evaluate grad(f) at quadrature points
        svect2d_init(&qp_grad_u);  svect2d_init(&qp_grad_v);
        for (i=0; i<NDONQUAD; i++) {
            qp_grad_u.x += v[i].x * qp->grad_shp[i].x;
            qp_grad_u.y += v[i].x * qp->grad_shp[i].y;
            qp_grad_v.x += v[i].y * qp->grad_shp[i].x;
            qp_grad_v.y += v[i].y * qp->grad_shp[i].y;
        }
        
        double t1 = c * qp->djac * qp->w;
        for (i=0; i<NDONQUAD; i++) {
            resultX[i] += t1 * qp->grad_shp[i].x * (T.xx * qp_grad_u.x + T.xy * qp_grad_u.y);
            resultY[i] += t1 * qp->grad_shp[i].y * (T.yx * qp_grad_v.x + T.yy * qp_grad_v.y);
        }
    }
    
    for (i=0; i<NDONQUAD; i++) {
        integral_x[i] += resultX[i];
        integral_y[i] += resultY[i];
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following quadrilateral integration:
 *  \f$ \int_{-1}^{1} \int_{-1}^{1} \big( c \, \phi_i(\widehat{x},\widehat{y}) \, f(\widehat{x},\widehat{y})\, g(\widehat{x},\widehat{y})\big) d\widehat{y}\,d\widehat{x}   \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_quadrilateral_phi_f_g(SVECT *nd, double c, double *f, double *g, double *integral) {
    
    double f1, f2, f3, f4, g1, g2, g3, g4, x1, x2, x3, x4, y1, y2, y3, y4;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
    double t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28;
    double d1, d2, d3, d4, d5, d6, d7, d8, p1, p2;
    double gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg8, gg9, gg10, gg11, gg12, gg13, gg14, gg15, gg16, gg17, gg18, gg19, gg20;
    double con;
    
    extractFunctionQuad(f,&f1,&f2,&f3,&f4);
    extractFunctionQuad(g,&g1,&g2,&g3,&g4);
    extractNodesQuad(nd,&x1,&x2,&x3,&x4,&y1,&y2,&y3,&y4);
    
    t1 =  (12*f1 +  4*f2 +    f3 +  3*f4); t2 = ( 4*f1 +  4*f2 +   f3 +    f4); t3 =   (3*f1 + 3*f2 +  2*f3 +  2*f4);
    t4 =  ( 9*f1 +  3*f2 +  2*f3 +  6*f4); t5 = (12*f1 +  3*f2 +   f3 +  4*f4); t6 =   (3*f1 + 2*f2 +  2*f3 +  3*f4);
    t7 =  ( 9*f1 +  6*f2 +  2*f3 +  3*f4); t8 = ( 4*f1 +    f2 +   f3 +  4*f4); t9 =   (6*f1 + 9*f2 +  3*f3 +  2*f4);
    t10 = ( 2*f1 +  2*f2 +  3*f3 +  3*f4); t11 = (2*f1 +  3*f2 + 3*f3 +  2*f4); t12 = (27*f1 + 6*f2 +    f3 +  6*f4);
    t13 = ( 4*f1 + 12*f2 +  3*f3 +    f4); t14 = (6*f1 + 27*f2 + 6*f3 +    f4); t15 =  (3*f1 + 9*f2 +  6*f3 +  2*f4);
    t16 = ( 3*f1 + 12*f2 +  4*f3 +    f4); t17 = (  f1 +  4*f2 + 4*f3 +    f4); t18 =  (2*f1 + 6*f2 +  9*f3 +  3*f4);
    t19 = ( 2*f1 +  3*f2 +  9*f3 +  6*f4); t20 = (3*f1 +  2*f2 + 6*f3 +  9*f4); t21 =    (f1 + 4*f2 + 12*f3 +  3*f4);
    t22 = (   f1 +  3*f2 + 12*f3 +  4*f4); t23 = (  f1 +    f2 + 4*f3 +  4*f4); t24 =    (f1 + 6*f2 + 27*f3 +  6*f4);
    t25 = ( 6*f1 +  2*f2 +  3*f3 +  9*f4); t26 = (6*f1 +    f2 + 6*f3 + 27*f4); t27 =  (3*f1 +   f2 +  4*f3 + 12*f4);
    t28 = ( 4*f1 +    f2 +  3*f3 + 12*f4);
    
    d1 = (3*f1 + 6*f2 +   f3); d2 = (3*f1 +   f3 + 6*f4); d3 = (6*f1 + 3*f2 + f4); d4 = (6*f1 +   f2 + 3*f4);
    d5 = (  f1 + 6*f2 + 3*f3); d6 = (  f1 + 3*f3 + 6*f4); d7 = (3*f2 + 6*f3 + f4); d8 = (  f2 + 6*f3 + 3*f4);
    
    p1 = (f2 - f4); p2 = (f1 - f3);
    
    gg1 = (3*t1*g1 +  3*t2*g2 +    t3*g3 +    t4*g4); gg2 = (3*p1*g1 +    d1*g2 +    p1*g3 -    d2*g4); gg3 = (3*t5*g1 +    t7*g2 +    t6*g3 +  3*t8*g4);
    gg4 = (  t7*g1 +    t9*g2 +   t11*g3 +    t6*g4); gg5 = ( t12*g1 +    d3*g2 +    p2*g3 +    d4*g4); gg6 = (  t4*g1 +    t3*g2 +   t10*g3 +   t25*g4);
    gg7 = (3*t2*g1 + 3*t13*g2 +   t15*g3 +    t3*g4); gg8 = (  d1*g1 +   t14*g2 +    d5*g3 +    p1*g4); gg9 = (  t9*g1 + 3*t16*g2 + 3*t17*g3 +   t11*g4);
    gg10 =  (d3*g1 +  3*p2*g2 -    d7*g3 +    p2*g4); gg11 = ( t3*g1 +   t15*g2 +   t18*g3 +   t10*g4); gg12 = ( p1*g1 +    d5*g2 +  3*p1*g3 -    d6*g4);
    gg13 = ( t6*g1 +   t11*g2 +   t19*g3 +   t20*g4); gg14 = (t11*g1 + 3*t17*g2 + 3*t21*g3 +   t19*g4); gg15 = ( p2*g1 -    d7*g2 -   t24*g3 -    d8*g4);
    gg16 = (t10*g1 +   t18*g2 + 3*t22*g3 + 3*t23*g4); gg17 = ( d2*g1 -    p1*g2 +    d6*g3 +   t26*g4); gg18 = (3*t8*g1 +    t6*g2 +   t20*g3 + 3*t28*g4);
    gg19 = ( d4*g1 +    p2*g2 -    d8*g3 +  3*p2*g4); gg20 = (t25*g1 +   t10*g2 + 3*t23*g3 + 3*t27*g4);
    
    con = c * one_720;
    integral[0] += con * (-( gg1*x2 -  gg2*x3 -  gg3*x4)*y1 + ( gg1*x1 -  gg4*x3 -  gg5*x4)*y2 - ( gg2*x1 -  gg4*x2 +  gg6*x4)*y3 - ( gg3*x1 -  gg5*x2 -  gg6*x3)*y4);
    integral[1] += con * (-( gg7*x2 -  gg8*x3 -  gg4*x4)*y1 + ( gg7*x1 -  gg9*x3 - gg10*x4)*y2 - ( gg8*x1 -  gg9*x2 + gg11*x4)*y3 - ( gg4*x1 - gg10*x2 - gg11*x3)*y4);
    integral[2] += con * (-(gg11*x2 - gg12*x3 - gg13*x4)*y1 + (gg11*x1 - gg14*x3 - gg15*x4)*y2 - (gg12*x1 - gg14*x2 + gg16*x4)*y3 - (gg13*x1 - gg15*x2 - gg16*x3)*y4);
    integral[3] += con * (-( gg6*x2 + gg17*x3 - gg18*x4)*y1 + ( gg6*x1 - gg13*x3 - gg19*x4)*y2 + (gg17*x1 + gg13*x2 - gg20*x4)*y3 - (gg18*x1 - gg19*x2 - gg20*x3)*y4);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following quadrilateral integration:
 *  \f{eqnarray*}{
 *   integral_x &=& \int_{-1}^{1} \int_{-1}^{1} \bigg(c \, \frac{\partial \phi_i(\widehat{x},\widehat{y})}{x} \, f(\widehat{x},\widehat{y})^2\bigg) d\widehat{y}\,d\widehat{x} \\
 *   integral_y &=& \int_{-1}^{1} \int_{-1}^{1} \bigg(c \, \frac{\partial \phi_i(\widehat{x},\widehat{y})}{y} \, f(\widehat{x},\widehat{y})^2\bigg) d\widehat{y}\,d\widehat{x}
 *  \f}
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_quadrilateral_dphi_f_f(SVECT *nd, double c, double *f, double *integral_x, double *integral_y) {
    double x1, x2, x3, x4, y1, y2, y3, y4;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, g1, g2, g3, g4, g5, g6, f1, f2, f3, f4;
    double con;
    
    extractFunctionQuad(f,&f1,&f2,&f3,&f4);
    extractNodesQuad(nd,&x1,&x2,&x3,&x4,&y1,&y2,&y3,&y4);
    
    double f1_2 = f1 * f1, f2_2 = f2 * f2, f3_2 = f3 * f3, f4_2 = f4 * f4;
    
    con = c * one_36;
    
    t1 = f1 + 2*f2;  t2 = f1 + f3; t3 = f1 + 3*f2; t4 = f1 - f3; t9 = 2*f1 + f2 + 3*f3;
    t5 = 2*f1 + f2 + f3;  t6 = f1 + f2; t7 = 3*f1 + f2 + 2*f3; t8 = f1 + f2 + 2*f3;
    
    g1 = f1*f2 + 2*f2_2 + f2*f3 - t2*f4 - 2*f4_2;
    g2 = 2*f1_2 + f1*f2 - f2*f3 - 2*f3_2 + t4*f4;
    g3 = 3*f1_2 + 3*f1*f2 + 3*f2_2 + t1*f3 + f3_2 + t5*f4 + f4_2;
    g4 = 3*f1_2 + 2*f1*f2 + f2_2 + t6*f3 + f3_2 + t7*f4 + 3*f4_2;
    g5 = f1_2 + 2*f1*f2 + 3*f2_2 + t3*f3 + 3*f3_2 + t8*f4 + f4_2;
    g6 = f1_2 + f1*f2 + f2_2 + t1*f3 + 3*f3_2 + t9*f4 + 3*f4_2;
    
    integral_x[0] += con * ( g3*y2 - g1*y3 - g4*y4);
    integral_x[1] += con * (-g3*y1 + g5*y3 + g2*y4);
    integral_x[2] += con * ( g1*y1 - g5*y2 + g6*y4);
    integral_x[3] += con * ( g4*y1 - g2*y2 - g6*y3);
    
    integral_y[0] += con * (-g3*x2 + g1*x3 + g4*x4);
    integral_y[1] += con * ( g3*x1 - g5*x3 - g2*x4);
    integral_y[2] += con * (-g1*x1 + g5*x2 - g6*x4);
    integral_y[3] += con * (-g4*x1 + g2*x2 + g6*x3);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following quadrilateral integration:
 *  \f{eqnarray*}{
 *   integral_x &=& \int_{-1}^{1} \int_{-1}^{1} \bigg(c \, \frac{\partial \phi_i(\widehat{x},\widehat{y})}{x} \, f(\widehat{x},\widehat{y})\, h(\widehat{x},\widehat{y})^2\bigg) d\widehat{y}\,d\widehat{x} \\
 *   integral_y &=& \int_{-1}^{1} \int_{-1}^{1} \bigg(c \, \frac{\partial \phi_i(\widehat{x},\widehat{y})}{y} \, f(\widehat{x},\widehat{y})\, h(\widehat{x},\widehat{y})^2\bigg) d\widehat{y}\,d\widehat{x}
 *  \f}
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_quadrilateral_dphi_f_h_h(SVECT *nd, double c, double *f, double *h, double *integral_x, double *integral_y) {

    double x1, x2, x3, x4, y1, y2, y3, y4, f1, f2, f3, f4, h1, h2, h3, h4;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
    double t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28;
    double d1, d2, d3, d4, d5, d6, d7, d8, p1, p2, gg1, gg2, gg3, gg4, gg5, gg6;
    double con;

    extractNodesQuad(nd,&x1,&x2,&x3,&x4,&y1,&y2,&y3,&y4);
    extractFunctionQuad(f,&f1,&f2,&f3,&f4);
    extractFunctionQuad(h,&h1,&h2,&h3,&h4);
    
    double h1_2 = h1 * h1, h2_2 = h2 * h2, h3_2 = h3 * h3, h4_2 = h4 * h4;
    
    t1 =  (12*f1 +  4*f2 +    f3 +  3*f4); t2 = ( 4*f1 +  4*f2 +   f3 +    f4); t3 =   (3*f1 + 3*f2 +  2*f3 +  2*f4);
    t4 =  ( 9*f1 +  3*f2 +  2*f3 +  6*f4); t5 = (12*f1 +  3*f2 +   f3 +  4*f4); t6 =   (3*f1 + 2*f2 +  2*f3 +  3*f4);
    t7 =  ( 9*f1 +  6*f2 +  2*f3 +  3*f4); t8 = ( 4*f1 +    f2 +   f3 +  4*f4); t9 =   (6*f1 + 9*f2 +  3*f3 +  2*f4);
    t10 = ( 2*f1 +  2*f2 +  3*f3 +  3*f4); t11 = (2*f1 +  3*f2 + 3*f3 +  2*f4); t12 = (27*f1 + 6*f2 +    f3 +  6*f4);
    t13 = ( 4*f1 + 12*f2 +  3*f3 +    f4); t14 = (6*f1 + 27*f2 + 6*f3 +    f4); t15 =  (3*f1 + 9*f2 +  6*f3 +  2*f4);
    t16 = ( 3*f1 + 12*f2 +  4*f3 +    f4); t17 = (  f1 +  4*f2 + 4*f3 +    f4); t18 =  (2*f1 + 6*f2 +  9*f3 +  3*f4);
    t19 = ( 2*f1 +  3*f2 +  9*f3 +  6*f4); t20 = (3*f1 +  2*f2 + 6*f3 +  9*f4); t21 =    (f1 + 4*f2 + 12*f3 +  3*f4);
    t22 = (   f1 +  3*f2 + 12*f3 +  4*f4); t23 = (  f1 +    f2 + 4*f3 +  4*f4); t24 =    (f1 + 6*f2 + 27*f3 +  6*f4);
    t25 = ( 6*f1 +  2*f2 +  3*f3 +  9*f4); t26 = (6*f1 +    f2 + 6*f3 + 27*f4); t27 =  (3*f1 +   f2 +  4*f3 + 12*f4);
    t28 = ( 4*f1 +    f2 +  3*f3 + 12*f4);
    
    d1 = (3*f1 + 6*f2 +   f3); d2 = (3*f1 +   f3 + 6*f4); d3 = (6*f1 + 3*f2 + f4); d4 = (6*f1 +   f2 + 3*f4);
    d5 = (  f1 + 6*f2 + 3*f3); d6 = (  f1 + 3*f3 + 6*f4); d7 = (3*f2 + 6*f3 + f4); d8 = (  f2 + 6*f3 + 3*f4);
    
    p1 = (f2 - f4); p2 = (f1 - f3);
    
    gg1 = (3*t1*h1_2 + 6*t2*h1*h2 + 3*t13*h2_2 + t18*h3_2 + t25*h4_2 + 2*(t3*h1 + t15*h2)*h3 + 2*(t4*h1 + t3*h2 + t10*h3)*h4);
    gg2 = (3*p1*h1_2 + 2*d1*h1*h2 + t14*h2_2 + 3*p1*h3_2 - t26*h4_2 + 2*(p1*h1 + d5*h2)*h3 - 2*(d2*h1 - p1*h2 + d6*h3)*h4);
    gg3 = (3*t5*h1_2 + 2*t7*h1*h2 + t9*h2_2 + t19*h3_2 + 3*t28*h4_2 + 2*(t6*h1 + t11*h2)*h3 + 2*(3*t8*h1 + t6*h2 + t20*h3)*h4);
    gg4 = (t7*h1_2 + 2*t9*h1*h2 + 3*t16*h2_2 + 3*t21*h3_2 + t20*h4_2 + 2*(t11*h1 + 3*t17*h2)*h3 + 2*(t6*h1 + t11*h2 + t19*h3)*h4);
    gg5 = (t12*h1_2 + 2*d3*h1*h2 + 3*p2*h2_2 - t24*h3_2 + 3*p2*h4_2 + 2*(p2*h1 - d7*h2)*h3 + 2*(d4*h1 + p2*h2 - d8*h3)*h4);
    gg6 = (t4*h1_2 + 2*t3*h1*h2 + t15*h2_2 + 3*t22*h3_2 + 3*t27*h4_2 + 2*(t10*h1 + t18*h2)*h3 + 2*(t25*h1 + t10*h2 + 3*t23*h3)*h4);
    
    
    con = one_720 * c;
    integral_x[0] += con * ( gg1*y2 - gg2*y3 - gg3*y4);
    integral_x[1] += con * (-gg1*y1 + gg4*y3 + gg5*y4);
    integral_x[2] += con * ( gg2*y1 - gg4*y2 + gg6*y4);
    integral_x[3] += con * ( gg3*y1 - gg5*y2 - gg6*y3);
    
    
    integral_y[0] += con * (-gg1*x2 + gg2*x3 + gg3*x4);
    integral_y[1] += con * ( gg1*x1 - gg4*x3 - gg5*x4);
    integral_y[2] += con * (-gg2*x1 + gg4*x2 - gg6*x4);
    integral_y[3] += con * (-gg3*x1 + gg5*x2 + gg6*x3);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following quadrilateral integration:
 *  \f{eqnarray*}{
 *   integral_x &=& \int_{-1}^{1} \int_{-1}^{1} \bigg(c \, \phi_i(\widehat{x},\widehat{y}) \, \frac{\partial \f(\widehat{x},\widehat{y})}{x} \, g(\widehat{x},\widehat{y})\, h(\widehat{x},\widehat{y}) \bigg) d\widehat{y}\,d\widehat{x} \\
 *   integral_y &=& \int_{-1}^{1} \int_{-1}^{1} \bigg(c \, \phi_i(\widehat{x},\widehat{y}) \, \frac{\partial \f(\widehat{x},\widehat{y})}{y} \, g(\widehat{x},\widehat{y})\, h(\widehat{x},\widehat{y}) \bigg) d\widehat{y}\,d\widehat{x}
 *  \f}
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: SW 2D pressure density body integral
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_quadrilateral_phi_h_g_df(SVECT *nd, double c, double *df, double *h, double *g, double *integral_x, double *integral_y) {
    
    double x1, x2, x3, x4, y1, y2, y3, y4, f1, f2, f3, f4, h1, h2, h3, h4, g1, g2, g3, g4;
    double con = 0.;
    
    extractNodesQuad(nd,&x1,&x2,&x3,&x4,&y1,&y2,&y3,&y4);
    extractFunctionQuad(df,&f1,&f2,&f3,&f4);
    extractFunctionQuad(h,&h1,&h2,&h3,&h4);
    extractFunctionQuad(g,&g1,&g2,&g3,&g4);
    
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21;
    double t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42;
    double p1, p2, p3, p4, p5, p6, p7, p8;
    double d1, d2, d3, d4, d5, d6, d7, d8;
    
    
    t1 = (f2 - f4); t2 = (f1 - f3); t3 = (4*f2 - f3 - 3*f4); t4 = (3*f2 + f3 - 4*f4); t5 = (2*f2 - f3 - f4); t6 = (3*f2 - f3 - 2*f4);
    t7 = (2*f2 + f3 - 3*f4); t8 = (f2 + f3 - 2*f4); t9 = (4*f1 - f3 - 3*f4); t10 = (2*f1 - f3 - f4); t11 = (3*f1 - 2*f3 - f4);
    t12 = (3*f1 - f3 - 2*f4); t13 = (4*f1 - 3*f3 - f4); t14 = (2*f1 - 3*f3 + f4); t15 = (f1 - 2*f2 + f4); t16 = (f1 + f2 - 2*f4);
    t17 = (2*f1 - 3*f2 + f4); t18 = (f1 - 3*f2 + 2*f4); t19 = (f1 + 2*f2 - 3*f4);t20 = (2*f1 + f2 - 3*f4);t21 = (3*f1 - 2*f2 - f3);
    t22 = (4*f1 - 3*f2 - f3); t23 = (3*f1 - f2 - 2*f3); t24 = (2*f1 - f2 - f3); t25 = (2*f1 + f2 - 3*f3);t26 = (4*f1 - f2 - 3*f3);
    t27 = (4*f2 - 3*f3 - f4); t28 = (3*f2 - 2*f3 - f4); t29 = (3*f1 - 4*f3 + f4); t30 = (f1 - 2*f3 + f4);t31 = (3*f1 - 4*f2 + f4);
    t32 = (f1 - 4*f2 + 3*f4); t33 = (f1 + f2 - 2*f3); t34 = (f1 + 2*f2 - 3*f3); t35 = (f2 + 2*f3 - 3*f4); t36 = (f1 - 4*f3 + 3*f4);
    t37 = (f1 - 3*f3 + 2*f4); t38 = (f1 + 3*f2 - 4*f3); t39 = (3*f1 + f2 - 4*f3); t40 = (f1 + 3*f2 - 4*f4); t41 = (f2 + 3*f3 - 4*f4);
    t42 = (3*f1 + f2 - 4*f4);
    
    p1 = (3*(12*t1*g1 + t3*g2 + t1*g3 + t4*g4)*h1 + (3*t3*g1 + 6*t5*g2 + t6*g3 + 3*t1*g4)*h2 + (3*t1*g1 + t6*g2 + 2*t1*g3 + t7*g4)*h3 + (3*t4*g1 + 3*t1*g2 + t7*g3 + 6*t8*g4)*h4);
    p2 = ((9*t9*g1 + 6*t10*g2 + t11*g3 + 3*t12*g4)*h1 + (6*t10*g1 + 3*t13*g2 + 3*t2*g3 + t11*g4)*h2 + (t11*g1 + 3*t2*g2 + t14*g3 + 2*t2*g4)*h3 + (3*t12*g1 + t11*g2 + 2*t2*g3 + 3*t10*g4)*h4);
    p3 = ((9*t1*g1 - 3*t15*g2 + 2*t1*g3 + 3*t16*g4)*h1 - (3*t15*g1 + 3*t17*g2 + t18*g3 - 2*t1*g4)*h2 + (2*t1*g1 - t18*g2 + 3*t1*g3 + t19*g4)*h3 + (3*t16*g1 + 2*t1*g2 + t19*g3 + 3*t20*g4)*h4);
    p4 = ((9*t22*g1 + 3*t21*g2 + t23*g3 + 6*t24*g4)*h1 + (3*t21*g1 + 3*t24*g2 + 2*t2*g3 + t23*g4)*h2 + (t23*g1 + 2*t2*g2 + t25*g3 + 3*t2*g4)*h3 + (6*t24*g1 + t23*g2 + 3*t2*g3 + 3*t26*g4)*h4);
    
    p5 = ((3*t3*g1 + 6*t5*g2 + t6*g3 + 3*t1*g4)*h1 + (6*t5*g1 + 9*t27*g2 + 3*t28*g3 + t6*g4)*h2 + (t6*g1 + 3*t28*g2 + 3*t5*g3 + 2*t1*g4)*h3 + (3*t1*g1 + t6*g2 + 2*t1*g3 + t7*g4)*h4);
    p6 = ((6*t10*g1 + 3*t13*g2 + 3*t2*g3 + t11*g4)*h1 + 3*(t13*g1 + 12*t2*g2 + t29*g3 + t2*g4)*h2 + (3*t2*g1 + 3*t29*g2 + 6*t30*g3 + t14*g4)*h3 + (t11*g1 + 3*t2*g2 + t14*g3 + 2*t2*g4)*h4);
    p7 = ((3*t15*g1 + 3*t17*g2 + t18*g3 - 2*t1*g4)*h1 + (3*t17*g1 + 9*t31*g2 + 6*t15*g3 + t18*g4)*h2 + (t18*g1 + 6*t15*g2 + 3*t32*g3 - 3*t1*g4)*h3 - (2*t1*g1 - t18*g2 + 3*t1*g3 + t19*g4)*h4);
    p8 = ((3*t21*g1 + 3*t24*g2 + 2*t2*g3 + t23*g4)*h1 + (3*t24*g1 + 9*t2*g2 + 3*t33*g3 + 2*t2*g4)*h2 + (2*t2*g1 + 3*t33*g2 + 3*t34*g3 + t25*g4)*h3 + (t23*g1 + 2*t2*g2 + t25*g3 + 3*t2*g4)*h4);
    
    d1 = ((3*t1*g1 + t6*g2 + 2*t1*g3 + t7*g4)*h1 + (t6*g1 + 3*t28*g2 + 3*t5*g3 + 2*t1*g4)*h2 + (2*t1*g1 + 3*t5*g2 + 9*t1*g3 + 3*t8*g4)*h3 + (t7*g1 + 2*t1*g2 + 3*t8*g3 + 3*t35*g4)*h4);
    d2 = ((t11*g1 + 3*t2*g2 + t14*g3 + 2*t2*g4)*h1 + (3*t2*g1 + 3*t29*g2 + 6*t30*g3 + t14*g4)*h2 + (t14*g1 + 6*t30*g2 + 9*t36*g3 + 3*t37*g4)*h3 + (2*t2*g1 + t14*g2 + 3*t37*g3 + 3*t30*g4)*h4);
    d3 = ((2*t1*g1 - t18*g2 + 3*t1*g3 + t19*g4)*h1 - (t18*g1 + 6*t15*g2 + 3*t32*g3 - 3*t1*g4)*h2 + 3*(t1*g1 - t32*g2 + 12*t1*g3 + t40*g4)*h3 + (t19*g1 + 3*t1*g2 + 3*t40*g3 + 6*t16*g4)*h4);
    d4 = ((t23*g1 + 2*t2*g2 + t25*g3 + 3*t2*g4)*h1 + (2*t2*g1 + 3*t33*g2 + 3*t34*g3 + t25*g4)*h2 + (t25*g1 + 3*t34*g2 + 9*t38*g3 + 6*t33*g4)*h3 + (3*t2*g1 + t25*g2 + 6*t33*g3 + 3*t39*g4)*h4);
    
    d5 = ((3*t4*g1 + 3*t1*g2 + t7*g3 + 6*t8*g4)*h1 + (3*t1*g1 + t6*g2 + 2*t1*g3 + t7*g4)*h2 + (t7*g1 + 2*t1*g2 + 3*t8*g3 + 3*t35*g4)*h3 + (6*t8*g1 + t7*g2 + 3*t35*g3 + 9*t41*g4)*h4);
    d6 = ((3*t12*g1 + t11*g2 + 2*t2*g3 + 3*t10*g4)*h1 + (t11*g1 + 3*t2*g2 + t14*g3 + 2*t2*g4)*h2 + (2*t2*g1 + t14*g2 + 3*t37*g3 + 3*t30*g4)*h3 + (3*t10*g1 + 2*t2*g2 + 3*t30*g3 + 9*t2*g4)*h4);
    d7 = ((3*t16*g1 + 2*t1*g2 + t19*g3 + 3*t20*g4)*h1 + (2*t1*g1 - t18*g2 + 3*t1*g3 + t19*g4)*h2 + (t19*g1 + 3*t1*g2 + 3*t40*g3 + 6*t16*g4)*h3 + (3*t20*g1 + t19*g2 + 6*t16*g3 + 9*t42*g4)*h4);
    d8 = ((6*t24*g1 + t23*g2 + 3*t2*g3 + 3*t26*g4)*h1 + (t23*g1 + 2*t2*g2 + t25*g3 + 3*t2*g4)*h2 + (3*t2*g1 + t25*g2 + 6*t33*g3 + 3*t39*g4)*h3 + 3*(t26*g1 + t2*g2 + t39*g3 + 12*t2*g4)*h4);
    
    con = one_720 * c;
    integral_x[0] += con * (-p1*y1 + p2*y2 + p3*y3 - p4*y4);
    integral_x[1] += con * (-p5*y1 + p6*y2 - p7*y3 - p8*y4);
    integral_x[2] += con * (-d1*y1 + d2*y2 + d3*y3 - d4*y4);
    integral_x[3] += con * (-d5*y1 + d6*y2 + d7*y3 - d8*y4);
    
    integral_y[0] += con * ( p1*x1 - p2*x2 - p3*x3 + p4*x4);
    integral_y[1] += con * ( p5*x1 - p6*x2 + p7*x3 + p8*x4);
    integral_y[2] += con * ( d1*x1 - d2*x2 - d3*x3 + d4*x4);
    integral_y[3] += con * ( d5*x1 - d6*x2 - d7*x3 + d8*x4);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following quadrilateral integration:
 *  \f$ \int_{-1}^{1} \int_{-1}^{1} \big( c \, \mathbf{\nabla} \phi_i(\widehat{x},\widehat{y}) \cdot \mathbf{\overline{v}}  \big) d\widehat{y}\,d\widehat{x}   \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_quadrilateral_gradPhi_dot_vbar(SVECT *nd, double c, SVECT2D vbar, double *integral) {
    
    double con = 6. * c;
    integral[0] += con * (-vbar.y * nd[1].x + vbar.y * nd[3].x + vbar.x * nd[1].y - vbar.x * nd[3].y);
    integral[1] += con * ( vbar.y * nd[0].x - vbar.y * nd[2].x - vbar.x * nd[0].y + vbar.x * nd[2].y);
    integral[2] += con * ( vbar.y * nd[1].x - vbar.y * nd[3].x - vbar.x * nd[1].y + vbar.x * nd[3].y);
    integral[3] += con * (-vbar.y * nd[0].x + vbar.y * nd[2].x + vbar.x * nd[0].y - vbar.x * nd[2].y);
    
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following quadrilateral integration:
 *  \f{eqnarray*}{
 *   integral_x &=& \int_{-1}^{1} \int_{-1}^{1} \bigg(c \, \phi_i(\widehat{x},\widehat{y}) \, \frac{\partial \f(\widehat{x},\widehat{y})}{x} \, h(\widehat{x},\widehat{y}) \bigg) d\widehat{y}\,d\widehat{x} \\
 *   integral_y &=& \int_{-1}^{1} \int_{-1}^{1} \bigg(c \, \phi_i(\widehat{x},\widehat{y}) \, \frac{\partial \f(\widehat{x},\widehat{y})}{y} \, h(\widehat{x},\widehat{y}) \bigg) d\widehat{y}\,d\widehat{x}
 *  \f}
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: SW 2D pressure body integral
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_quadrilateral_phi_h_df(SVECT *nd, double c, double *df, double *h, double *integral_x, double *integral_y) {
    
    double x1, x2, x3, x4, y1, y2, y3, y4, f1, f2, f3, f4, h1, h2, h3, h4;
    double g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14, g15, g16;
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26;
    double con = 0.;
    
    extractNodesQuad(nd,&x1,&x2,&x3,&x4,&y1,&y2,&y3,&y4);
    extractFunctionQuad(df,&f1,&f2,&f3,&f4);
    extractFunctionQuad(h,&h1,&h2,&h3,&h4);
    
    t1 = f2 - f4; t2 = 3*f2 - f3 - 2*f4; t3 = 2*f2 + f3 - 3*f4; t4 = 3*f1 - f3 - 2*f4; t5 = 3*f1 - 2*f3 - f4;
    t6 = f1 - f3; t7 = 2*f1 - f3 - f4; t8 = f1 - 2*f2 + f4; t9 = f1 + f2 - 2*f4; t10 = 3*f1 - 2*f2 - f3;
    t11 = 2*f1 - f2 - f3; t12 = 3*f2 - 2*f3 - f4; t13 = 2*f2 - f3 - f4; t14 = 2*f1 - 3*f3 + f4; t15 = 2*f1 - 3*f2 + f4;
    t16 = f1 - 3*f2 + 2*f4; t17 = f1 + f2 - 2*f3; t18 = f2 + f3 - 2*f4; t19 = f1 - 3*f3 + 2*f4; t20 = f1 - 2*f3 + f4;
    t21 = f1 + 2*f2 - 3*f4; t22 = f1 + 2*f2 - 3*f3; t23 = 2*f1 + f2 - 3*f3; t24 = f2 + 2*f3 - 3*f4; t25 = 2*f1 + f2 - 3*f4;
    t26 = 3*f1 - f2 - 2*f3;
    
    g1 = (6*t1*h1 + t2*h2 + t1*h3 + t3*h4);    g5 = (2*t4*h1 + t5*h2 + t6*h3 + t7*h4);
    g2 = (t2*h1 + 2*t12*h2 + t13*h3 + t1*h4);  g6 = (t5*h1 + 6*t6*h2 + t14*h3 + t6*h4);
    g3 = (t1*h1 + t13*h2 + 2*t1*h3 + t18*h4);  g7 = (t6*h1 + t14*h2 + 2*t19*h3 + t20*h4);
    g4 = (t3*h1 + t1*h2 + t18*h3 + 2*t24*h4);  g8 = (t7*h1 + t6*h2 + t20*h3 + 2*t6*h4);
    
    g9 = (2*t1*h1 - t8*h2 + t1*h3 + t9*h4);    g13 = (2*t10*h1 + t11*h2 + t6*h3 + t26*h4);
    g10 = (t8*h1 + 2*t15*h2 + t16*h3 - t1*h4); g14 = (t11*h1 + 2*t6*h2 + t17*h3 + t6*h4);
    g11 = (t1*h1 - t16*h2 + 6*t1*h3 + t21*h4); g15 = (t6*h1 + t17*h2 + 2*t22*h3 + t23*h4);
    g12 = (t9*h1 + t1*h2 + t21*h3 + 2*t25*h4); g16 = (t26*h1 + t6*h2 + t23*h3 + 6*t6*h4);
    
    con = one_72 * c;
    if (integral_x != NULL) {
        integral_x[0] += con * ( -g1*y1 + g5*y2 + g9*y3  - g13*y4);
        integral_x[1] += con * ( -g2*y1 + g6*y2 - g10*y3 - g14*y4);
        integral_x[2] += con * ( -g3*y1 + g7*y2 + g11*y3 - g15*y4);
        integral_x[3] += con * ( -g4*y1 + g8*y2 + g12*y3 - g16*y4);
    }
    
    if (integral_y != NULL) {
        integral_y[0] += con * ( g1*x1 - g5*x2 - g9*x3  + g13*x4);
        integral_y[1] += con * ( g2*x1 - g6*x2 + g10*x3 + g14*x4);
        integral_y[2] += con * ( g3*x1 - g7*x2 - g11*x3 + g15*x4);
        integral_y[3] += con * ( g4*x1 - g8*x2 - g12*x3 + g16*x4);
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                             VERTICAL QUADRILATERALS IN 3D SPACE                          */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following *vertical quadrilateral integration:
 *  \f$ \int_{-1}^{1} \int_{-1}^{1} \big( c \, f(\widehat{s},\widehat{t}) \big) d\widehat{s}\,d\widehat{t}   \f$
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: for this integration to be correct, grid must be in column format and numbered correctly.
 *                (i.e. x3 = x2, y3 = y2, x4 = x1, y4 = y1)
 *  \note CJT \:: negative on t1 here because of ordering of nodes
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline double integrate_quadZ_f(SVECT *nd, double c, double *f) {
#if _DEBUG
    assert(nd[0].x == nd[3].x);
    assert(nd[1].x == nd[2].x);
    assert(nd[0].y == nd[3].y);
    assert(nd[1].y == nd[2].y);
#endif
    double t1 = one_2 * sqrt(pow(nd[0].x - nd[1].x,2) + pow(nd[0].y - nd[1].y,2));
    double t2 = -t1 * one_6 * c; // negative was not here in old code ...
    return (t2 * ((2*f[0] +   f[1] +   f[2] + 2*f[3])*nd[0].z +
                  (  f[0] + 2*f[1] + 2*f[2] +   f[3])*nd[1].z -
                  (  f[0] + 2*f[1] + 2*f[2] +   f[3])*nd[2].z -
                  (2*f[0] +   f[1] +   f[2] + 2*f[3])*nd[3].z));
}

