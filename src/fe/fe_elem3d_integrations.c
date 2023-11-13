/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  \brief fe_elem3d_integrations.c This file collects functions that integrate over 3d AdH elements. */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                       TETRAHEDRONS                                       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following tetrahedral integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{0}^{1-\widehat{x}-\widehat{y}}
 c \,
 *       \, d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline double integrate_tetrahedron_c(double djac, double c) {
    return (c * djac);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following tetrahedral integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{0}^{1-\widehat{x}-\widehat{y}}
 *          c\, f_i \,
 *       \, d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_tetrahedron_fi(double djac, double c, double *f, double *integral) {
    int i;
    for (i=0; i<NDONTET; i++) {
        integral[i] = djac * c * f[i];
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following tetrahedral integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{0}^{1-\widehat{x}-\widehat{y}}
 *          c \, f(\widehat{x},\widehat{y},\widehat{z}) \,
 *      \, d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline double integrate_tetrahedron_f(double djac, double c, double *f) {
    return (c * djac * one_4 * sarray_sum_dbl(f,NDONTET));
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following tetrahedral integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{0}^{1-\widehat{x}-\widehat{y}}
 *          c \, f(\widehat{x},\widehat{y},\widehat{z}) \, g(\widehat{x},\widehat{y},\widehat{z}) \,
 *      \, d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline double integrate_tetrahedron_f_g(double djac, double c, double *f, double *g) {
    double t1 = djac * one_20 * c;
    
    //double sum_f = f[0] + f[1] + f[2] + f[3];
    //double sum_g = g[0] + g[1] + g[2] + g[3];
    //double t2 = sum_f * sum_g + f[0] * g[0] + f[1] * g[1] + f[2] * g[2] + f[3] * g[3];
    
    // CJT :: there is nothing wrong with this I can find, just significantly different than trunk results
    return (t1 * ( (2*f[0] + f[1] +   f[2] + f[3])*g[0] + (f[0] + 2*f[1] + f[2] +   f[3])*g[1] +
                     (f[0] + f[1] + 2*f[2] + f[3])*g[2] + (f[0] +   f[1] + f[2] + 2*f[3])*g[3] ) );
    
    //return (t1 * t2);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following tetrahedral integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{0}^{1-\widehat{x}-\widehat{y}}
 *          c \, \omega_i(\widehat{x},\widehat{y},\widehat{z}) \, f(\widehat{x},\widehat{y},\widehat{z}) \,
 *      \, d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_tetrahedron_phi_f(double djac, double c, double *f, double *integral) {
    int i;
    double t1 = djac * one_20 * c;
    double sum_f = sarray_sum_dbl(f,NDONTET);
    for (i=0; i<NDONTET; i++) {
        integral[i] = t1 * (sum_f + f[i]);
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following tetrahedral integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{0}^{1-\widehat{x}-\widehat{y}}
 *          c \, f(\widehat{x},\widehat{y},\widehat{z}) \, g(\widehat{x},\widehat{y},\widehat{z}) \,
 *      \, d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_tetrahedron_phi_f_g(double djac, double c, double *f, double *g, double *integral) {
    int i;
    double t1 = djac * one_120 * c;
    double a1 = 2*f[0] +   f[1] +   f[2] + 2*f[3], a2 = 2*f[0] + 2*f[1] +   f[2] +   f[3];
    double a3 = 2*f[0] +   f[1] + 2*f[2] +   f[3], a4 =   f[0] + 2*f[1] +   f[2] + 2*f[3];
    double a5 =   f[0] + 2*f[1] + 2*f[2] +   f[3], a6 =   f[0] +   f[1] + 2*f[2] + 2*f[3];
    
    integral[0] = t1 * ( 2*(3*f[0] + f[1] + f[2] + f[3])*g[0] + a2*g[1] + a3*g[2] + a1*g[3] );
    integral[1] = t1 * ( a2*g[0] + 2*(f[0] + 3*f[1] + f[2] + f[3])*g[1] + a5*g[2] + a4*g[3] );
    integral[2] = t1 * ( a3*g[0] + a5*g[1] + 2*(f[0] + f[1] + 3*f[2] + f[3])*g[2] + a6*g[3] );
    integral[3] = t1 * ( a1*g[0] + a4*g[1] + a6*g[2] + 2*(f[0] + f[1] + f[2] + 3*f[3])*g[3] );
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following tetrahedral integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{0}^{1-\widehat{x}-\widehat{y}}
 *          c_i \, f(\widehat{x},\widehat{y},\widehat{z}) \, g(\widehat{x},\widehat{y},\widehat{z}) \,
 *      \, d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_tetrahedron_ci_f_g(double djac, double *c, double *f, double *g, double *integral) {
    int i;
    double sum_f = sarray_sum_dbl(f,NDONTET);
    double t1 = djac * one_20 * ((sum_f + f[0])*g[0]+ (sum_f + f[1])*g[1]+ (sum_f + f[2])*g[2]+ (sum_f + f[3])*g[3]);
    for (i=0; i<NDONTET; i++) {
        integral[i] = t1 * c[i];
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following tetrahedral integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{0}^{1-\widehat{x}-\widehat{y}}
 *          c \, \mathbf{\nabla} \omega_i(\widehat{x},\widehat{y},\widehat{z}) \cdot \mathbf{v}(\widehat{x},\widehat{y},\widehat{z}) \,
 *      \, d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_tetrahedron_gradPhi_dot_v(SVECT *grad_phi, double djac, double c, SVECT *vector, double *integral) {
    
    int i;
    double u[NDONTET], v[NDONTET], w[NDONTET];
    for (i=0; i<NDONTET; i++) {
        u[i] = vector[i].x;
        v[i] = vector[i].y;
        w[i] = vector[i].z;
    }
    for (i=0; i<NDONTET; i++) {
        integral[i] = integrate_tetrahedron_f(djac, c * grad_phi[i].x, u) +
                      integrate_tetrahedron_f(djac, c * grad_phi[i].y, v) +
                      integrate_tetrahedron_f(djac, c * grad_phi[i].z, w);
    }
}

inline void integrate_tetrahedron_gradPhi_dot_vcon(SVECT *grad_phi, double djac, double c, SVECT vector, double *integral) {
    int i;
    double t1 = djac * c;
    for (i=0; i<NDONTET; i++) {
        integral[i] = t1 * (grad_phi[i].x * vector.x + grad_phi[i].y * vector.y + grad_phi[i].z * vector.z);
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following tetrahedral integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{0}^{1-\widehat{x}-\widehat{y}}
 *          c \, \mathbf{\nabla} \omega_i(\widehat{r}) \cdot \bigg(f(\widehat{r}) \mathbf{v}(\widehat{r})\bigg) \,
 *       d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_tetrahedron_gradPhi_dot_f_v(SVECT *grad_phi, double djac, double c, double *f, SVECT *vector, double *integral) {
    
    int i;
    double u[NDONTET], v[NDONTET], w[NDONTET];
    for (i=0; i<NDONTET; i++) {
        u[i] = vector[i].x;
        v[i] = vector[i].y;
        w[i] = vector[i].z;
    }
    double t1 = djac * c;
    for (i=0; i<NDONTET; i++) {
        integral[i] = t1 * (
        grad_phi[i].x * integrate_tetrahedron_f_g(1., 1., f, u) +
        grad_phi[i].y * integrate_tetrahedron_f_g(1., 1., f, v) +
        grad_phi[i].z * integrate_tetrahedron_f_g(1., 1., f, w));
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                     TRIANGULAR PRISMS                                    */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Checks the location so prism nodes for compatibility with integration routines
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: nodes must be in columns for these routines
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void check_prism_nodes(SVECT *nd, int linenumber, char *filename) {
    int flag = 0;
    if (nd[0].x != nd[3].x || nd[0].y != nd[3].y) {
        printf("\nERROR IN PRISM NODE ORDER :: FILE: %s LINE: %d\n",filename,linenumber);
        printf("-- Node 1 and node 4 are not in a column!\n");
        flag = 1;
    }
    if (nd[1].x != nd[4].x || nd[1].y != nd[4].y) {
        printf("\nERROR IN PRISM NODE ORDER :: FILE: %s LINE: %d\n",filename,linenumber);
        printf("-- Node 2 and node 5 are not in a column!\n");
        flag = 1;
    }
    if (nd[2].x != nd[5].x || nd[2].y != nd[5].y) {
        printf("\nERROR IN PRISM NODE ORDER :: FILE: %s LINE: %d\n",filename,linenumber);
        printf("-- Node 3 and node 6 are not in a column!\n");
        flag = 1;
    }
    if (flag == 1) {
        svect_printScreen(nd[0],"node 1: ");
        svect_printScreen(nd[1],"node 2: ");
        svect_printScreen(nd[2],"node 3: ");
        svect_printScreen(nd[3],"node 4: ");
        svect_printScreen(nd[4],"node 5: ");
        svect_printScreen(nd[5],"node 6: ");
        exit(-1);
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular prism integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1}
 *          c \,
 *       \, d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: assumed top and bottom triangles have same x,y coordinates
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline double integrate_triPrism_area(SVECT *nd, double c) {
#if _DEBUG
    check_prism_nodes(nd,__LINE__,__FILE__);
#endif
    double t1 = one_6 * (nd[1].x*nd[0].y - nd[2].x*nd[0].y - (nd[0].x - nd[2].x)*nd[1].y + (nd[0].x - nd[1].x)*nd[2].y);
    double t2 = t1 * c;
    return ( t2 * (nd[0].z + nd[1].z + nd[2].z - nd[3].z - nd[4].z - nd[5].z) );
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular prism integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1}
 *          c \, f(\widehat{x},\widehat{y},\widehat{z})
 *       \, d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: assumed top and bottom triangles have same x,y coordinates
 *  \note CJT \:: the 3D coriolis test case is agitatingly sensitive to the form of this calculation.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline double integrate_triPrism_f(SVECT *nd, double c, double *f) {
#if _DEBUG
    check_prism_nodes(nd,__LINE__,__FILE__);
#endif
    double x1 = nd[0].x, x2 = nd[1].x, x3 = nd[2].x;
    double y1 = nd[0].y, y2 = nd[1].y, y3 = nd[2].y;
    double z1 = nd[0].z, z2 = nd[1].z, z3 = nd[2].z, z4 = nd[3].z, z5 = nd[4].z, z6 = nd[5].z;
    
    double f211, f_211, f121, f_121, f112, f_112, a1, a2, a3;
    
    f211 = 2*f[0]+   f[1]+   f[2]; f_211 = 2*f[3]+   f[4]+   f[5];
    f121 =   f[0]+ 2*f[1]+   f[2]; f_121 =   f[3]+ 2*f[4]+   f[5];
    f112 =   f[0]+   f[1]+ 2*f[2]; f_112 =   f[3]+   f[4]+ 2*f[5];
    a1 = f211 + f_211;
    a2 = f121 + f_121;
    a3 = f112 + f_112;
    
    double t1 = c * one_48 * ((x2 - x3)*y1 - (x1 - x3)*y2 + (x1 - x2)*y3);
    return (t1 * (a1*z1 + a2*z2 + a3*z3 - a1*z4 - a2*z5 - a3*z6));
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular prism integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1}
 *          c \, \omega_i(\widehat{x},\widehat{y},\widehat{z})
 *       \, d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: assumed top and bottom triangles have same x,y coordinates
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triPrism_phi(SVECT *nd, double c, double *integral) {
#if _DEBUG
    check_prism_nodes(nd,__LINE__,__FILE__);
#endif
    double t1 = one_6 * (nd[1].x*nd[0].y - nd[2].x*nd[0].y - (nd[0].x - nd[2].x)*nd[1].y + (nd[0].x - nd[1].x)*nd[2].y);
    double t2 = t1 * c * one_8;
    
    double p1, p2, p3;
    p1 = 2*nd[0].z +   nd[1].z +   nd[2].z - 2*nd[3].z -   nd[4].z -   nd[5].z;
    p2 =   nd[0].z + 2*nd[1].z +   nd[2].z -   nd[3].z - 2*nd[4].z -   nd[5].z;
    p3 =   nd[0].z +   nd[1].z + 2*nd[2].z -   nd[3].z -   nd[4].z - 2*nd[5].z;
    
    integral[0] = t2 * p1;
    integral[1] = t2 * p2;
    integral[2] = t2 * p3;
    integral[3] = integral[0];
    integral[4] = integral[1];
    integral[5] = integral[2];
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular prism integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1}
 *          c \, \omega_i(\widehat{x},\widehat{y},\widehat{z}) \, f(\widehat{x},\widehat{y},\widehat{z})
 *       \, d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: assumed top and bottom triangles have same x,y coordinates
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triPrism_phi_f(SVECT *nd, double c, double *f, double *integral) {
#if _DEBUG
    check_prism_nodes(nd,__LINE__,__FILE__);
#endif
    double t1 = one_6 * (nd[1].x*nd[0].y - nd[2].x*nd[0].y - (nd[0].x - nd[2].x)*nd[1].y + (nd[0].x - nd[1].x)*nd[2].y);
    double t2 = t1 * c * one_120;
    
    double p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12;
    double f311,f_311,f131,f_131,f113,f_113,f221,f_221,f212,f_212,f122,f_122;
    
    f311 = 3*f[0]+   f[1]+   f[2]; f_311 = 3*f[3]+   f[4]+   f[5];
    f131 =   f[0]+ 3*f[1]+   f[2]; f_131 =   f[3]+ 3*f[4]+   f[5];
    f113 =   f[0]+   f[1]+ 3*f[2]; f_113 =   f[3]+   f[4]+ 3*f[5];
    f221 = 2*f[0]+ 2*f[1]+   f[2]; f_221 = 2*f[3]+ 2*f[4]+   f[5];
    f212 = 2*f[0]+   f[1]+ 2*f[2]; f_212 = 2*f[3]+   f[4]+ 2*f[5];
    f122 =   f[0]+ 2*f[1]+ 2*f[2]; f_122 =   f[3]+ 2*f[4]+ 2*f[5];
    
    p1  = 2*f311 + f_311; p7  = f311 + 2*f_311;
    p4  = 2*f131 + f_131; p10 = f131 + 2*f_131;
    p6  = 2*f113 + f_113; p12 = f113 + 2*f_113;
    p2  = 2*f221 + f_221; p8  = f221 + 2*f_221;
    p3  = 2*f212 + f_212; p9  = f212 + 2*f_212;
    p5  = 2*f122 + f_122; p11 = f122 + 2*f_122;
    
    integral[0] = t2 * (2*p1*nd[0].z +    p2*nd[1].z +    p3*nd[2].z - 2*p1*nd[3].z -    p2*nd[4].z -    p3*nd[5].z);
    integral[1] = t2 * (  p2*nd[0].z +  2*p4*nd[1].z +    p5*nd[2].z -   p2*nd[3].z -  2*p4*nd[4].z -    p5*nd[5].z);
    integral[2] = t2 * (  p3*nd[0].z +    p5*nd[1].z +  2*p6*nd[2].z -   p3*nd[3].z -    p5*nd[4].z -  2*p6*nd[5].z);
    integral[3] = t2 * (2*p7*nd[0].z +    p8*nd[1].z +    p9*nd[2].z - 2*p7*nd[3].z -    p8*nd[4].z -    p9*nd[5].z);
    integral[4] = t2 * (  p8*nd[0].z + 2*p10*nd[1].z +   p11*nd[2].z -   p8*nd[3].z - 2*p10*nd[4].z -   p11*nd[5].z);
    integral[5] = t2 * (  p9*nd[0].z +   p11*nd[1].z + 2*p12*nd[2].z -   p9*nd[3].z -   p11*nd[4].z - 2*p12*nd[5].z);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular prism integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1}
 *          c \, \frac{\partial f(\widehat{x},\widehat{y},\widehat{z})}{\partial x,y,z}
 *       \, d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: only one directional derivative is returned here
 *  \note CJT \:: assumed top and bottom triangles have same x,y coordinates
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline double integrate_triPrism_df(SVECT *nd, double c, double *f, int direction) {
#if _DEBUG
    check_prism_nodes(nd,__LINE__,__FILE__);
#endif
    double t1, t2, df = 0.;
    
    if (direction == 3) { // df/dz
        t1 = one_6 * (nd[1].x*nd[0].y - nd[2].x*nd[0].y - (nd[0].x - nd[2].x)*nd[1].y + (nd[0].x - nd[1].x)*nd[2].y);
        df = t1 * c * (f[0] + f[1] + f[2] - f[3] - f[4] - f[5]);
    } else if (direction == 1 || direction == 2) { // df/dx,y
        double a1, a2, a3, a4, a5, a6, a7, a8, a9;
        a1 =   f[1] -   f[2] +   f[4] -   f[5];
        a2 =   f[1] + 2*f[2] - 2*f[3] -   f[4];
        a3 = 2*f[1] +   f[2] - 2*f[3] -   f[5];
        a4 =   f[0] + 2*f[2] -   f[3] - 2*f[4];
        a5 =   f[0] -   f[2] +   f[3] -   f[5];
        a6 = 2*f[0] +   f[2] - 2*f[4] -   f[5];
        a7 =   f[0] + 2*f[1] -   f[3] - 2*f[5];
        a8 = 2*f[0] +   f[1] -   f[4] - 2*f[5];
        a9 =   f[0] -   f[1] +   f[3] -   f[4];
        
        double c1, c2, c3, c4, c5, c6;
        if (direction == 2) { // y-direction
            t1 = -one_12 * c;
            c1 = nd[0].x; c2 = nd[1].x; c3 = nd[2].x;
        } else {              // x-direction
            t1 =  one_12 * c;
            c1 = nd[0].y; c2 = nd[1].y; c3 = nd[2].y;
        }
        df = t1 * ((a1*c1 + a2*c2 - a3*c3)*nd[0].z - (a4*c1 + a5*c2 - a6*c3)*nd[1].z + (a7*c1 - a8*c2 + a9*c3)*nd[2].z -
        (a1*c1 - a8*c2 + a6*c3)*nd[3].z - (a7*c1 - a5*c2 - a3*c3)*nd[4].z + (a4*c1 - a2*c2 - a9*c3)*nd[5].z);
    } else {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> direction must be 1 = x, 2 = y, or 3 = z\n");
    }
    
    return df;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular prism integration:
 *  \f{eqnarray*}{
 * integral_x &=& \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1} c \, \deriv{f(\widehat{x},\widehat{y},\widehat{z})}{x} \, d\Omega \\
 * integral_y &=& \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1} c \, \deriv{f(\widehat{x},\widehat{y},\widehat{z})}{y} \, d\Omega \\
 * integral_z &=& \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1} c \, \deriv{f(\widehat{x},\widehat{y},\widehat{z})}{z} \, d\Omega \\
 *  \f}
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline SVECT integrate_triPrism_df_full(SVECT *nd, double c, double *f) {
#if _DEBUG
    check_prism_nodes(nd,__LINE__,__FILE__);
#endif
    
    double x1 = nd[0].x, x2 = nd[1].x, x3 = nd[2].x;
    double y1 = nd[0].y, y2 = nd[1].y, y3 = nd[2].y;
    double z1 = nd[0].z, z2 = nd[1].z, z3 = nd[2].z, z4 = nd[3].z, z5 = nd[4].z, z6 = nd[5].z;
    
    double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, t1;
    a1 = f[1]- f[2]+ f[4]- f[5];
    a2 = f[1]+ 2*f[2]- 2*f[3]- f[4];
    a3 = 2*f[1]+ f[2]- 2*f[3]- f[5];
    a4 = f[0]+ 2*f[2]- f[3]- 2*f[4];
    a5 = f[0]- f[2]+ f[3]- f[5];
    a6 = 2*f[0]+ f[2]- 2*f[4]- f[5];
    a7 = f[0]+ 2*f[1]- f[3]- 2*f[5];
    a8 = 2*f[0]+ f[1]- f[4]- 2*f[5];
    a9 = f[0]- f[1]+ f[3]- f[4];
    a10 = f[0]+ f[1]+ f[2]- f[3]- f[4]- f[5];
    
    SVECT integral;
    
    t1 = one_12 * c;
    integral.x =  t1*( (a1*y1 + a2*y2 - a3*y3)*z1 - (a4*y1 + a5*y2 - a6*y3)*z2 + (a7*y1 - a8*y2 + a9*y3)*z3 - (a1*y1 - a8*y2 + a6*y3)*z4 - (a7*y1 - a5*y2 - a3*y3)*z5 + (a4*y1 - a2*y2 - a9*y3)*z6);
    integral.y =  t1*(-(a1*x1 + a2*x2 - a3*x3)*z1 + (a4*x1 + a5*x2 - a6*x3)*z2 - (a7*x1 - a8*x2 + a9*x3)*z3 + (a1*x1 - a8*x2 + a6*x3)*z4 + (a7*x1 - a5*x2 - a3*x3)*z5 - (a4*x1 - a2*x2 - a9*x3)*z6);
    
    t1 = one_6 * ((nd[1].x - nd[2].x)*nd[0].y - (nd[0].x - nd[2].x)*nd[1].y + (nd[0].x - nd[1].x)*nd[2].y);
    //integral.z =  t1* 2 * a10 * ((x2 - x3)*y1 - 2*(x1 - x3)*y2 + 2*(x1 - x2)*y3);
    integral.z = t1 * c * a10;
    
    return integral;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular prism integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1}
 *          c \, \frac{\partial \omega_i(\widehat{x},\widehat{y},\widehat{z})}{\partial x,y,z}
 *       \, d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: derivative in all directions are returned here
 *  \note CJT \:: assumed top and bottom triangles have same x,y coordinates
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triPrism_dphi(SVECT *nd, double c, double *integral, int direction) {
#if _DEBUG
    check_prism_nodes(nd,__LINE__,__FILE__);
#endif
    
    double t1, t2;
    
    if (direction == 3) { // dphi_i/dz
        t1 = one_6 * (nd[1].x*nd[0].y - nd[2].x*nd[0].y - (nd[0].x - nd[2].x)*nd[1].y + (nd[0].x - nd[1].x)*nd[2].y);
        t2 = c * t1;
        integral[0] =  t2;
        integral[1] =  t2;
        integral[2] =  t2;
        integral[3] = -t2;
        integral[4] = -t2;
        integral[5] = -t2;
    } else if (direction == 1 || direction == 2) { // dphi_i/dx,y
        double c1, c2, c3, c4, c5, c6;
        if (direction == 2) {  // y-direction
            t2 = -one_12 * c;
            c1 = nd[0].x; c2 = nd[1].x; c3 = nd[2].x;
        } else {
            t2 =  one_12 * c;  // x-direction
            c1 = nd[0].y; c2 = nd[1].y; c3 = nd[2].y;
        }
        
        integral[0] = t2 * ( -(c1 + c2 - 2*c3)*nd[1].z + (c1 - 2*c2 + c3)*nd[2].z + 2*(c2 - c3)*nd[3].z -        (c1 - c2)*nd[4].z +        (c1 - c3)*nd[5].z );
        integral[1] = t2 * (  (c1 + c2 - 2*c3)*nd[0].z + (2*c1 - c2 - c3)*nd[2].z -   (c1 - c2)*nd[3].z -      2*(c1 - c3)*nd[4].z -        (c2 - c3)*nd[5].z );
        integral[2] = t2 * ( -(c1 - 2*c2 + c3)*nd[0].z - (2*c1 - c2 - c3)*nd[1].z +   (c1 - c3)*nd[3].z -        (c2 - c3)*nd[4].z +      2*(c1 - c2)*nd[5].z );
        integral[3] = t2 * (      -2*(c2 - c3)*nd[0].z +        (c1 - c2)*nd[1].z -   (c1 - c3)*nd[2].z + (c1 + c2 - 2*c3)*nd[4].z - (c1 - 2*c2 + c3)*nd[5].z );
        integral[4] = t2 * (         (c1 - c2)*nd[0].z +      2*(c1 - c3)*nd[1].z +   (c2 - c3)*nd[2].z - (c1 + c2 - 2*c3)*nd[3].z - (2*c1 - c2 - c3)*nd[5].z );
        integral[5] = t2 * (        -(c1 - c3)*nd[0].z +        (c2 - c3)*nd[1].z - 2*(c1 - c2)*nd[2].z + (c1 - 2*c2 + c3)*nd[3].z + (2*c1 - c2 - c3)*nd[4].z );
    } else {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> direction must be 1 = x, 2 = y, or 3 = z\n");
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular prism integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1}
 *          c \, \deriv{\omega_i(\widehat{r})}{x,y,z} \, f(\widehat{r}) \, d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: derivative in all directions are returned here
 *  \note CJT \:: assumed top and bottom triangles have same x,y coordinates
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triPrism_dphi_f(SVECT *nd, double c, double *f, double *integral, int direction) {
#if _DEBUG
    check_prism_nodes(nd,__LINE__,__FILE__);
#endif
    
    double t1, t2;
    
    double f211, f_211, f121, f_121, f112, f_112, f332, f_332, f323, f_323, f233, f_233;
    double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21;
    
    
    f211 = 2*f[0]+   f[1]+   f[2]; f_211 = 2*f[3]+   f[4]+   f[5];
    f121 =   f[0]+ 2*f[1]+   f[2]; f_121 =   f[3]+ 2*f[4]+   f[5];
    f112 =   f[0]+   f[1]+ 2*f[2]; f_112 =   f[3]+   f[4]+ 2*f[5];
    a6  = f211 + f_211;
    a10 = f121 + f_121;
    a12 = f112 + f_112;
    
    if (direction == 3) { // dphi_i/dz
        
        t1 = one_6 * (nd[1].x*nd[0].y - nd[2].x*nd[0].y - (nd[0].x - nd[2].x)*nd[1].y + (nd[0].x - nd[1].x)*nd[2].y);
        t2 = c * t1 * one_8;
        integral[0] = t2 * a6;
        integral[1] = t2 * a10;
        integral[2] = t2 * a12;
        integral[3] = -integral[0];
        integral[4] = -integral[1];
        integral[5] = -integral[2];
        
    } else if (direction == 1 || direction == 2) { // dphi_i/dx,y
        
        
        f332 = 3*f[0]+ 3*f[1]+ 2*f[2]; f_332 = 3*f[3]+ 3*f[4]+ 2*f[5];
        f323 = 3*f[0]+ 2*f[1]+ 3*f[2]; f_323 = 3*f[3]+ 2*f[4]+ 3*f[5];
        f233 = 2*f[0]+ 3*f[1]+ 3*f[2]; f_233 = 2*f[3]+ 3*f[4]+ 3*f[5];
        
        a1  = 2*f211 + f_211; a7  = f211 + 2*f_211;
        a2  = 2*f121 + f_121; a9  = f121 + 2*f_121;
        a5  = 2*f112 + f_112; a11 = f112 + 2*f_112;
        a3  = 2*f332 + f_332; a13 = f332 + 2*f_332;
        a4  = 2*f323 + f_323; a14 = f323 + 2*f_323;
        a8  = 2*f233 + f_233; a15 = f233 + 2*f_233;
        
        a16 =        3*f[1]+   f[2]- 3*f[3]-          f[5];
        a17 =          f[1]+ 3*f[2]- 3*f[3]-   f[4];
        a18 =   f[0]+        3*f[2]-   f[3]- 3*f[4];
        a19  =  f[0]+ 3*f[1]-          f[3]-        3*f[5];
        a20 = 3*f[0]+   f[1]-                 f[4]- 3*f[5];
        a21 = 3*f[0]+          f[2]-        3*f[4]-   f[5];
        
        
        double c1, c2, c3, c4, c5, c6;
        if (direction == 2) { // y-direction
            t2 = -one_144 * c;
            c1 = nd[0].x; c2 = nd[1].x; c3 = nd[2].x;
        } else {              // x-direction
            t2 =  one_144 * c;
            c1 = nd[0].y; c2 = nd[1].y; c3 = nd[2].y;;
        }
        
        integral[0] = t2 * (-(a1*c1 + a2*c2 - a3*c3)*nd[1].z + (a1*c1 - a4*c2 + a5*c3)*nd[2].z + 3*(a6*c2 - a6*c3)*nd[3].z -
                            (a7*c1 - a2*c2 + a16*c3)*nd[4].z + (a7*c1 + a17*c2 - a5*c3)*nd[5].z);
        integral[1] = t2 * ((a1*c1 + a2*c2 - a3*c3)*nd[0].z + (a8*c1 - a2*c2 - a5*c3)*nd[2].z - (a1*c1 - a9*c2 - a21*c3)*nd[3].z -
                            3*(a10*c1 - a10*c3)*nd[4].z - (a18*c1 + a9*c2 - a5*c3)*nd[5].z);
        integral[2] = t2 * (-(a1*c1 - a4*c2 + a5*c3)*nd[0].z - (a8*c1 - a2*c2 - a5*c3)*nd[1].z + (a1*c1 - a20*c2 - a11*c3)*nd[3].z +
                            (a19*c1 - a2*c2 + a11*c3)*nd[4].z + 3*(a12*c1 - a12*c2)*nd[5].z);
        integral[3] = t2 * (-3*(a6*c2 - a6*c3)*nd[0].z + (a1*c1 - a9*c2 - a21*c3)*nd[1].z - (a1*c1 - a20*c2 - a11*c3)*nd[2].z +
                            (a7*c1 + a9*c2 - a13*c3)*nd[4].z - (a7*c1 - a14*c2 + a11*c3)*nd[5].z);
        integral[4] = t2 * ((a7*c1 - a2*c2 + a16*c3)*nd[0].z + 3*(a10*c1 - a10*c3)*nd[1].z - (a19*c1 - a2*c2 + a11*c3)*nd[2].z -
                            (a7*c1 + a9*c2 - a13*c3)*nd[3].z - (a15*c1 - a9*c2 - a11*c3)*nd[5].z);
        integral[5] = t2 * (-(a7*c1 + a17*c2 - a5*c3)*nd[0].z + (a18*c1 + a9*c2 - a5*c3)*nd[1].z - 3*(a12*c1 - a12*c2)*nd[2].z +
                            (a7*c1 - a14*c2 + a11*c3)*nd[3].z + (a15*c1 - a9*c2 - a11*c3)*nd[4].z);
    } else {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> direction must be 1 = x, 2 = y, or 3 = z\n");
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular prism integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1}
 *          c \, \deriv{\omega_i(\widehat{r})}{x,y,z} \, f(\widehat{r}) \, g(\widehat{r})
 *       \, d\widehat{z}d\widehat{y}d\widehat{x} \f$
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: derivative in all directions are returned here
 *  \note CJT \:: assumed top and bottom triangles have same x,y coordinates
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triPrism_dphi_f_g(SVECT *nd, double c, double *f, double *g, double *integral, int direction) {
#if _DEBUG
    check_prism_nodes(nd,__LINE__,__FILE__);
#endif
    
    double t1, t2;
    
    if (direction == 3) { // dphi_i/dz
        
        double f311, f_311, f221, f_221, f212, f_212, f131, f_131, f122, f_122, f113, f_113;
        double a25, a26, a27, a28, a29, a30, a42, a43, a44, a45, a46, a47, p1, p2, p3;
        
        f311 = 3*f[0] +   f[1] +   f[2]; f_311 = 3*f[3] +   f[4] +   f[5];
        f221 = 2*f[0] + 2*f[1] +   f[2]; f_221 = 2*f[3] + 2*f[4] +   f[5];
        f212 = 2*f[0] +   f[1] + 2*f[2]; f_212 = 2*f[3] +   f[4] + 2*f[5];
        f131 =   f[0] + 3*f[1] +   f[2]; f_131 =   f[3] + 3*f[4] +   f[5];
        f122 =   f[0] + 2*f[1] + 2*f[2]; f_122 =   f[3] + 2*f[4] + 2*f[5];
        f113 =   f[0] +   f[1] + 3*f[2]; f_113 =   f[3] +   f[4] + 3*f[5];
        
        a25 = 2*f311 +   f_311;
        a26 = 2*f221 +   f_221;
        a27 = 2*f212 +   f_212;
        a28 =   f311 + 2*f_311;
        a29 =   f221 + 2*f_221;
        a30 =   f212 + 2*f_212;
        a42 = 2*f131 +   f_131;
        a43 = 2*f122 +   f_122;
        a44 =   f131 + 2*f_131;
        a45 =   f122 + 2*f_122;
        a46 = 2*f113 +   f_113;
        a47 =   f113 + 2*f_113;
        
        t1 = one_6 * (nd[1].x*nd[0].y - nd[2].x*nd[0].y - (nd[0].x - nd[2].x)*nd[1].y + (nd[0].x - nd[1].x)*nd[2].y);
        t2 = c * t1 * one_120;
        
        p1 = 2*a25*g[0] +   a26*g[1] +   a27*g[2] + 2*a28*g[3] +   a29*g[4] +   a30*g[5];
        p2 =   a26*g[0] + 2*a42*g[1] +   a43*g[2] +   a29*g[3] + 2*a44*g[4] +   a45*g[5];
        p3 =   a27*g[0] +   a43*g[1] + 2*a46*g[2] +   a30*g[3] +   a45*g[4] + 2*a47*g[5];
        
        integral[0] =  t2 * p1;
        integral[1] =  t2 * p2;
        integral[2] =  t2 * p3;
        integral[3] =  -integral[0];
        integral[4] =  -integral[1];
        integral[5] =  -integral[2];
        
    } else if (direction == 1 || direction == 2) { // dphi_i/dx,y
        
        double n1, n2, n3;
        double f113, f_113, f131, f_131, f311, f_311, f122, f_122, f221, f_221, f212, f_212;
        double f334, f_334, f443, f_443, f343, f_343, f348, f_348, f834, f_834, f843, f_843;
        double f384, f_384, f483, f_483, f438, f_438, f433, f_433;
        double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20;
        double a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a40;
        double a41, a42, a43, a44, a45, a46, a47, a48, a49, a50, a51, a52, a53, a54, a55, a56, a57;
        double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15;
        double c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27, c28, c29, c30;
        double b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15;
        double d1, d2, d3, d4, d5, d6;
        double e1, e2, e3, e4, e5, e6, e7, e8, e9, e10;
        double e11, e12, e13, e14, e15, e16, e17, e18, e19, e20;
        double e21, e22, e23, e24, e25, e26, e27, e28, e29, e30, e31;
        
        // operation count = 18 + 24 + 90 = 132
        f113 =   f[0]+   f[1]+ 3*f[2]; f_113 =   f[3]+   f[4]+ 3*f[5];
        f131 =   f[0]+ 3*f[1]+   f[2]; f_131 =   f[3]+ 3*f[4]+   f[5];
        f311 = 3*f[0]+   f[1]+   f[2]; f_311 = 3*f[3]+   f[4]+   f[5];
        
        f122 =   f[0]+ 2*f[1]+ 2*f[2]; f_122 =   f[3]+ 2*f[4]+ 2*f[5];
        f221 = 2*f[0]+ 2*f[1]+   f[2]; f_221 = 2*f[3]+ 2*f[4]+   f[5];
        f212 = 2*f[0]+   f[1]+ 2*f[2]; f_212 = 2*f[3]+   f[4]+ 2*f[5];
        
        f334 = 3*f[0]+ 3*f[1]+ 4*f[2]; f_334 = 3*f[3]+ 3*f[4]+ 4*f[5];
        f433 = 4*f[0]+ 3*f[1]+ 3*f[2]; f_433 = 4*f[3]+ 3*f[4]+ 3*f[5];
        f343 = 3*f[0]+ 4*f[1]+ 3*f[2]; f_343 = 3*f[3]+ 4*f[4]+ 3*f[5];
        
        f348 = 3*f[0]+ 4*f[1]+ 8*f[2]; f_348 = 3*f[3]+ 4*f[4]+ 8*f[5];
        f834 = 8*f[0]+ 3*f[1]+ 4*f[2]; f_834 = 8*f[3]+ 3*f[4]+ 4*f[5];
        f843 = 8*f[0]+ 4*f[1]+ 3*f[2]; f_843 = 8*f[3]+ 4*f[4]+ 3*f[5];
        f384 = 3*f[0]+ 8*f[1]+ 4*f[2]; f_384 = 3*f[3]+ 8*f[4]+ 4*f[5];
        f483 = 4*f[0]+ 8*f[1]+ 3*f[2]; f_483 = 4*f[3]+ 8*f[4]+ 3*f[5];
        f438 = 4*f[0]+ 3*f[1]+ 8*f[2]; f_438 = 4*f[3]+ 3*f[4]+ 8*f[5];
        
        // operation count = 99
        a1 =  3*f311 +   f_311; a2 =  3*f221 +   f_221; a3 =  3*f212 +   f_212;
        a4 =    f311 +   f_311; a5 =    f221 +   f_221; a6 =    f212 +   f_212;
        a7 =  3*f131 +   f_131; a8 =  3*f122 +   f_122; a9 =    f131 +   f_131;
        a10 =   f122 +   f_122; a11 = 3*f843 +   f_843; a12 = 3*f483 +   f_483;
        a13 = 3*f334 +   f_334; a14 =   f843 +   f_843; a15 =   f483 +   f_483;
        a16 =   f334 +   f_334; a17 = 3*f834 +   f_834; a18 = 3*f343 +   f_343;
        a19 = 3*f438 +   f_438; a20 =   f834 +   f_834; a21 =   f343 +   f_343;
        a22 =   f438 +   f_438; a23 = 3*f113 +   f_113; a24 =   f113 +   f_113;
        a25 = 2*f311 +   f_311; a26 = 2*f221 +   f_221; a27 = 2*f212 +   f_212;
        a28 =   f311 + 2*f_311; a29 =   f221 + 2*f_221; a30 =   f212 + 2*f_212;
        a31 =   f311 + 3*f_311; a32 =   f221 + 3*f_221; a33 =   f212 + 3*f_212;
        a34 = 3*f433 +   f_433; a35 = 3*f384 +   f_384; a36 = 3*f348 +   f_348;
        a37 =   f433 +   f_433; a38 =   f384 +   f_384; a39 =   f348 +   f_348;
        a40 =   f131 + 3*f_131; a41 =   f122 + 3*f_122; a42 = 2*f131 +   f_131;
        a43 = 2*f122 +   f_122; a44 =   f131 + 2*f_131; a45 =   f122 + 2*f_122;
        a46 = 2*f113 +   f_113; a47 =   f113 + 2*f_113; a48 =   f113 + 3*f_113;
        a49 =   f843 + 3*f_843; a50 =   f483 + 3*f_483; a51 =   f334 + 3*f_334;
        a52 =   f438 + 3*f_438; a53 =   f834 + 3*f_834; a54 =   f343 + 3*f_343;
        a55 =   f384 + 3*f_384; a56 =   f348 + 3*f_348; a57 =   f433 + 3*f_433;
        
        // operation count = 200
        c1 = 4*f[1]+ f[2]- 4*f[3]- f[5];
        c2 = 4*f[0]+ 16*f[1]+ 5*f[2]+ 4*f[4]+ f[5];
        c3 = f[0]+ 5*f[1]+ 4*f[2]- f[3]+ f[4];
        c4 = 4*f[0]+ f[2]+ 16*f[3]+ 4*f[4]+ 5*f[5];
        c5 = f[0]- f[1]+ 5*f[3]+ f[4]+ 4*f[5];
        c6 = f[1]+ 4*f[2]- 4*f[3]- f[4];
        c7 = f[0]+ 4*f[1]+ 5*f[2]- f[3]+ f[5];
        c8 = 4*f[0]+ 5*f[1]+ 16*f[2]+ f[4]+ 4*f[5];
        c9 = 4*f[0]+ f[1]+ 16*f[3]+ 5*f[4]+ 4*f[5];
        c10 = f[0]- f[2]+ 5*f[3]+ 4*f[4]+ f[5];
        c11 = 16*f[0]+ 4*f[1]+ 5*f[2]+ 4*f[3]+ f[5];
        c12 = 4*f[0]+ f[2]- 4*f[4]- f[5];
        c13 = 5*f[0]+ f[1]+ 4*f[2]+ f[3]- f[4];
        c14 = 4*f[1]+ f[2]+ 4*f[3]+ 16*f[4]+ 5*f[5];
        c15 = f[0]- f[1]- f[3]- 5*f[4]- 4*f[5];
        c16 = 4*f[0]+ f[1]+ 5*f[2]- f[4]+ f[5];
        c17 = f[0]+ 4*f[2]- f[3]- 4*f[4];
        c18 = 5*f[0]+ 4*f[1]+ 16*f[2]+ f[3]+ 4*f[5];
        c19 = f[1]- f[2]+ 4*f[3]+ 5*f[4]+ f[5];
        c20 = f[0]+ 4*f[1]+ 5*f[3]+ 16*f[4]+ 4*f[5];
        c21 = 16*f[0]+ 5*f[1]+ 4*f[2]+ 4*f[3]+ f[4];
        c22 = 5*f[0]+ 4*f[1]+ f[2]+ f[3]- f[5];
        c23 = 4*f[0]+ f[1]- f[4]- 4*f[5];
        c24 = f[0]- f[2]- f[3]- 4*f[4]- 5*f[5];
        c25 = f[1]+ 4*f[2]+ 4*f[3]+ 5*f[4]+ 16*f[5];
        c26 = 4*f[0]+ 5*f[1]+ f[2]+ f[4]- f[5];
        c27 = 5*f[0]+ 16*f[1]+ 4*f[2]+ f[3]+ 4*f[4];
        c28 = f[0]+ 4*f[1]- f[3]- 4*f[5];
        c29 = f[1]- f[2]- 4*f[3]- f[4]- 5*f[5];
        c30 = f[0]+ 4*f[2]+ 5*f[3]+ 4*f[4]+ 16*f[5];
        
        // operation count = 13 X 15 = 195
        b1 =  2*a1*g[0]+    a2*g[1]+    a3*g[2]+  2*a4*g[3]+    a5*g[4]+    a6*g[5];
        b2 =    a2*g[0]+  2*a7*g[1]+    a8*g[2]+    a5*g[3]+  2*a9*g[4]+   a10*g[5];
        b3 =   a11*g[0]+   a12*g[1]+   a13*g[2]+   a14*g[3]+   a15*g[4]+   a16*g[5];
        b4 =   a17*g[0]+   a18*g[1]+   a19*g[2]+   a20*g[3]+   a21*g[4]+   a22*g[5];
        b5 =    a3*g[0]+    a8*g[1]+ 2*a23*g[2]+    a6*g[3]+   a10*g[4]+ 2*a24*g[5];
        b6 = 2*a25*g[0]+   a26*g[1]+   a27*g[2]+ 2*a28*g[3]+   a29*g[4]+   a30*g[5];
        b7 =  2*a4*g[0]+    a5*g[1]+    a6*g[2]+ 2*a31*g[3]+   a32*g[4]+   a33*g[5];
        b8 =   a34*g[0]+   a35*g[1]+   a36*g[2]+   a37*g[3]+   a38*g[4]+   a39*g[5];
        b9 =    a5*g[0]+  2*a9*g[1]+   a10*g[2]+   a32*g[3]+ 2*a40*g[4]+   a41*g[5];
        b10 =  a26*g[0]+ 2*a42*g[1]+   a43*g[2]+   a29*g[3]+ 2*a44*g[4]+   a45*g[5];
        b11 =   a6*g[0]+   a10*g[1]+ 2*a24*g[2]+   a33*g[3]+   a41*g[4]+ 2*a48*g[5];
        b12 =  a27*g[0]+   a43*g[1]+ 2*a46*g[2]+   a30*g[3]+   a45*g[4]+ 2*a47*g[5];
        b13 =  a14*g[0]+   a15*g[1]+   a16*g[2]+   a49*g[3]+   a50*g[4]+   a51*g[5];
        b14 =  a20*g[0]+   a21*g[1]+   a22*g[2]+   a53*g[3]+   a54*g[4]+   a52*g[5];
        b15 =  a37*g[0]+   a38*g[1]+   a39*g[2]+   a57*g[3]+   a55*g[4]+   a56*g[5];
        
        // operation count = 66
        d1 =  c1*g[0]+  c2*g[1]+  c3*g[2]-  c4*g[3]+  c1*g[4]-  c5*g[5];
        d2 =  c6*g[0]+  c7*g[1]+  c8*g[2]-  c9*g[3]- c10*g[4]+  c6*g[5];
        d3 = c11*g[0]+ c12*g[1]+ c13*g[2]+ c12*g[3]- c14*g[4]+ c15*g[5];
        d4 = c16*g[0]+ c17*g[1]+ c18*g[2]- c19*g[3]- c20*g[4]+ c17*g[5];
        d5 = c21*g[0]+ c22*g[1]+ c23*g[2]+ c23*g[3]+ c24*g[4]- c25*g[5];
        d6 = c26*g[0]+ c27*g[1]+ c28*g[2]+ c29*g[3]+ c28*g[4]- c30*g[5];
        
        if (direction == 2) {
            c = -c;
            n1 = nd[0].x; n2 = nd[1].x; n3 = nd[2].x;
        } else {
            n1 = nd[0].y; n2 = nd[1].y; n3 = nd[2].y;;
        }
        
        e1 = b1*n1 + b2*n2 - b3*n3;    e2 = b1*n1 - b4*n2 + b5*n3;
        e3 = b7*n1 -  b2*n2 +  d1*n3;  e4 = b7*n1 +  d2*n2 -  b5*n3;
        e5 = b8*n1 - b2*n2 - b5*n3;    e6 = b1*n1 -  d5*n2 - b11*n3;
        e7 = d6*n1 -  b2*n2 + b11*n3;  e8 = b1*n1 - b9*n2 - d3*n3;
        e9 = b7*n1 - b2*n2 + d1*n3;    e10 = b7*n1 +  b9*n2 - b13*n3;
        e11 = b7*n1 - b14*n2 + b11*n3; e12 = b15*n1 -  b9*n2 - b11*n3;
        e13 = d4*n1 + b9*n2 - b5*n3;   e14 = b1*n1 -  b9*n2 -  d3*n3;
        e15 = d4*n1 +  b9*n2 -  b5*n3; e16 = b6*n2 -  b6*n3;
        e17 = b12*n1 - b12*n2;         e18 = b6*n2 - b6*n3;
        e19 = b10*n1 - b10*n3;         e20 = b7*n1 + d2*n2 - b5*n3;
        e21 = b10*n1 - b10*n3;
        
        t2 = one_1440 * c;
        integral[0] = t2 * (               -    e1*nd[1].z +    e2*nd[2].z +  2*e16*nd[3].z -    e3*nd[4].z +    e4*nd[5].z);
        integral[1] = t2 * (    e1*nd[0].z                 +    e5*nd[2].z -    e14*nd[3].z - 2*e19*nd[4].z -   e15*nd[5].z);
        integral[2] = t2 * (   -e2*nd[0].z -    e5*nd[1].z                 +     e6*nd[3].z +    e7*nd[4].z + 2*e17*nd[5].z);
        integral[3] = t2 * (-2*e18*nd[0].z +    e8*nd[1].z -    e6*nd[2].z                  +   e10*nd[4].z -   e11*nd[5].z);
        integral[4] = t2 * (    e9*nd[0].z + 2*e21*nd[1].z -    e7*nd[2].z -    e10*nd[3].z                 -   e12*nd[5].z);
        integral[5] = t2 * (  -e20*nd[0].z +   e13*nd[1].z - 2*e17*nd[2].z +    e11*nd[3].z +   e12*nd[4].z                );
        
    } else {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> direction must be 1 = x, 2 = y, or 3 = z\n");
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following quadrilateral integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1}
 *        \big(
 *           c \, \mathbf{\nabla}\phiddd{i}(\widehat{x},\widehat{y},\widehat{z}) \cdot (\mathbf{D}(\mathbf{\nabla}f(\widehat{x},\widehat{y},\widehat{z})))
 *         \big) d\widehat{y}\,d\widehat{x}   \f$
 *  \n where D is a 3 X 3 tensor
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: Used for scalar A/D equations w/ potential anisotropic diffusion
 *  \note CJT \:: Use quadrature here, since dphi * dphi integrals will not solve analytically
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triPrism_gradPhi_dot_DgradF(SVECT *nd, SQUAD *quad, double c, STENSOR T, double *f, double *integral) {
    
#if _DEBUG
    check_prism_nodes(nd,__LINE__,__FILE__);
#endif
    
    int i, iqp, quad_order = 3; // quadrature order
    SVECT qp_grad_f;
    SQUAD_PT *qp = NULL;
    
    sarray_init_dbl(integral,NDONPRISM);
    
    for (iqp=0; iqp<quad[quad_order].n; iqp++) {
        qp = &(quad[quad_order].pt[iqp]);
        
        // evaluate triangular prism and cartesian shape function gradients at quadrature points
        qp->djac = get_triprism_linear_djac_gradPhi(qp->xhat, qp->yhat, qp->zhat, nd, qp->grad_shp);
        
        // evaluate grad(f) at quadrature points
        svect_init(&qp_grad_f);
        for (i=0; i<NDONPRISM; i++) {
            qp_grad_f.x += f[i] * qp->grad_shp[i].x;
            qp_grad_f.y += f[i] * qp->grad_shp[i].y;
            qp_grad_f.z += f[i] * qp->grad_shp[i].z;
        }
        double t1 = c * qp->djac * qp->w;
        for (i=0; i<NDONPRISM; i++) {
            integral[i] += t1 * qp->grad_shp[i].x * (T.xx * qp_grad_f.x + T.xy * qp_grad_f.y + T.xz * qp_grad_f.z);
            integral[i] += t1 * qp->grad_shp[i].y * (T.xy * qp_grad_f.x + T.yy * qp_grad_f.y + T.yz * qp_grad_f.z);
            integral[i] += t1 * qp->grad_shp[i].z * (T.xz * qp_grad_f.x + T.yz * qp_grad_f.y + T.zz * qp_grad_f.z);
        }
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following quadrilateral integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1}
 *        \big(
 *           c \, \mathbf{\nabla}\phiddd{i}(\widehat{x},\widehat{y},\widehat{z}) \cdot (f(\widehat{x},\widehat{y},\widehat{z}) \, \mathbf{v}(\widehat{x},\widehat{y},\widehat{z}))
 *         \big) d\widehat{y}\,d\widehat{x}\,d\widehat{z}   \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triPrism_gradPhi_dot_f_v(SVECT *nd, double c, double *f, SVECT *v, double *integral) {
    
#if _DEBUG
    check_prism_nodes(nd,__LINE__,__FILE__);
#endif
    
    int i;
    
    // new stuff *******************
    double ut[NDONPRISM],vt[NDONPRISM],wt[NDONPRISM],integral_X[NDONPRISM],integral_Y[NDONPRISM],integral_Z[NDONPRISM];
    for (i=0; i<NDONPRISM; i++) {
        ut[i] = v[i].x; vt[i] = v[i].y; wt[i] = v[i].z;
        integral_X[i] = 0.; integral_Y[i] = 0.; integral_Z[i] = 0.;
    }
    integrate_triPrism_dphi_f_g(nd, c, ut, f, integral_X, 1);
    integrate_triPrism_dphi_f_g(nd, c, vt, f, integral_Y, 2);
    integrate_triPrism_dphi_f_g(nd, c, wt, f, integral_Z, 3);
    for (i=0; i<NDONPRISM; i++) {integral[i] = integral_X[i] + integral_Y[i] + integral_Z[i];}
    return;
    // ********************************
    
    double t1, t2, g1, g2, g3, g4, g5, g6, n1, n2, n3;
    double f113, f_113, f131, f_131, f311, f_311, f122, f_122, f221, f_221, f212, f_212;
    double f334, f_334, f443, f_443, f343, f_343, f348, f_348, f834, f_834, f843, f_843;
    double f384, f_384, f483, f_483, f438, f_438, f433, f_433;
    double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20;
    double a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a40;
    double a41, a42, a43, a44, a45, a46, a47, a48, a49, a50, a51, a52, a53, a54, a55, a56, a57;
    double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15;
    double c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27, c28, c29, c30;
    double b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15;
    double d1, d2, d3, d4, d5, d6, p1, p2, p3;
    double e1, e2, e3, e4, e5, e6, e7, e8, e9, e10;
    double e11, e12, e13, e14, e15, e16, e17, e18, e19, e20;
    double e21, e22, e23, e24, e25, e26, e27, e28, e29, e30, e31;
    
    f113 =   f[0]+   f[1]+ 3*f[2]; f_113 =   f[3]+   f[4]+ 3*f[5];
    f131 =   f[0]+ 3*f[1]+   f[2]; f_131 =   f[3]+ 3*f[4]+   f[5];
    f311 = 3*f[0]+   f[1]+   f[2]; f_311 = 3*f[3]+   f[4]+   f[5];
    
    f122 =   f[0]+ 2*f[1]+ 2*f[2]; f_122 =   f[3]+ 2*f[4]+ 2*f[5];
    f221 = 2*f[0]+ 2*f[1]+   f[2]; f_221 = 2*f[3]+ 2*f[4]+   f[5];
    f212 = 2*f[0]+   f[1]+ 2*f[2]; f_212 = 2*f[3]+   f[4]+ 2*f[5];
    
    f334 = 3*f[0]+ 3*f[1]+ 4*f[2]; f_334 = 3*f[3]+ 3*f[4]+ 4*f[5];
    f433 = 4*f[0]+ 3*f[1]+ 3*f[2]; f_433 = 4*f[3]+ 3*f[4]+ 3*f[5];
    f343 = 3*f[0]+ 4*f[1]+ 3*f[2]; f_343 = 3*f[3]+ 4*f[4]+ 3*f[5];
    
    f348 = 3*f[0]+ 4*f[1]+ 8*f[2]; f_348 = 3*f[3]+ 4*f[4]+ 8*f[5];
    f834 = 8*f[0]+ 3*f[1]+ 4*f[2]; f_834 = 8*f[3]+ 3*f[4]+ 4*f[5];
    f843 = 8*f[0]+ 4*f[1]+ 3*f[2]; f_843 = 8*f[3]+ 4*f[4]+ 3*f[5];
    f384 = 3*f[0]+ 8*f[1]+ 4*f[2]; f_384 = 3*f[3]+ 8*f[4]+ 4*f[5];
    f483 = 4*f[0]+ 8*f[1]+ 3*f[2]; f_483 = 4*f[3]+ 8*f[4]+ 3*f[5];
    f438 = 4*f[0]+ 3*f[1]+ 8*f[2]; f_438 = 4*f[3]+ 3*f[4]+ 8*f[5];
    
    // operation count = 99
    a1 =  3*f311 +   f_311; a2 =  3*f221 +   f_221; a3 =  3*f212 +   f_212;
    a4 =    f311 +   f_311; a5 =    f221 +   f_221; a6 =    f212 +   f_212;
    a7 =  3*f131 +   f_131; a8 =  3*f122 +   f_122; a9 =    f131 +   f_131;
    a10 =   f122 +   f_122; a11 = 3*f843 +   f_843; a12 = 3*f483 +   f_483;
    a13 = 3*f334 +   f_334; a14 =   f843 +   f_843; a15 =   f483 +   f_483;
    a16 =   f334 +   f_334; a17 = 3*f834 +   f_834; a18 = 3*f343 +   f_343;
    a19 = 3*f438 +   f_438; a20 =   f834 +   f_834; a21 =   f343 +   f_343;
    a22 =   f438 +   f_438; a23 = 3*f113 +   f_113; a24 =   f113 +   f_113;
    a25 = 2*f311 +   f_311; a26 = 2*f221 +   f_221; a27 = 2*f212 +   f_212;
    a28 =   f311 + 2*f_311; a29 =   f221 + 2*f_221; a30 =   f212 + 2*f_212;
    a31 =   f311 + 3*f_311; a32 =   f221 + 3*f_221; a33 =   f212 + 3*f_212;
    a34 = 3*f433 +   f_433; a35 = 3*f384 +   f_384; a36 = 3*f348 +   f_348;
    a37 =   f433 +   f_433; a38 =   f384 +   f_384; a39 =   f348 +   f_348;
    a40 =   f131 + 3*f_131; a41 =   f122 + 3*f_122; a42 = 2*f131 +   f_131;
    a43 = 2*f122 +   f_122; a44 =   f131 + 2*f_131; a45 =   f122 + 2*f_122;
    a46 = 2*f113 +   f_113; a47 =   f113 + 2*f_113; a48 =   f113 + 3*f_113;
    a49 =   f843 + 3*f_843; a50 =   f483 + 3*f_483; a51 =   f334 + 3*f_334;
    a52 =   f438 + 3*f_438; a53 =   f834 + 3*f_834; a54 =   f343 + 3*f_343;
    a55 =   f384 + 3*f_384; a56 =   f348 + 3*f_348; a57 =   f433 + 3*f_433;
    
    // operation count = 200
    c1 = 4*f[1]+ f[2]- 4*f[3]- f[5];
    c2 = 4*f[0]+ 16*f[1]+ 5*f[2]+ 4*f[4]+ f[5];
    c3 = f[0]+ 5*f[1]+ 4*f[2]- f[3]+ f[4];
    c4 = 4*f[0]+ f[2]+ 16*f[3]+ 4*f[4]+ 5*f[5];
    c5 = f[0]- f[1]+ 5*f[3]+ f[4]+ 4*f[5];
    c6 = f[1]+ 4*f[2]- 4*f[3]- f[4];
    c7 = f[0]+ 4*f[1]+ 5*f[2]- f[3]+ f[5];
    c8 = 4*f[0]+ 5*f[1]+ 16*f[2]+ f[4]+ 4*f[5];
    c9 = 4*f[0]+ f[1]+ 16*f[3]+ 5*f[4]+ 4*f[5];
    c10 = f[0]- f[2]+ 5*f[3]+ 4*f[4]+ f[5];
    c11 = 16*f[0]+ 4*f[1]+ 5*f[2]+ 4*f[3]+ f[5];
    c12 = 4*f[0]+ f[2]- 4*f[4]- f[5];
    c13 = 5*f[0]+ f[1]+ 4*f[2]+ f[3]- f[4];
    c14 = 4*f[1]+ f[2]+ 4*f[3]+ 16*f[4]+ 5*f[5];
    c15 = f[0]- f[1]- f[3]- 5*f[4]- 4*f[5];
    c16 = 4*f[0]+ f[1]+ 5*f[2]- f[4]+ f[5];
    c17 = f[0]+ 4*f[2]- f[3]- 4*f[4];
    c18 = 5*f[0]+ 4*f[1]+ 16*f[2]+ f[3]+ 4*f[5];
    c19 = f[1]- f[2]+ 4*f[3]+ 5*f[4]+ f[5];
    c20 = f[0]+ 4*f[1]+ 5*f[3]+ 16*f[4]+ 4*f[5];
    c21 = 16*f[0]+ 5*f[1]+ 4*f[2]+ 4*f[3]+ f[4];
    c22 = 5*f[0]+ 4*f[1]+ f[2]+ f[3]- f[5];
    c23 = 4*f[0]+ f[1]- f[4]- 4*f[5];
    c24 = f[0]- f[2]- f[3]- 4*f[4]- 5*f[5];
    c25 = f[1]+ 4*f[2]+ 4*f[3]+ 5*f[4]+ 16*f[5];
    c26 = 4*f[0]+ 5*f[1]+ f[2]+ f[4]- f[5];
    c27 = 5*f[0]+ 16*f[1]+ 4*f[2]+ f[3]+ 4*f[4];
    c28 = f[0]+ 4*f[1]- f[3]- 4*f[5];
    c29 = f[1]- f[2]- 4*f[3]- f[4]- 5*f[5];
    c30 = f[0]+ 4*f[2]+ 5*f[3]+ 4*f[4]+ 16*f[5];
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // f * dphi_i/dz * w
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    g1 = v[0].z; g2 = v[1].z; g3 = v[2].z; g4 = v[3].z; g5 = v[4].z; g6 = v[5].z;
    
    a25 = 2*f311 +   f_311;
    a26 = 2*f221 +   f_221;
    a27 = 2*f212 +   f_212;
    a28 =   f311 + 2*f_311;
    a29 =   f221 + 2*f_221;
    a30 =   f212 + 2*f_212;
    a42 = 2*f131 +   f_131;
    a43 = 2*f122 +   f_122;
    a44 =   f131 + 2*f_131;
    a45 =   f122 + 2*f_122;
    a46 = 2*f113 +   f_113;
    a47 =   f113 + 2*f_113;
    
    p1 = 2*a25*g1 +   a26*g2 +   a27*g3 + 2*a28*g4 +   a29*g5 +   a30*g6;
    p2 =   a26*g1 + 2*a42*g2 +   a43*g3 +   a29*g4 + 2*a44*g5 +   a45*g6;
    p3 =   a27*g1 +   a43*g2 + 2*a46*g3 +   a30*g4 +   a45*g5 + 2*a47*g6;
    
    t1 = one_6 * (nd[1].x * nd[0].y - nd[2].x * nd[0].y - (nd[0].x - nd[2].x) * nd[1].y + (nd[0].x - nd[1].x) * nd[2].y);
    t2 = c * t1 * one_120;
    integral[0] =  t2 * p1;
    integral[1] =  t2 * p2;
    integral[2] =  t2 * p3;
    integral[3] =  -integral[0];
    integral[4] =  -integral[1];
    integral[5] =  -integral[2];
    
    for (i=0; i<2; i++) {
        if (i==0) { // f * dphi_i/dx * u ++++++++++++++++++++++++++++++++++++++++++++++++
            n1 = nd[0].y; n2 = nd[1].y; n3 = nd[2].y;
            g1 = v[0].x; g2 = v[1].x; g3 = v[2].x; g4 = v[3].x; g5 = v[4].x; g6 = v[5].x;
        } else {    // f * dphi_i/dy * v ++++++++++++++++++++++++++++++++++++++++++++++++
            c = -c;
            n1 = nd[0].x; n2 = nd[1].x; n3 = nd[2].x;
            g1 = v[0].y; g2 = v[1].y; g3 = v[2].y; g4 = v[3].y; g5 = v[4].y; g6 = v[5].y;
        }
        
        // operation count = 13 X 15 = 195
        b1 =  2*a1*g1+    a2*g2+    a3*g3+  2*a4*g4+    a5*g5+    a6*g6;
        b2 =    a2*g1+  2*a7*g2+    a8*g3+    a5*g4+  2*a9*g5+   a10*g6;
        b3 =   a11*g1+   a12*g2+   a13*g3+   a14*g4+   a15*g5+   a16*g6;
        b4 =   a17*g1+   a18*g2+   a19*g3+   a20*g4+   a21*g5+   a22*g6;
        b5 =    a3*g1+    a8*g2+ 2*a23*g3+    a6*g4+   a10*g5+ 2*a24*g6;
        b6 = 2*a25*g1+   a26*g2+   a27*g3+ 2*a28*g4+   a29*g5+   a30*g6;
        b7 =  2*a4*g1+    a5*g2+    a6*g3+ 2*a31*g4+   a32*g5+   a33*g6;
        b8 =   a34*g1+   a35*g2+   a36*g3+   a37*g4+   a38*g5+   a39*g6;
        b9 =    a5*g1+  2*a9*g2+   a10*g3+   a32*g4+ 2*a40*g5+   a41*g6;
        b10 =  a26*g1+ 2*a42*g2+   a43*g3+   a29*g4+ 2*a44*g5+   a45*g6;
        b11 =   a6*g1+   a10*g2+ 2*a24*g3+   a33*g4+   a41*g5+ 2*a48*g6;
        b12 =  a27*g1+   a43*g2+ 2*a46*g3+   a30*g4+   a45*g5+ 2*a47*g6;
        b13 =  a14*g1+   a15*g2+   a16*g3+   a49*g4+   a50*g5+   a51*g6;
        b14 =  a20*g1+   a21*g2+   a22*g3+   a53*g4+   a54*g5+   a52*g6;
        b15 =  a37*g1+   a38*g2+   a39*g3+   a57*g4+   a55*g5+   a56*g6;
        
        // operation count = 66
        d1 =  c1*g1+  c2*g2+  c3*g3-  c4*g4+  c1*g5-  c5*g6;
        d2 =  c6*g1+  c7*g2+  c8*g3-  c9*g4- c10*g5+  c6*g6;
        d3 = c11*g1+ c12*g2+ c13*g3+ c12*g4- c14*g5+ c15*g6;
        d4 = c16*g1+ c17*g2+ c18*g3- c19*g4- c20*g5+ c17*g6;
        d5 = c21*g1+ c22*g2+ c23*g3+ c23*g4+ c24*g5- c25*g6;
        d6 = c26*g1+ c27*g2+ c28*g3+ c29*g4+ c28*g5- c30*g6;
        
        e1 = b1*n1 + b2*n2 - b3*n3; e2 = b1*n1 - b4*n2 + b5*n3;
        e3 = b7*n1 -  b2*n2 +  d1*n3; e4 = b7*n1 +  d2*n2 -  b5*n3;
        e5 = b8*n1 - b2*n2 - b5*n3; e6 = b1*n1 -  d5*n2 - b11*n3;
        e7 = d6*n1 -  b2*n2 + b11*n3; e8 = b1*n1 - b9*n2 - d3*n3;
        e9 = b7*n1 - b2*n2 + d1*n3; e10 = b7*n1 +  b9*n2 - b13*n3;
        e11 = b7*n1 - b14*n2 + b11*n3; e12 = b15*n1 -  b9*n2 - b11*n3;
        e13 = d4*n1 + b9*n2 - b5*n3; e14 = b1*n1 -  b9*n2 -  d3*n3;
        e15 = d4*n1 +  b9*n2 -  b5*n3; e16 = b6*n2 -  b6*n3;
        e17 = b12*n1 - b12*n2; e18 = b6*n2 - b6*n3;
        e19 = b10*n1 - b10*n3;e20 = b7*n1 + d2*n2 - b5*n3;
        e21 = b10*n1 - b10*n3;
        
        t2 = one_1440 * c;
        integral[0] += t2 * (               -    e1*nd[1].z +    e2*nd[2].z +  2*e16*nd[3].z -    e3*nd[4].z +    e4*nd[5].z);
        integral[1] += t2 * (    e1*nd[0].z                 +    e5*nd[2].z -    e14*nd[3].z - 2*e19*nd[4].z -   e15*nd[5].z);
        integral[2] += t2 * (   -e2*nd[0].z -    e5*nd[1].z                 +     e6*nd[3].z +    e7*nd[4].z + 2*e17*nd[5].z);
        integral[3] += t2 * (-2*e18*nd[0].z +    e8*nd[1].z -    e6*nd[2].z                  +   e10*nd[4].z -   e11*nd[5].z);
        integral[4] += t2 * (    e9*nd[0].z + 2*e21*nd[1].z -    e7*nd[2].z -    e10*nd[3].z                 -   e12*nd[5].z);
        integral[5] += t2 * (  -e20*nd[0].z +   e13*nd[1].z - 2*e17*nd[2].z +    e11*nd[3].z +   e12*nd[4].z                );
        
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following quadrilateral integration:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1}
 *        \big(
 *           c \, \mathbf{\nabla}\phiddd{i}(\widehat{x},\widehat{y},\widehat{z}) \cdot (f(\widehat{x},\widehat{y},\widehat{z}) \, \mathbf{v}(\widehat{x},\widehat{y},\widehat{z}))
 *         \big) d\widehat{y}\,d\widehat{x}\,d\widehat{z}   \f$
 *  \author  Corey Trahan, Ph.D.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triPrism_gradPhi_dot_v(SVECT *nd, double c, SVECT *v, double *integral) {
    
#if _DEBUG
    check_prism_nodes(nd,__LINE__,__FILE__);
#endif
    
    // new stuff *******************
    int i;
    double ut[NDONPRISM],vt[NDONPRISM],wt[NDONPRISM],integral_X[NDONPRISM],integral_Y[NDONPRISM],integral_Z[NDONPRISM];
    for (i=0; i<NDONPRISM; i++) {
        ut[i] = v[i].x; vt[i] = v[i].y; wt[i] = v[i].z;
        integral_X[i] = 0.; integral_Y[i] = 0.; integral_Z[i] = 0.;
    }
    integrate_triPrism_dphi_f(nd, c, ut, integral_X, 1);
    integrate_triPrism_dphi_f(nd, c, vt, integral_Y, 2);
    integrate_triPrism_dphi_f(nd, c, wt, integral_Z, 3);
    for (i=0; i<NDONPRISM; i++) integral[i] = integral_X[i] + integral_Y[i] + integral_Z[i];
    return;
    // ********************************
    
    double t1, x12, x13, x23, d1, d2, d3, d4, d5, d6, d7, d8, e1, e2, e3, e4, e5, e6, e7, f1, f2, f3, f4, f5, f6, g1, g2, g3;
    double u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6, w1, w2, w3, w4, w5, w6;
    double w211211, w121121, w112112;
    double v422211, v242121, v664332, v646323, v224112, v211211, v211422, v466233, v121242, v121121, v112224, v112112, v332664, v323646, v233466;
    double u422211, u242121, u664332, u646323, u224112, u211211, u211422, u466233, u121242, u121121, u112224, u112112, u332664, u323646, u233466;
    double v031_301, v013_310, v301_031, v103_130, v310_013, v130_103;
    double u031_301, u013_310, u301_031, u103_130, u310_013, u130_103;
    
    double x1 = nd[0].x, x2 = nd[1].x, x3 = nd[2].x;
    double y1 = nd[0].y, y2 = nd[1].y, y3 = nd[2].y;
    double z1 = nd[0].z, z2 = nd[1].z, z3 = nd[2].z, z4 = nd[3].z, z5 = nd[4].z, z6 = nd[5].z;
    
    u1 = v[0].x; u2 = v[1].x; u3 = v[2].x; u4 = v[3].x; u5 = v[4].x; u6 = v[5].x;
    v1 = v[0].y; v2 = v[1].y; v3 = v[2].y; v4 = v[3].y; v5 = v[4].y; v6 = v[5].y;
    w1 = v[0].z; w2 = v[1].z; w3 = v[2].z; w4 = v[3].z; w5 = v[4].z; w6 = v[5].z;
    
    x12 = (nd[0].x - nd[1].x); x23 = (nd[1].x - nd[2].x); x13 = (nd[0].x - nd[2].x);
    
    // can be broke down further ... I'm just going fucking insane doing this ...
    w211211 = (2*w1 + w2 + w3 + 2*w4 + w5 + w6); w121121 = (w1 + 2*w2 + w3 + w4 + 2*w5 + w6);
    w112112 =(w1 + w2 + 2*w3 + w4 + w5 + 2*w6);
    
    v422211 = (4*v1 + 2*v2 + 2*v3 + 2*v4 + v5 + v6); v242121 = (2*v1 + 4*v2 + 2*v3 + v4 + 2*v5 + v6);
    v664332 = (6*v1 + 6*v2 + 4*v3 + 3*v4 + 3*v5 + 2*v6); v646323 = (6*v1 + 4*v2 + 6*v3 + 3*v4 + 2*v5 + 3*v6);
    v224112 = (2*v1 + 2*v2 + 4*v3 + v4 + v5 + 2*v6); v211211 = (2*v1 + v2 + v3 + 2*v4 + v5 + v6);
    v211422 = (2*v1 + v2 + v3 + 4*v4 + 2*v5 + 2*v6); v466233 = (4*v1 + 6*v2 + 6*v3 + 2*v4 + 3*v5 + 3*v6);
    v121242 = (v1 + 2*v2 + v3 + 2*v4 + 4*v5 + 2*v6); v121121 = (v1 + 2*v2 + v3 + v4 + 2*v5 + v6);
    v112224 = (v1 + v2 + 2*v3 + 2*v4 + 2*v5 + 4*v6); v112112 = (v1 + v2 + 2*v3 + v4 + v5 + 2*v6);
    v332664 = (3*v1 + 3*v2 + 2*v3 + 6*v4 + 6*v5 + 4*v6); v323646 = (3*v1 + 2*v2 + 3*v3 + 6*v4 + 4*v5 + 6*v6);
    v233466 = (2*v1 + 3*v2 + 3*v3 + 4*v4 + 6*v5 + 6*v6);
    
    u422211 = (4*u1 + 2*u2 + 2*u3 + 2*u4 + u5 + u6); u242121 = (2*u1 + 4*u2 + 2*u3 + u4 + 2*u5 + u6);
    u664332 = (6*u1 + 6*u2 + 4*u3 + 3*u4 + 3*u5 + 2*u6); u646323 = (6*u1 + 4*u2 + 6*u3 + 3*u4 + 2*u5 + 3*u6);
    u224112 = (2*u1 + 2*u2 + 4*u3 + u4 + u5 + 2*u6); u211211 = (2*u1 + u2 + u3 + 2*u4 + u5 + u6);
    u211422 = (2*u1 + u2 + u3 + 4*u4 + 2*u5 + 2*u6); u466233 = (4*u1 + 6*u2 + 6*u3 + 2*u4 + 3*u5 + 3*u6);
    u121242 = (u1 + 2*u2 + u3 + 2*u4 + 4*u5 + 2*u6); u121121 = (u1 + 2*u2 + u3 + u4 + 2*u5 + u6);
    u112224 = (u1 + u2 + 2*u3 + 2*u4 + 2*u5 + 4*u6); u112112 = (u1 + u2 + 2*u3 + u4 + u5 + 2*u6);
    u332664 = (3*u1 + 3*u2 + 2*u3 + 6*u4 + 6*u5 + 4*u6); u323646 = (3*u1 + 2*u2 + 3*u3 + 6*u4 + 4*u5 + 6*u6);
    u233466 = (2*u1 + 3*u2 + 3*u3 + 4*u4 + 6*u5 + 6*u6);
    
    v031_301 = (3*v2 + v3 - 3*v4 - v6); u031_301 = (3*u2 + u3 - 3*u4 - u6);
    v013_310 = (v2 + 3*v3 - 3*v4 - v5); u013_310 = (u2 + 3*u3 - 3*u4 - u5);
    v301_031 = (3*v1 + v3 - 3*v5 - v6); u301_031 = (3*u1 + u3 - 3*u5 - u6);
    v103_130 = (v1 + 3*v3 - v4 - 3*v5); u103_130 = (u1 + 3*u3 - u4 - 3*u5);
    v310_013 = (3*v1 + v2 - v5 - 3*v6); u310_013 = (3*u1 + u2 - u5 - 3*u6);
    v130_103 = (v1 + 3*v2 - v4 - 3*v6); u130_103 = (u1 + 3*u2 - u4 - 3*u6);
    
    d1 = w211211*x23; d2 = w211211*x13; d3 = w211211*x12;
    d4 = v422211*x1 + v242121*x2 - v664332*x3 - u422211*y1 - u242121*y2 + u664332*y3;
    d5 = v422211*x1 - v646323*x2 + v224112*x3 - u422211*y1 + u646323*y2 - u224112*y3;
    d6 = v211211*x2 - v211211*x3 - u211211*y2 + u211211*y3;
    d7 = v211422*x1 - v242121*x2 + v031_301*x3 - u211422*y1 + u242121*y2 - u031_301*y3;
    d8 = v211422*x1 + v013_310*x2 - v224112*x3 - u211422*y1 - u013_310*y2 + u224112*y3;
    
    e1 = w121121*x23; e2 = w121121*x13; e3 = w121121*x12;
    e4 = v466233*x1 - v242121*x2 - v224112*x3 - u466233*y1 + u242121*y2 + u224112*y3;
    e5 = v422211*x1 - v121242*x2 - v301_031*x3 - u422211*y1 + u121242*y2 + u301_031*y3;
    e6 = v121121*x1 - v121121*x3 - u121121*y1 + u121121*y3;
    e7 = v103_130*x1 + v121242*x2 - v224112*x3 - u103_130*y1 - u121242*y2 + u224112*y3;
    
    f1 = w112112*x23; f2 = w112112*x13; f3 = w112112*x12;
    f4 = v422211*x1 - v310_013*x2 - v112224*x3 - u422211*y1 + u310_013*y2 + u112224*y3;
    f5 = v130_103*x1 - v242121*x2 + v112224*x3 - u130_103*y1 + u242121*y2 - u112224*y3;
    f6 = v112112*x1 - v112112*x2 - u112112*y1 + u112112*y2;
    
    g1 = v211422*x1 + v121242*x2 - v332664*x3 - u211422*y1 - u121242*y2 + u332664*y3;
    g2 = v211422*x1 - v323646*x2 + v112224*x3 - u211422*y1 + u323646*y2 - u112224*y3;
    g3 = v233466*x1 - v121242*x2 - v112224*x3 - u233466*y1 + u121242*y2 + u112224*y3;

    t1 = one_144 * c;
    integral[0] = t1 * (  3*( d1*y1 - d2*y2 + d3*y3) +   d4*z2 -   d5*z3 - 3*d6*z4 +   d7*z5 -   d8*z6 );
    integral[1] = t1 * (  3*( e1*y1 - e2*y2 + e3*y3) -   d4*z1 -   e4*z3 +   e5*z4 + 3*e6*z5 +   e7*z6 );
    integral[2] = t1 * (  3*( f1*y1 - f2*y2 + f3*y3) +   d5*z1 +   e4*z2 -   f4*z4 -   f5*z5 - 3*f6*z6 );
    integral[3] = t1 * (  3*(-d1*y1 + d2*y2 - d3*y3) + 3*d6*z1 -   e5*z2 +   f4*z3 -   g1*z5 +   g2*z6 );
    integral[4] = t1 * (  3*(-e1*y1 + e2*y2 - e3*y3) -   d7*z1 - 3*e6*z2 +   f5*z3 +   g1*z4 +   g3*z6 );
    integral[5] = t1 * (  3*(-f1*y1 + f2*y2 - f3*y3) +   d8*z1 -   e7*z2 + 3*f6*z3 -   g2*z4 -   g3*z5 );
}
    
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular prism integration using quadrature:
 *  \f{eqnarray*}{
 *       integral_x &=&
 *        \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1}
 *          c \, (\mathbf{\nabla} \omega_i(\widehat{x},\widehat{y},\widehat{z}) \cdot T_x \, \mathbf{\nabla} u(\widehat{x},\widehat{y},\widehat{z}))
 *       \, d\widehat{z}d\widehat{y}d\widehat{x}  \\
 *       integral_y &=&
 *        \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1}
 *          c \, (\mathbf{\nabla} \omega_i(\widehat{x},\widehat{y},\widehat{z}) \cdot T_y \,\mathbf{\nabla} v(\widehat{x},\widehat{y},\widehat{z}))
 *       \, d\widehat{z}d\widehat{y}d\widehat{x}
 *  \f}
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: uses quadrature, cannot directly integrate in sage
 *  \note CJT \:: assumed top and bottom triangles have same x,y coordinates
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triPrism_gradPhi_dot_DgradV(SVECT *nd, SQUAD *quad, double c, SVECT tx, SVECT ty, double *u, double *v, double *integral_x, double *integral_y) {
    
#if _DEBUG
    check_prism_nodes(nd,__LINE__,__FILE__);
#endif
    
    int i, iqp, quad_order = 4; // quadrature order (cjt :: I think this can be order 3)
    SQUAD_PT *qp = NULL;
    SVECT qp_grad_u, qp_grad_v;
    
    sarray_init_dbl(integral_x,NDONPRISM);
    sarray_init_dbl(integral_y,NDONPRISM);
    for (iqp=0; iqp<quad[quad_order].n; iqp++) {
        qp = &(quad[quad_order].pt[iqp]);
        
        // evaluate triangular prism and cartesian shape function gradients at quadrature points
        qp->djac = get_triprism_linear_djac_gradPhi(qp->xhat, qp->yhat, qp->zhat, nd, qp->grad_shp);
        
        // evalute grad(vel) at quadrature points
        svect_init(&qp_grad_u); svect_init(&qp_grad_v);
        for (i=0; i<NDONPRISM; i++) {
            qp_grad_u.x += u[i] * qp->grad_shp[i].x;
            qp_grad_u.y += u[i] * qp->grad_shp[i].y;
            qp_grad_u.z += u[i] * qp->grad_shp[i].z;
            qp_grad_v.x += v[i] * qp->grad_shp[i].x;
            qp_grad_v.y += v[i] * qp->grad_shp[i].y;
            qp_grad_v.z += v[i] * qp->grad_shp[i].z;
        }
        qp_grad_u.x *= tx.x;  qp_grad_u.y *= tx.y;  qp_grad_u.z *= tx.z;
        qp_grad_v.x *= ty.x;  qp_grad_v.y *= ty.y;  qp_grad_v.z *= ty.z;
        
        double t1 = c * qp->djac * qp->w;
        for (i=0; i<NDONPRISM; i++) {
            integral_x[i] += t1 * svect_dotp(qp->grad_shp[i], qp_grad_u);
            integral_y[i] += t1 * svect_dotp(qp->grad_shp[i], qp_grad_v);
        }
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief  Peforms following triangular prism integration using quadrature:
 *  \f{eqnarray*}{
 *       integral_x &=& \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1} c \, \mathbf{\nabla} \deriv{\phiddd{i}}{x} \, f(\widehat{x},\widehat{y},\widehat{z}) \, d\widehat{z}d\widehat{y}d\widehat{x}  \\
 *       integral_y &=& \int_{0}^{1} \int_{0}^{1-\widehat{x}} \int_{-1}^{1} c \, \mathbf{\nabla} \deriv{\phiddd{i}}{y} \, f(\widehat{x},\widehat{y},\widehat{z}) \, d\widehat{z}d\widehat{y}d\widehat{x}
 *  \f}
 *  \author  Corey Trahan, Ph.D.
 *  \note CJT \:: assumed top and bottom triangles have same x,y coordinates
 *  \note CJT \:: total operation count = about 500 * 8 quadrature evaluations = 4000
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
inline void integrate_triPrism_dphi_f_quadrature(SVECT *nd, SQUAD *quad, double c, double *f, double *integral_x, double *integral_y, int quad_order) {
    
#if _DEBUG
    check_prism_nodes(nd,__LINE__,__FILE__);
#endif
    
    int i, iqp;
    double qp_f = 0.;
    SQUAD_PT *qp = NULL;
    
    sarray_init_dbl(integral_x,NDONPRISM);
    sarray_init_dbl(integral_y,NDONPRISM);
    for (iqp=0; iqp<quad[quad_order].n; iqp++) {
        qp = &(quad[quad_order].pt[iqp]);
        
        // evaluate triangular prism and cartesian shape function gradients at quadrature points
        qp->djac = get_triprism_linear_djac_gradPhi(qp->xhat, qp->yhat, qp->zhat, nd, qp->grad_shp);
        
        // evaluate function at quadrature point :: 14 operations
        qp_f = SQUAD_get_function(qp, f, NDONPRISM);
        
        // 16 operations
        double t1 = c * qp_f * qp->djac * qp->w;
        for (i=0; i<NDONPRISM; i++) {
            integral_x[i] += t1 * qp->grad_shp[i].x;
            integral_y[i] += t1 * qp->grad_shp[i].y;
        }
    }
}










