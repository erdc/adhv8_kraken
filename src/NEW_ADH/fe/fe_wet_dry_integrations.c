/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_wet_dry_integrations.c This file collects the 2D shallow water wet/dry integrations.
 *   All function arguements must be consistent with fe_sw2_wet_dry_wrapper.                */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Wet/Dry routine for calculating the wet-dry integral of shallow water a generic function, f, and/or vector, v:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \big( c \, f(\widehat{x},\widehat{y})\big) d\widehat{y}\,d\widehat{x} \f$
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gary Brown
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \date      October, 2017
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[inout] rhs the 2D elemental residual array contribution with wet/dry inclusion
 *  @param[in]  elem_nds not used
 *  @param[in]  elem_head the 2d element nodal wet-dry depths
 *  @param[in]  grad_phi not used
 *  @param[in]  v not used
 *  @param[in]  v_wet_dry wet-dry velocities (averaged here)
 *  @param[in]  f not used
 *  @param[in]  f_wet_dry the scalar function to be averaged (averaged here)
 *  @param[in]  djac the 2d element area
 *  @param[in]  dt not used
 *  @param[in]  c a constant
 *
 * \note CJT \:: must have same call arguements as all generic wet-dry wrapped functions
 * \note CJT \:: djac is used here because later, in the wet-dry wrapper, it is divided by the wet area (djac)
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void inline fe_sw2_wd_average(SVECT *elem_nds, double *depth, SVECT2D *v_wd, double *f_wd, double djac, double *f_wd_avg, SVECT2D *v_wd_avg) {
    
    double averages[NDONTRI*3]; sarray_init_dbl(averages,NDONTRI*3);
    int redistribute_flag = OFF; // do not redistribute in SUPG
    double vars[2];
    vars[0] = 1.;
    vars[1] = 1.;
    double factor = fe_sw2_wet_dry_wrapper(averages, elem_nds, depth, NULL, NULL, v_wd, NULL, f_wd, djac, redistribute_flag, OFF, vars, fe_sw2_wd_integrate_triangle_f);
    //double factor = fe_sw2_wet_dry_wrapper(averages, elem_nds, depth, NULL, NULL, v_wd, NULL, f_wd, djac, 1., 1., redistribute_flag, OFF, fe_sw2_wd_average_tri);
    *f_wd_avg = 0.;
    if (v_wd_avg != NULL) {v_wd_avg->x = 0.; v_wd_avg->y = 0.;}
    if (factor > 0.0) {
        double wet_area = djac * factor;
        *f_wd_avg   = averages[0] / wet_area;
        if (v_wd_avg != NULL) {
            v_wd_avg->x = averages[1] / wet_area; // cjt :: note, same as (*v_avg).x
            v_wd_avg->y = averages[2] / wet_area; // cjt :: note, same as (*v_avg).y
        }
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Add the 2D shallow water wet-dry wrapped convection terms to the FE residual.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[inout] rhs the 2D elemental residual array contribution with wet/dry inclusion
 *  @param[in]  elem_nds not used
 *  @param[in]  elem_head the 2d element nodal wet-dry depths
 *  @param[in]  grad_phi the elemental shape function gradients in configuration space
 *  @param[in]  v fully wet velocites
 *  @param[in]  v_wet_dry the wet-dry velocities
 *  @param[in]  f not used
 *  @param[in]  f_wet_dry not used
 *  @param[in]  djac the 2d element area
 *  @param[in]  dt not used
 *  @param[in]  c a constant
 *
 *  \details Solves the following weak, wet/dry discrete body terms of the 2D shallow water equation: \n
 *  \f{eqnarray*}{
 *   \residDA{i}{}{c}   &=& -  \bodyConv{\,2d}{e}{\phidd{i}}{(\velb{h} \, \depth{h})} \\
 *   \residDA{i}{}{mx}  &=& -  \bodyConv{\,2d}{e}{\phidd{i}}{(\velb{h} \, \ub{h}\depth{h})} \\
 *   \residDA{i}{}{my}  &=& -  \bodyConv{\,2d}{e}{\phidd{i}}{(\velb{h} \, \vb{h}\depth{h})}
 *  \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void inline fe_sw2_wd_convection_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *elem_rhs) {
    
#ifdef _DEBUG
    assert(f == NULL);
    assert(f_wet_dry == NULL);
#endif
    
    double dt = vars[0];
    double c = vars[1];
    
    int i;
    c *= dt;
    
    double rhs_c[NDONTRI] = {0., 0., 0.};
    double rhs_x[NDONTRI] = {0., 0., 0.};
    double rhs_y[NDONTRI] = {0., 0., 0.};
    
    double uu[NDONTRI], vv[NDONTRI];
    dumpVector2D(v_wet_dry, NDONTRI, uu, vv);
    
    // mass contibution :: d(phi)/dx uh + d(phi)/dy vh
    integrate_triangle_gradPhi_dot_f_v(grad_phi, djac, c, elem_head, v_wet_dry, rhs_c);
    
    // momentum contribution :: dphi/dx * uvh, dphi/dy * uvh
    integrate_triangle_dphi_f_g_h(grad_phi, djac, c, uu, vv, elem_head, rhs_y, rhs_x);
    
    // momentum contribution :: dphi/dx * uuh, dphi/dy * vvh
    integrate_triangle_dphi_f_f_h(grad_phi, djac, c, uu, elem_head, rhs_x, NULL); // integer here is direction of derivative
    integrate_triangle_dphi_f_f_h(grad_phi, djac, c, vv, elem_head, NULL, rhs_y); // integer here is direction of derivative
    
    for (i=0; i<NDONTRI; i++) {
        elem_rhs[i*3] += rhs_c[i];
        elem_rhs[i*3+1] += rhs_x[i];
        elem_rhs[i*3+2] += rhs_y[i];
    }
    
    //    // REMOVE LATER ***********
    //    // old adh way
    //    c = -dt;
    //    double const_contrib;     /* the contributions */
    //    double sum, sum_u, sum_v, sum_h;
    //
    //    sum_u = (v_wet_dry[0].x + v_wet_dry[1].x + v_wet_dry[2].x);
    //    sum_v = (v_wet_dry[0].y + v_wet_dry[1].y + v_wet_dry[2].y);
    //    sum_h = (elem_head[0] + elem_head[1] + elem_head[2]);
    //
    //    /* contributions to the continuity equation */
    //    /* d(phi)/dx uh + d(phi)/dy vh  */
    //    const_contrib = djac * c / 12.;
    //
    //    sum = sum_u * sum_h + elem_head[0] * v_wet_dry[0].x + elem_head[1] * v_wet_dry[1].x + elem_head[2] * v_wet_dry[2].x;
    //    elem_rhs[0].c_eq = const_contrib * sum * grad_phi[0].x;
    //    elem_rhs[1].c_eq = const_contrib * sum * grad_phi[1].x;
    //    elem_rhs[2].c_eq = const_contrib * sum * grad_phi[2].x;
    //
    //    sum = sum_v * sum_h + elem_head[0] * v_wet_dry[0].y + elem_head[1] * v_wet_dry[1].y + elem_head[2] * v_wet_dry[2].y;
    //    elem_rhs[0].c_eq += const_contrib * sum * grad_phi[0].y;
    //    elem_rhs[1].c_eq += const_contrib * sum * grad_phi[1].y;
    //    elem_rhs[2].c_eq += const_contrib * sum * grad_phi[2].y;
    //    // ************************
    //
    //
    //
    //    /* the contribution to the y equation of */
    //    /* d(phi)/dx uvh  */
    //    sum =
    //    elem_head[0] * (2. * v_wet_dry[0].x * (3. * v_wet_dry[0].y + v_wet_dry[1].y + v_wet_dry[2].y) +
    //            v_wet_dry[1].x * (2. * v_wet_dry[0].y + 2. * v_wet_dry[1].y + v_wet_dry[2].y) + v_wet_dry[2].x * (2. * v_wet_dry[0].y + v_wet_dry[1].y +
    //                                                                      2. * v_wet_dry[2].y)) +
    //    elem_head[2] * (v_wet_dry[0].x * (2. * v_wet_dry[0].y + v_wet_dry[1].y + 2. * v_wet_dry[2].y) +
    //            2. * v_wet_dry[2].x * (v_wet_dry[0].y + v_wet_dry[1].y + 3. * v_wet_dry[2].y) + v_wet_dry[1].x * (v_wet_dry[0].y +
    //                                                                      2. * (v_wet_dry[1].y +
    //                                                                            v_wet_dry[2].y))) +
    //    elem_head[1] * (v_wet_dry[0].x * (2. * v_wet_dry[0].y + 2. * v_wet_dry[1].y + v_wet_dry[2].y) +
    //            2. * v_wet_dry[1].x * (v_wet_dry[0].y + 3. * v_wet_dry[1].y + v_wet_dry[2].y) + v_wet_dry[2].x * (v_wet_dry[0].y +
    //                                                                      2. * (v_wet_dry[1].y +
    //                                                                            v_wet_dry[2].y)));
    //    const_contrib = sum * djac * c / 60.;
    //
    //    elem_rhs[0].y_eq += const_contrib * grad_phi[0].x;
    //    elem_rhs[1].y_eq += const_contrib * grad_phi[1].x;
    //    elem_rhs[2].y_eq += const_contrib * grad_phi[2].x;
    //
    //    /* the contribution to the x equation of */
    //    /* d(phi)/dy uvh */
    //    elem_rhs[0].x_eq += const_contrib * grad_phi[0].y;
    //    elem_rhs[1].x_eq += const_contrib * grad_phi[1].y;
    //    elem_rhs[2].x_eq += const_contrib * grad_phi[2].y;
    //
    //    /* the contribution to the x equation of */
    //    /* d(phi)/dx uuh */
    //    sum =
    //    elem_head[1] * (v_wet_dry[0].x * v_wet_dry[0].x + 2. * v_wet_dry[0].x * v_wet_dry[1].x + 3. * v_wet_dry[1].x * v_wet_dry[1].x +
    //            v_wet_dry[0].x * v_wet_dry[2].x + 2. * v_wet_dry[1].x * v_wet_dry[2].x + v_wet_dry[2].x * v_wet_dry[2].x) +
    //    elem_head[2] * (v_wet_dry[0].x * v_wet_dry[0].x + v_wet_dry[0].x * v_wet_dry[1].x + v_wet_dry[1].x * v_wet_dry[1].x + 2. * v_wet_dry[0].x * v_wet_dry[2].x +
    //            2. * v_wet_dry[1].x * v_wet_dry[2].x + 3. * v_wet_dry[2].x * v_wet_dry[2].x) + elem_head[0] * (3. * v_wet_dry[0].x * v_wet_dry[0].x +
    //                                                                   v_wet_dry[1].x * v_wet_dry[1].x +
    //                                                                   v_wet_dry[1].x * v_wet_dry[2].x +
    //                                                                   v_wet_dry[2].x * v_wet_dry[2].x +
    //                                                                   2. * v_wet_dry[0].x * (v_wet_dry[1].x +
    //                                                                                  v_wet_dry[2].x));
    //    const_contrib = sum * djac * c / 30.;
    //
    //    elem_rhs[0].x_eq += const_contrib * grad_phi[0].x;
    //    elem_rhs[1].x_eq += const_contrib * grad_phi[1].x;
    //    elem_rhs[2].x_eq += const_contrib * grad_phi[2].x;
    //
    //    /* the contribution to the y equation of */
    //    /* d(phi)/dy vvh */
    //    sum =
    //    elem_head[1] * (v_wet_dry[0].y * v_wet_dry[0].y + 2. * v_wet_dry[0].y * v_wet_dry[1].y + 3. * v_wet_dry[1].y * v_wet_dry[1].y +
    //            v_wet_dry[0].y * v_wet_dry[2].y + 2. * v_wet_dry[1].y * v_wet_dry[2].y + v_wet_dry[2].y * v_wet_dry[2].y) +
    //    elem_head[2] * (v_wet_dry[0].y * v_wet_dry[0].y + v_wet_dry[0].y * v_wet_dry[1].y + v_wet_dry[1].y * v_wet_dry[1].y + 2. * v_wet_dry[0].y * v_wet_dry[2].y +
    //            2. * v_wet_dry[1].y * v_wet_dry[2].y + 3. * v_wet_dry[2].y * v_wet_dry[2].y) + elem_head[0] * (3. * v_wet_dry[0].y * v_wet_dry[0].y +
    //                                                                   v_wet_dry[1].y * v_wet_dry[1].y +
    //                                                                   v_wet_dry[1].y * v_wet_dry[2].y +
    //                                                                   v_wet_dry[2].y * v_wet_dry[2].y +
    //                                                                   2. * v_wet_dry[0].y * (v_wet_dry[1].y +
    //                                                                                  v_wet_dry[2].y));
    //
    //    const_contrib = sum * djac * c / 30.;
    //
    //    elem_rhs[0].y_eq += const_contrib * grad_phi[0].y;
    //    elem_rhs[1].y_eq += const_contrib * grad_phi[1].y;
    //    elem_rhs[2].y_eq += const_contrib * grad_phi[2].y;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \brief     Routine for wet-dry integrating the SW 2D continuity temporal term on a triangle.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[inout] rhs the 2D elemental residual array contribution with wet/dry inclusion
 *  @param[in]  elem_nds not used
 *  @param[in]  elem_head the 2d element nodal wet-dry depths
 *  @param[in]  grad_phi not used
 *  @param[in]  v not used
 *  @param[in]  v_wet_dry not used
 *  @param[in]  f not used
 *  @param[in]  f_wet_dry not used
 *  @param[in]  djac the 2d element area
 *  @param[in]  dt not used
 *  @param[in]  c the second order time-derivative weight
 *
 *  \f{eqnarray*}{ \residDA{i}{}{c} &=& dt * \alpha * \bodyTime{\,2d}{e}{\phidd{i}}{\depth{h}}  \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void inline fe_sw2_wd_continuity_temporal_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *rhs) {
    
    double dt = vars[0];
    double c = vars[1];
    
    int i;
    double rhs_c_eq[NDONTRI]; sarray_init_dbl(rhs_c_eq, NDONTRI);
    integrate_triangle_phi_f(djac, c, elem_head, rhs_c_eq);
    for (i=0; i<NDONTRI; i++) {
        rhs[i*3] += rhs_c_eq[i];;
    }
    
    // old adh way
    //    double uu[3],vv[3];
    //    uu[0] = 1.; uu[1] = 1.; uu[2] = 1.;
    //    vv[0] = elem_head[0]; vv[1] = elem_head[1]; vv[2] = elem_head[2];
    //    double contrib = djac * c * one_60;
    //    double a = 2. * (vv[0] + vv[1] + vv[2]);
    //     rhs[0].c_eq += contrib * (uu[0] * (a + 4. * vv[0]) + uu[1] * (a - vv[2]) + uu[2] * (a - vv[1]));
    //     rhs[1].c_eq += contrib * (uu[0] * (a - vv[2]) + uu[1] * (a + 4. * vv[1]) + uu[2] * (a - vv[0]));
    //     rhs[2].c_eq += contrib * (uu[0] * (a - vv[1]) + uu[1] * (a - vv[0]) + uu[2] * (a + 4. * vv[2]));
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Add the 2D shallow water wet-dry wrapped pressure terms to the FE residual.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[inout] rhs the 2D elemental residual array contribution with wet/dry inclusion
 *  @param[in]  elem_nds not used
 *  @param[in]  elem_head the 2d element nodal wet-dry depths
 *  @param[in]  grad_phi the elemental shape function gradients in configuration space
 *  @param[in]  v not used
 *  @param[in]  v_wet_dry not used
 *  @param[in]  f not used
 *  @param[in]  f_wet_dry not used
 *  @param[in]  djac the 2d element area
 *  @param[in]  dt not used
 *  @param[in]  c a constant
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void inline fe_sw2_wd_pressure_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *rhs) {
    
    double dt = vars[0];
    double c = vars[1];
    
    int i;
    double rhs_x[NDONTRI] = {0., 0., 0.};
    double rhs_y[NDONTRI] = {0., 0., 0.};
    integrate_triangle_dphi_f_f(grad_phi, djac, c, elem_head, rhs_x, rhs_y);
    for (i=0; i<NDONTRI; i++) {
        rhs[i*3+1] += rhs_x[i];
        rhs[i*3+2] += rhs_y[i];
    }
    
    
    // old adh way
    //    double sum_h = (elem_head[0] + elem_head[1] + elem_head[2]);
    //    double sum = sum_h * sum_h - elem_head[0] * elem_head[1] - elem_head[1] * elem_head[2] - elem_head[2] * elem_head[0];
    //    double const_contrib = sum * djac * c / 6.;
    //
    //    /* contribution to the x eq. of d(phi)/dx 1/2 g h h */
    //    rhs[0].x_eq += const_contrib * grad_phi[0].x;
    //    rhs[1].x_eq += const_contrib * grad_phi[1].x;
    //    rhs[2].x_eq += const_contrib * grad_phi[2].x;
    //
    //    /* contribution to the y eq. of d(phi)/dy 1/2 g h h */
    //    rhs[0].y_eq += const_contrib * grad_phi[0].y;
    //    rhs[1].y_eq += const_contrib * grad_phi[1].y;
    //    rhs[2].y_eq += const_contrib * grad_phi[2].y;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Add the 2D shallow water wet-dry wrapped body-force terms to the FE residual.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[inout] rhs the 2D elemental residual array contribution with wet/dry inclusion
 *  @param[in]  elem_nds not used
 *  @param[in]  elem_head the 2d element nodal wet-dry depths
 *  @param[in]  grad_phi not used
 *  @param[in]  v the fully wet bathymetry gradient
 *  @param[in]  v_wet_dry not used
 *  @param[in]  f not used
 *  @param[in]  f_wet_dry not used
 *  @param[in]  djac the 2d element area
 *  @param[in]  dt not used
 *  @param[in]  c a constant
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void inline fe_sw2_wd_bodyForce_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *rhs) {
    int i;
    
    double dt = vars[0];
    double c = vars[1];
    
    double rhs_x[NDONTRI] = {0., 0., 0.};
    double rhs_y[NDONTRI] = {0., 0., 0.};
    integrate_triangle_phi_h_df(djac, c, elem_head, *v, rhs_x, rhs_y);
    for (i=0; i<NDONTRI; i++) {
        rhs[i*3+1]+= rhs_x[i];
        rhs[i*3+2]+= rhs_y[i];
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Add the 2D shallow water wet-dry wrapped boundary pressure terms to the FE residual.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[inout] rhs the 2D elemental residual array contribution with wet/dry inclusion
 *  @param[in]  elem_nds the wet-dry nodal positions
 *  @param[in]  elem_head the 2d element nodal wet-dry depths
 *  @param[in]  grad_phi not used
 *  @param[in]  v not used
 *  @param[in]  v_wet_dry not used
 *  @param[in]  f not used
 *  @param[in]  f_wet_dry not used
 *  @param[in]  djac the 2d element area
 *  @param[in]  dt not used
 *  @param[in]  c a constant
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void inline fe_sw2_wd_boundaryPressure_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *elem_rhs) {
    double dt = vars[0];
    double c = vars[1];
    
    int nd1, nd2, iedge;
    double delx, dely;      /* the signed edge x and y lengths */
    double rhs[2] = {0., 0.};
    int nd1s[3] = {0,1,2};
    int nd2s[3] = {1,2,0};
    
    // loop over edges to evaulate edge integrals
    for (iedge=0; iedge<NDONTRI; iedge++) {
        // cjt :: this must agree with what's in the structs/sgrid.c file
        //if (iedge == 0)      {nd1=0; nd2=1;}
        //else if (iedge == 1) {nd1=1; nd2=2;}
        //else if (iedge == 2) {nd1=2; nd2=0;}
        //replacing if()
        nd1 = nd1s[iedge];
        nd2 = nd2s[iedge];
        
        delx = elem_nds[nd1].x - elem_nds[nd2].x;
        dely = elem_nds[nd2].y - elem_nds[nd1].y;
        
        rhs[0] = 0.; rhs[1] = 0.;
        integrate_line_phi_h_h(1., c, elem_head[nd1], elem_head[nd2], rhs); // set djac to 1 here so that we can multiply in this routine
        elem_rhs[nd1*3+1] += dely * rhs[0];
        elem_rhs[nd2*3+1] += dely * rhs[1];
        elem_rhs[nd1*3+2] += delx * rhs[0];
        elem_rhs[nd2*3+2] += delx * rhs[1];
    }
    
    //    // old adh way
    //    double one6 = 1. / 6.;
    //    double one12 = 1. / 12.;
    //    double delx, dely;        /* the signed edge x and y lengths */
    //    double sum, sum_h;
    //    double fact_start, fact_end;
    //
    //    double s = c;
    //    double h[3]; h[0] = elem_head[0]; h[1] = elem_head[1]; h[2] = elem_head[2];
    //
    //    int i;
    //    for (i=0; i<3; i++) {
    //        elem_rhs[i].x_eq = 0.;
    //        elem_rhs[i].y_eq = 0.;
    //    }
    //
    //    /* the first edge */
    //    delx = elem_nds[0].x - elem_nds[1].x;
    //    dely = elem_nds[1].y - elem_nds[0].y;
    //    sum_h = (h[0] + h[1]);
    //    sum = sum_h * sum_h * one12;
    //    fact_start = (sum + h[0] * h[0] * one6) * s;
    //    fact_end = (sum + h[1] * h[1] * one6) * s;
    //    elem_rhs[0].x_eq += dely * fact_start;
    //    elem_rhs[0].y_eq += delx * fact_start;
    //    elem_rhs[1].x_eq += dely * fact_end;
    //    elem_rhs[1].y_eq += delx * fact_end;
    //    //printf("1 :: delx: %30.20e dely: %30.20e fact_start: %30.20e fact_end: %30.20e :: h: %30.20e %30.20e %30.20e :: %30.20e %30.20e %30.20e %30.20e\n",delx,dely,fact_start,fact_end,h[0],h[1],h[2],elem_nds[0].x,elem_nds[1].x,elem_nds[2].x,elem_rhs[1].x_eq);
    //
    //    /* the second edge */
    //    delx = elem_nds[1].x - elem_nds[2].x;
    //    dely = elem_nds[2].y - elem_nds[1].y;
    //    sum_h = (h[1] + h[2]);
    //    sum = sum_h * sum_h * one12;
    //    fact_start = (sum + h[1] * h[1] * one6) * s;
    //    fact_end = (sum + h[2] * h[2] * one6) * s;
    //
    //    elem_rhs[1].x_eq += dely * fact_start;
    //    elem_rhs[1].y_eq += delx * fact_start;
    //    elem_rhs[2].x_eq += dely * fact_end;
    //    elem_rhs[2].y_eq += delx * fact_end;
    //    //printf("2 :: delx: %30.20e dely: %30.20e fact_start: %30.20e fact_end: %30.20e :: h: %30.20e %30.20e %30.20e :: %30.20e %30.20e %30.20e %30.20e\n",delx,dely,fact_start,fact_end,h[0],h[1],h[2],elem_nds[0].y,elem_nds[1].y,elem_nds[2].y,elem_rhs[1].x_eq);
    //
    //    /* the third edge */
    //    delx = elem_nds[2].x - elem_nds[0].x;
    //    dely = elem_nds[0].y - elem_nds[2].y;
    //    sum_h = (h[2] + h[0]);
    //    sum = sum_h * sum_h * one12;
    //    fact_start = (sum + h[2] * h[2] * one6) * s;
    //    fact_end = (sum + h[0] * h[0] * one6) * s;
    //
    //    elem_rhs[2].x_eq += dely * fact_start;
    //    elem_rhs[2].y_eq += delx * fact_start;
    //    elem_rhs[0].x_eq += dely * fact_end;
    //    elem_rhs[0].y_eq += delx * fact_end;
    //    //printf("3 :: delx: %30.20e dely: %30.20e fact_start: %30.20e fact_end: %30.20e :: h: %30.20e %30.20e %30.20e :: %30.20e %30.20e %30.20e %30.20e\n",delx,dely,fact_start,fact_end,h[0],h[1],h[2],elem_nds[0].z,elem_nds[1].z,s,elem_rhs[1].x_eq);   
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Add the 2D shallow water wet-dry wrapped density pressure terms to the FE residual.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[inout] rhs the 2D elemental residual array contribution with wet/dry inclusion
 *  @param[in]  elem_nds not used
 *  @param[in]  elem_head the 2d element nodal wet-dry depths
 *  @param[in]  grad_phi the elemental shape function gradients in configuration space
 *  @param[in]  v not used
 *  @param[in]  v_wet_dry not used
 *  @param[in]  f fully wet density
 *  @param[in]  f_wet_dry not used
 *  @param[in]  djac the 2d element area
 *  @param[in]  dt not used
 *  @param[in]  c a constant
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void inline fe_sw2_wd_densityPressure_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *rhs) {
    int i;
    double dt = vars[0];
    double c = vars[1];
    
    double rhs_x[NDONTRI] = {0., 0., 0.};
    double rhs_y[NDONTRI] = {0., 0., 0.};
    integrate_triangle_dphi_f_f_h(grad_phi, djac, c, elem_head, f, rhs_x, rhs_y);
    for (i=0; i<NDONTRI; i++) {
        rhs[i*3+1] += rhs_x[i];
        rhs[i*3+2] += rhs_y[i];
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Add the 2D shallow water wet-dry wrapped density body-force terms to the FE residual.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[inout] rhs the 2D elemental residual array contribution with wet/dry inclusion
 *  @param[in]  elem_nds not used
 *  @param[in]  elem_head the 2d element nodal wet-dry depths
 *  @param[in]  grad_phi not used
 *  @param[in]  v the fully wet bathymetry gradient
 *  @param[in]  v_wet_dry not used
 *  @param[in]  f fully we density
 *  @param[in]  f_wet_dry not used
 *  @param[in]  djac the 2d element area
 *  @param[in]  dt not used
 *  @param[in]  c a constant
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void inline fe_sw2_wd_densityBodyForce_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *rhs) {
    int i;
    double dt = vars[0];
    double c = vars[1];
    
    double rhs_x[NDONTRI] = {0., 0., 0.};
    double rhs_y[NDONTRI] = {0., 0., 0.};
    integrate_triangle_phi_h_g_df(djac, c, elem_head, f, *v, rhs_x, rhs_y);
    for (i=0; i<NDONTRI; i++) {
        rhs[i*3+1] += rhs_x[i];
        rhs[i*3+2] += rhs_y[i];
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Add the 2D shallow water wet-dry wrapped density boundary pressure terms to the FE residual.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[inout] rhs the 2D elemental residual array contribution with wet/dry inclusion
 *  @param[in]  elem_nds the wet-dry nodal positions
 *  @param[in]  elem_head the 2d element nodal wet-dry depths
 *  @param[in]  grad_phi not used
 *  @param[in]  v not used
 *  @param[in]  v_wet_dry not used
 *  @param[in]  f fully wet density
 *  @param[in]  f_wet_dry not used
 *  @param[in]  djac the 2d element area
 *  @param[in]  dt not used
 *  @param[in]  c a constant
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void inline fe_sw2_wd_densityBoundaryPressure_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *rhs) {
    
    double dt = vars[0];
    double c = vars[1];
    
    int nd1, nd2, iedge;
    double delx, dely;      /* the signed edge x and y lengths */
    int nd1s[3] = {0,1,2};
    int nd2s[3] = {1,2,0};
    
    // loop over edges to evaluate edge integrals
    for (iedge=0; iedge<NDONTRI; iedge++) {
        // cjt :: this must agree with what's in the structs/sgrid.c file
        //if (iedge == 0)      {nd1=0; nd2=1;}
        //else if (iedge == 1) {nd1=1; nd2=2;}
        //else if (iedge == 2) {nd1=2; nd2=0;}
        nd1 = nd1s[iedge];
        nd2 = nd2s[iedge];
        
        delx = elem_nds[nd1].x - elem_nds[nd2].x;
        dely = elem_nds[nd2].y - elem_nds[nd1].y;
        
        double rhst[2] = {0., 0.};
        integrate_line_phi_h_h_g(1., c, elem_head[nd1], elem_head[nd2], f[nd1], f[nd2], rhst); // set djac to 1 here so that we can multiply in this routine
        rhs[nd1*3+1] += dely * rhst[0];
        rhs[nd2*3+1] += dely * rhst[1];
        rhs[nd1*3+2] += delx * rhst[0];
        rhs[nd2*3+2] += delx * rhst[1];
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Wet/Dry routine for calculating the wet-dry integral of shallow water H, u and v:
 *  \f$  \int_{0}^{1} \int_{0}^{1-\widehat{x}} \big( c \, f(\widehat{x},\widehat{y})\big) d\widehat{y}\,d\widehat{x} \f$
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gary Brown
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \date      October, 2017
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[inout] rhs the 2D elemental residual array contribution with wet/dry inclusion
 *  @param[in]  elem_nds not used
 *  @param[in]  elem_head the 2d element nodal wet-dry depths
 *  @param[in]  grad_phi not used
 *  @param[in]  v not used
 *  @param[in]  v_wet_dry wet-dry velocities (integrated here)
 *  @param[in]  f not used
 *  @param[in]  f_wet_dry the scalar function (integrated here)
 *  @param[in]  djac the 2d element area
 *  @param[in]  dt not used
 *  @param[in]  c a constant
 *
 * \note CJT \:: must have same call arguements as all generic wet-dry wrapped functions
 * \note CJT \:: djac is used here because later, in the wet-dry wrapper, it is divided by the wet area (djac)
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void inline fe_sw2_wd_integrate_triangle_f(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *rhs) {
    
    double dt = vars[0];
    double c = vars[1];
    
    // cjt :: doing this crap to match the trunk while debugging
    //double t1 = (f_wet_dry[0]   + f_wet_dry[1]   + f_wet_dry[2]) / 3.0; // gives slightly different results
    double t1 = (elem_head[0]   + elem_head[1]   + elem_head[2]) / 3.0;
    double t2 = 0.0;
    double t3 = 0.0;
    
    if (v_wet_dry != NULL) {
        t2 = (v_wet_dry[0].x + v_wet_dry[1].x + v_wet_dry[2].x) / 3.0;
        t3 = (v_wet_dry[0].y + v_wet_dry[1].y + v_wet_dry[2].y) / 3.0;
    }

    //why all node 0???
    rhs[0] += djac * t1;
    rhs[1] += djac * t2;
    rhs[2] += djac * t3;
    
    return;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Add the 2D shallow water wet-dry wrapped *nonconservative convection terms to the FE residual.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[inout] rhs the 2D elemental residual array contribution with wet/dry inclusion
 *  @param[in]  elem_nds not used
 *  @param[in]  elem_head the 2d element nodal wet-dry depths
 *  @param[in]  grad_phi the elemental shape function gradients in configuration space
 *  @param[in]  v not used
 *  @param[in]  v_wet_dry wet-dry velocities
 *  @param[in]  f not used
 *  @param[in]  f_wet_dry wet-dry velocities
 *  @param[in]  djac the 2d element area
 *  @param[in]  dt time-step
 *  @param[in]  c a constant
 *
 *  \details Solves the following weak, non-conservative wet/dry discrete equation: \n
 *  \f{eqnarray*}{
 *   \residDA{0}{}{c}   &=& \intg{\overline{u} \deriv{H}{x} + \overline{H} \deriv{u}{x} + \overline{v} \deriv{H}{y} + \overline{H} \deriv{v}{y} }{2d}{e}  \\
 *   \residDA{0}{}{mx}  &=& \intg{\overline{u} \deriv{u}{x} + \overline{v} \deriv{u}{x}}{2d}{e} \\
 *   \residDA{0}{}{my}  &=& \intg{\overline{u} \deriv{v}{x} + \overline{v} \deriv{v}{x}}{2d}{e}
 *  \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void inline fe_sw2_wd_gls_convection_triangle(SVECT *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, double *rhs) {
    
    double dt = vars[0];
    double c = vars[1];
    
    // calculate wet-dry elemental averages
    double h_avg = one_3 * (f_wet_dry[0]   + f_wet_dry[1]   + f_wet_dry[2]);
    double u_avg = one_3 * (v_wet_dry[0].x + v_wet_dry[1].x + v_wet_dry[2].x);
    double v_avg = one_3 * (v_wet_dry[0].y + v_wet_dry[1].y + v_wet_dry[2].y);
    
    // calculate wet-dry gradients (constant on triangles)
    double dhdx = grad_phi[0].x * f_wet_dry[0]   + grad_phi[1].x * f_wet_dry[1]   + grad_phi[2].x * f_wet_dry[2];
    double dhdy = grad_phi[0].y * f_wet_dry[0]   + grad_phi[1].y * f_wet_dry[1]   + grad_phi[2].y * f_wet_dry[2];
    double dudx = grad_phi[0].x * v_wet_dry[0].x + grad_phi[1].x * v_wet_dry[1].x + grad_phi[2].x * v_wet_dry[2].x;
    double dudy = grad_phi[0].y * v_wet_dry[0].x + grad_phi[1].y * v_wet_dry[1].x + grad_phi[2].y * v_wet_dry[2].x;
    double dvdx = grad_phi[0].x * v_wet_dry[0].y + grad_phi[1].x * v_wet_dry[1].y + grad_phi[2].x * v_wet_dry[2].y;
    double dvdy = grad_phi[0].y * v_wet_dry[0].y + grad_phi[1].y * v_wet_dry[1].y + grad_phi[2].y * v_wet_dry[2].y;
    
    // the convection terms chain-ruled out
    rhs[0] += djac * dt * c * (u_avg * dhdx + h_avg * dudx + v_avg * dhdy + h_avg * dvdy);
    rhs[1] += djac * dt * c * (u_avg * dudx + v_avg * dudy);
    rhs[2] += djac * dt * c * (u_avg * dvdx + v_avg * dvdy);
}
