/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_sw2_boundary_integrals.c This file collections functions responsible for
 *          the 2D shallow water body contributions to the elemental residual.              */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG_LOCAL = OFF;
static int DEBUG_PICKETS = OFF;
static int DEBUG_NODE_ID = UNSET_INT;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 1D element explicit flow contributions to 2D-SW equations
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] elem_rhs    the 1D elemental residual array
 * @param[in]  ie            the 1D element ID
 * @param[in]  nrml          the 1D element normal
 * @param[in]  djac          the 1D element jacobian
 * @param[in]  dt            the current model time-step
 * @param[in]  flux          a user-supplied flux
 * @param[in]  h_avg         the elementally averaged depth
 * @param[in]  u             x-velocity
 * @param[in]  v             y-velocity
 *
 * \f{eqnarray*}{  \residDA{i}{}{c}   &=& dt * \intbe{\phidd{i} \, q_n}{e}    \f}
 *
 * For the momentum contributions: \n
 * -- if the flux is into the domain, the normal flux is split into x,y velocities, so that
 *  \f{eqnarray*}{
 *     u^* = (q_n / \overline{\depth{h}}) \, n_x \\
 *     v^* = (q_n / \overline{\depth{h}}) \, n_y
 *  \f}
 * and added to the residual as
 *  \f{eqnarray*}{
 *   \residDA{i}{}{mx}  &=& dt * \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \ub{h} \, \depth{h})} = dt * \intbe{\phidd{i} \, q_n \, u^*}{e} \\
 *   \residDA{i}{}{my}  &=& dt * \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \vb{h} \, \depth{h})} = dt * \intbe{\phidd{i} \, q_n \, v^*}{e}
 *  \f}
 * -- if the flux is out of the domain, an implicity u,v is used so that
 *  \f{eqnarray*}{
 *   \residDA{i}{}{mx}  &=& dt * \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \ub{h} \, \depth{h})} = dt * \intbe{\phidd{i} \, q_n \, u}{e} \\
 *   \residDA{i}{}{my}  &=& dt 8 \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \vb{h} \, \depth{h})} = dt * \intbe{\phidd{i} \, q_n \, v}{e}
 *  \f}
 *
 *  \details  Forms the 1d element matrix for the flux terms in the shallow water equations. \n
 *            This is the unit discharge out of the boundary. Positive is out, Negative is in.
 * NOTE: CJT \:: all units here must be equivalent to Int [d(uh)/dt * dxdy] = m^4/s^2
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_sw2_1d_explicit_flow(double *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, double flux, double h_avg, double *u, double *v, int DEBUG, int DEBUG_LOCAL) {
    int i;
    double rhs_c_eq[NDONSEG], rhs_x_eq[NDONSEG], rhs_y_eq[NDONSEG]; // local rhs
    sarray_init_dbl(rhs_c_eq, NDONSEG);
    sarray_init_dbl(rhs_x_eq, NDONSEG);
    sarray_init_dbl(rhs_y_eq, NDONSEG);
    integrate_line_phi(djac, 1., rhs_c_eq);
    if (flux < 0) {    // flow is into the model
        double u_star = 0., v_star = 0.;
        if (h_avg > 0) {
            u_star = flux * nrml.x / h_avg; // average directional velocity
            v_star = flux * nrml.y / h_avg; // average directional velocity
        }
        integrate_line_phi(djac, u_star, rhs_x_eq);
        integrate_line_phi(djac, v_star, rhs_y_eq);
        
    } else {           // flow is out of the model
        integrate_line_phi_f(djac, 1., u, rhs_x_eq);
        integrate_line_phi_f(djac, 1., v, rhs_y_eq);
    }
    
    double t1 = dt * flux;
    for (i=0; i<NDONSEG; i++) {
        elem_rhs[i*3] += t1 * rhs_c_eq[i];
        elem_rhs[i*3+1] += t1 * rhs_x_eq[i];
        elem_rhs[i*3+2] += t1 * rhs_y_eq[i];
    }
    
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        double rhs[NDONSEG*3];
        for (i=0; i<NDONSEG; i++) {
            rhs[i*3] = t1 * rhs_c_eq[i];
            rhs[i*3+1] = t1 * rhs_x_eq[i];
            rhs[i*3+2] = t1 * rhs_y_eq[i];
        }
        printScreen_rhs_3dof("1D SW EXPLICIT FLOW: ", NDONSEG, ie, elem1d.nodes, rhs);
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 1D element implicit flow contributions to 2D-SW equations
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] elem_rhs    the 1D elemental residual array
 * @param[in]  ie            the 1D element ID
 * @param[in]  nrml          the 1D element normal
 * @param[in]  djac          the 1D element jacobian
 * @param[in]  dt            the current model time-step
 * @param[in]  v             the hydrodynamic velocity
 * @param[in]  h             the water depth
 *
 * \f{eqnarray*}{  \residDA{i}{}{c}   &=& dt * \intbe{\phidd{i} \, q_n}{e} \\
 *                 \residDA{i}{}{mx}  &=& dt * \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \ub{h} \, \depth{h})} \\
 *                 \residDA{i}{}{my}  &=& dt * \bcConv{}{e}{\phidd{i}}{(\velb{h} \, \vb{h} \, \depth{h})}
 * \f}
 *
 * \note CJT\::  Gaurav *thinks the "r" variable is a prox for v so he can eliminate eddies forming on the boundary
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_sw2_1d_implicit_flow(double *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, SVECT2D *v, double *h, int DEBUG, int DEBUG_LOCAL) {
    int i;
    double rhs_c_eq[NDONSEG], rhs_x_eq[NDONSEG], rhs_y_eq[NDONSEG]; // local rhs
    sarray_init_dbl(rhs_c_eq, NDONSEG);
    sarray_init_dbl(rhs_x_eq, NDONSEG);
    sarray_init_dbl(rhs_y_eq, NDONSEG);
    
    double u1 =  v[0].x, u2 = v[1].x;
    double v1 =  v[0].y, v2 = v[1].y;
    double h1 =  h[0],   h2 = h[1];
    double ur1 = v[0].x, ur2 = v[1].x;
    double vr1 = v[0].y, vr2 = v[1].y;
    
    
    double dir_vect[NDONSEG];
    dir_vect[0] = v[0].x * nrml.x + v[0].y * nrml.y;
    dir_vect[1] = v[1].x * nrml.x + v[1].y * nrml.y;
    if (dir_vect[0] < SMALL) {
        ur1 = 0.;
        vr1 = 0.;
    }
    if (dir_vect[1] < SMALL) {
        ur2 = 0.;
        vr2 = 0.;
    }
    
    integrate_line_phi_h_v_dot_n(djac, 1., h, v, nrml, rhs_c_eq); // phi_i * u * h * nx + phi_i * v * h * ny
    integrate_line_phi_f_g_h(djac, nrml.x, ur1, ur2, u1, u2, h1, h2, rhs_x_eq); // phi_i * ur * u * h * nx
    integrate_line_phi_f_g_h(djac, nrml.y, vr1, vr2, v1, v2, h1, h2, rhs_y_eq); // phi_i * vr * v * h * ny
    integrate_line_phi_f_g_h(djac, nrml.y, ur1, ur2, v1, v2, h1, h2, rhs_x_eq); // phi_i * ur * v * h * ny
    integrate_line_phi_f_g_h(djac, nrml.x, vr1, vr2, u1, u2, h1, h2, rhs_y_eq); // phi_i * vr * u * h * nx
    
    for (i=0; i<NDONSEG; i++) {
        elem_rhs[i*3] += dt * rhs_c_eq[i];
        elem_rhs[i*3+1] += dt * rhs_x_eq[i];
        elem_rhs[i*3+2] += dt * rhs_y_eq[i];
    }
    
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        double rhs[NDONSEG*3];
        for (i=0; i<NDONSEG; i++) {
            rhs[i*3] = dt * rhs_c_eq[i];
            rhs[i*3+1] = dt * rhs_x_eq[i];
            rhs[i*3+2] = dt * rhs_y_eq[i];
        }
        printScreen_rhs_3dof("1D SW IMPLICIT FLOW: ", NDONSEG, ie, elem1d.nodes, rhs);
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 1D element wall friction contributions to 2D-SW equations
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] elem_rhs    the 1D elemental residual array
 * @param[in]  ie            the 1D element ID
 * @param[in]  nrml          the 1D element normal
 * @param[in]  djac          the 1D element jacobian
 * @param[in]  dt            the current model time-step
 * @param[in]  u             x-velocity
 * @param[in]  v             y-velocity
 * @param[in]  resistance    resistance
 *
 *  \details Integrates the weak, discrete outflow boundary terms: \n
 *  \f{eqnarray*}{
 *   \residDA{i}{}{mx}  &=& dt * \bcFriction{}{e}{\phidd{i}}{\ub{h}}{\velb{h}} \\
 *   \residDA{i}{}{my}  &=& dt * \bcFriction{}{e}{\phidd{i}}{\vb{h}}{\velb{h}}
 *  \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_sw2_1d_wall_friction(double *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, double *u, double *v, double *resistance, int DEBUG, int DEBUG_LOCAL) {
    
    int i;
    double rhs_x_eq[NDONSEG], rhs_y_eq[NDONSEG]; // local rhs
    sarray_init_dbl(rhs_x_eq, NDONSEG);
    sarray_init_dbl(rhs_y_eq, NDONSEG);
    
    integrate_line_phi_f(djac, fabs(nrml.y), u, rhs_x_eq);
    integrate_line_phi_f(djac, fabs(nrml.x), v, rhs_y_eq);
    
    for (i=0; i<NDONSEG; i++) {
        elem_rhs[i*3+1] += dt * resistance[i] * rhs_x_eq[i];
        elem_rhs[i*3+2] += dt * resistance[i] * rhs_y_eq[i];
    }
    
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        double rhs[NDONSEG*3];
        for (i=0; i<NDONSEG; i++) {
            rhs[i*3] = 0;
            rhs[i*3+1] = dt * resistance[i] * rhs_x_eq[i];
            rhs[i*3+2] = dt * resistance[i] * rhs_y_eq[i];
        }
        printScreen_rhs_3dof("1D SW WALL FRICTION ", NDONSEG, ie, elem1d.nodes, rhs);
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Returns the 1D element pressure contributions to 2D-SW equations
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[inout] elem_rhs    the 1D elemental residual array
 * @param[in]  ie            the 1D element ID
 * @param[in]  nrml          the 1D element normal
 * @param[in]  djac          the 1D element jacobian
 * @param[in]  dt            the current model time-step
 * @param[in]  h             the water depth
 * @param[in]  g             gravity
 *
 * \f{eqnarray*}{
 *                 \residDA{i}{}{mx}  &=& dt * \bcPressureDD{}{e}{\phidd{i}}{(\depth{h})}{x} \\
 *                 \residDA{i}{}{my}  &=& dt * \bcPressureDD{}{e}{\phidd{i}}{(\depth{h})}{y}
 * \f}
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_sw2_1d_pressure(double *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, double *h, double g, int DEBUG, int DEBUG_LOCAL) {
    
    int i;
    double rhs[NDONSEG] = {0., 0.};
    
    
    integrate_line_phi_h_h(djac, 1., h[0], h[1], rhs);
    
    double t1 = dt * one_2 * g;
    for (i=0; i<NDONSEG; i++) {
        elem_rhs[i*3+1] += t1 * nrml.x * rhs[i];
        elem_rhs[i*3+2] += t1 * nrml.y * rhs[i];
    }
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        double rhs2[NDONSEG*3];
        for (i=0; i<NDONSEG; i++) {
            rhs2[i*3] = 0;
            rhs2[i*3+1] = t1 * nrml.x * rhs[i];
            rhs2[i*3+2] = t1 * nrml.y * rhs[i];
        }
        printScreen_rhs_3dof("1D SW PRESSURE: ", NDONSEG, ie, elem1d.nodes, rhs2);
    }
#endif
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Returns the 1D element density pressure contributions to 2D-SW equations
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Gary Brown, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[inout] elem_rhs    the 1D elemental residual array
 * @param[in]  ie            the 1D element ID
 * @param[in]  nrml          the 1D element normal
 * @param[in]  djac          the 1D element jacobian
 * @param[in]  dt            the current model time-step
 * @param[in]  h             the water depth
 * @param[in]  d             water density
 * @param[in]  g             gravity
 *
 * \f{eqnarray*}{
 *                 \residDA{i}{}{mx}  &=& dt * \bcPressureDD{}{e}{\phidd{i}}{(\rho^h \, \depth{h})}{x} \\
 *                 \residDA{i}{}{my}  &=& dt * \bcPressureDD{}{e}{\phidd{i}}{(\rho^h \, \depth{h})}{y}
 * \f}
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_sw2_1d_density_pressure(double *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, double *h, double *d, double g, int DEBUG, int DEBUG_LOCAL) {
    
    int i;
    double rhs[NDONSEG] = {0., 0.};
    
    integrate_line_phi_h_h_g(djac, 1., h[0], h[1], d[0], d[1], rhs);
    
    double t1 = dt * one_2 * g;
    for (i=0; i<NDONSEG; i++) {
        elem_rhs[i*3+1] += t1 * nrml.x * rhs[i];
        elem_rhs[i*3+2] += t1 * nrml.y * rhs[i];
    }
    
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        double rhs2[NDONSEG*3];
        for (i=0; i<NDONSEG; i++) {
            rhs2[i*3] = 0;
            rhs2[i*3+1] = t1 * nrml.x * rhs[i];
            rhs2[i*3+2] = t1 * nrml.y * rhs[i];
        }
        printScreen_rhs_3dof("1D SW DENITY PRESSURE: ", NDONSEG, ie, elem1d.nodes, rhs2);
    }
#endif
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Returns the 1D element explicit flow and pressure contributions to 2D-SW equations
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] elem_rhs    the 1D elemental residual array
 * @param[in]  ie            the 1D element ID
 * @param[in]  nrml          the 1D element normal
 * @param[in]  djac          the 1D element jacobian
 * @param[in]  dt            the current model time-step
 * @param[in]  v             the water velocity
 * @param[in]  user_velocity the water velocity prescribed by user
 * @param[in]  h             the water depth
 *
 * \f{eqnarray*}{  \residDA{i}{}{c}   &=& dt * \intbe{\phidd{i} \, (\velb{h} \cdot \nrml) \, \depth{h}}{e} \\
 *                 \residDA{i}{}{mx}  &=& dt * \intbe{\phidd{i} \, q_n \, \ub{}}{e}\\
 *                 \residDA{i}{}{my}  &=& dt * \intbe{\phidd{i} \, q_n \, \vb{}}{e}
 * \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_sw2_1d_hvel(double *elem_rhs, SELEM_1D elem1d, int ie, SVECT2D nrml, double djac, double dt, SVECT2D *v, SVECT2D user_velocity, double *h, int DEBUG, int DEBUG_LOCAL){
    
    int i;
    double rhs[NDONSEG] = {0., 0.};
    integrate_line_phi_h_v_dot_n(djac, 1., h, v, nrml, rhs); // phi_i * u * h * nx + phi_i * v * h * ny
    
    for (i=0; i<NDONSEG; i++) {
        elem_rhs[i*3] += dt * rhs[i];
        elem_rhs[i*3+1] += dt * rhs[i] * user_velocity.x;
        elem_rhs[i*3+2] += dt * rhs[i] * user_velocity.y;
    }
    
#ifdef _DEBUG
    if (DEBUG == ON || DEBUG_LOCAL == ON) {
        double rhs2[NDONSEG*3];
        for (i=0; i<NDONSEG; i++) {
            rhs2[i*3] = dt * rhs[i];
            rhs2[i*3+1] = dt * rhs[i] * user_velocity.x;
            rhs2[i*3+2] = dt * rhs[i] * user_velocity.y;
        }
        printScreen_rhs_3dof("1D HVEL ", NDONSEG, ie, elem1d.nodes, rhs2);
    }
#endif
}

