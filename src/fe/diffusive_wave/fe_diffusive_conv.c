/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates the 2D TRIANGULAR element matrix contributions for the Galerkin *AND*
 *             SUPG convective terms for the diffusive wave equations.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \date      October, 2017
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[inout] elem_rhs the 2D elemental residual array contribution with wet/dry inclusion
 *  @param[in]  elem_nds
 *  @param[in]  depth the 2d element nodal wet-dry depths
 *  @param[in]  grad_phi
 *  @param[in]  v non-wet-dry velocities
 *  @param[in]  wd_vel  wet-dry velocities
 *  @param[in]  f non wet-dry function
 *  @param[in]  wd_head wet-dry head
 *  @param[in]  djac the 2d element area
 *  @param[in]  vars  an array of doubles to be passed such at dt, etc.
 *
 *
 * \note   CJT \::  what
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#include "global_header.h"
static int DEBUG = OFF;

void fe_diffusive_conv(SVECT *elem_nds, double *depth, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *wd_vel, double *f, double *wd_head, double djac, double *vars, DOF_3 *elem_rhs) {
    
    int i;
    double rhs[NDONTRI]; sarray_init_dbl(rhs,NDONTRI);
    double dt = vars[0];
    double tau_pg = vars[1];
    double mannings = vars[2];
    double roughness = vars[3];
    double z[3];
    SVECT2D wdvel[NDONTRI];
    z[0] = elem_nds[0].z;
    z[1] = elem_nds[1].z;
    z[2] = elem_nds[2].z;
    
    
    // get wet-dry elementally averaged depth
    double elem_avg_depth = sarray_avg_dbl(wd_head, NDONTRI);
    
    // get wet-dry nodal and elementally averaged velocities (based on wet-dry depth)
    SVECT2D elem_avg_vel = getDiffusiveWaveVelocities(NDONTRI, wd_head, elem_avg_depth, z, grad_phi, roughness, mannings, wdvel);
    
    // Galerkin contribution with elementally averaged depth and nodal velocities :: note cjt :: results not as good as elemental velocities
    //integrate_triangle_gradPhi_dot_v(grad_phi, djac, elem_avg_depth, wd_vel, rhs);
    // Galerkin contribution with elementally averaged depth and elementally averaged velocities
    integrate_triangle_gradPhi_dot_vbar(grad_phi, djac, elem_avg_depth, elem_avg_vel, rhs);
    for (i=0; i<NDONTRI; i++) {
        elem_rhs[i].c_eq -= dt * rhs[i];
    }
    
    // Petrov contribution -- averaged depth used here
    double g = 9.8;
    int tau_method_flag = 0; // not used
    int le_method_flag = 2; // area based
    double tau =  fe_get_supg_tau_sw(NDONTRI, elem_nds, g, elem_avg_depth, elem_avg_vel.x, elem_avg_vel.y,0., NULL, NULL, NULL, djac, tau_pg, 2, tau_method_flag, le_method_flag);
    
    // -- strong advection term is constant for triangles
    double convection_bar = elem_avg_depth * (grad_phi[0].x * wd_vel[0].x + grad_phi[1].x * wd_vel[1].x + grad_phi[2].x * wd_vel[2].x +
                                              grad_phi[0].y * wd_vel[0].y + grad_phi[1].y * wd_vel[1].y + grad_phi[2].y * wd_vel[2].y);
    
    sarray_init_dbl(rhs,NDONTRI);
    integrate_triangle_gradPhi_dot_v(grad_phi, djac, convection_bar, wd_vel, rhs);
    for (i=0; i<NDONTRI; i++) {
        elem_rhs[i].c_eq += tau * dt * rhs[i];
    }
    
#ifdef _DEBUG
    if (DEBUG) {
        Is_Double_Inf_or_NaN(elem_avg_depth ,__FILE__ ,__LINE__);
        Is_DoubleArray_Inf_or_NaN(z, 3 ,__FILE__ ,__LINE__);
        Is_Double_Inf_or_NaN(elem_avg_vel.x,__FILE__ ,__LINE__);
        Is_Double_Inf_or_NaN(elem_avg_vel.y,__FILE__ ,__LINE__);
        Is_DoubleArray_Inf_or_NaN(wd_head, 3 ,__FILE__ ,__LINE__);
        Is_Double_Inf_or_NaN(djac ,__FILE__ ,__LINE__);
        Is_Double_Inf_or_NaN(tau ,__FILE__ ,__LINE__);
        
        Is_Double_Inf_or_NaN(elem_rhs[0].c_eq ,__FILE__ ,__LINE__);
        Is_Double_Inf_or_NaN(elem_rhs[1].c_eq ,__FILE__ ,__LINE__);
        Is_Double_Inf_or_NaN(elem_rhs[2].c_eq ,__FILE__ ,__LINE__);
    }
#endif
}
