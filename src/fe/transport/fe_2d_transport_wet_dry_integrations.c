/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_2d_transport_wet_dry_integrations.c This file collects the 2D shallow water transport wet/dry integrations.
 *   All function arguements must be consistent with fe_sw2_wet_dry_wrapper.
*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Wet/Dry routine for calculating the wet-dry convection terms to the FE residual.
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
 *  @param[in]  v fully no used
 *  @param[in]  v_wet_dry the wet-dry velocities
 *  @param[in]  f not used
 *  @param[in]  f_wet_dry the wet-dry constituent concentrations (elem_c)
 *  @param[in]  djac the 2d element area
 *  @param[in]  dt not used
 *  @param[in]  c a constant
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_2d_transport_wd_convection_triangle(SVECT2D *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, DOF_3 *elem_rhs) {
    
    double dt = vars[0];
    double c = vars[1];
    
    int i;
    double rhs[NDONTRI]; sarray_init_dbl(rhs,NDONTRI);
//    printf("c: %20.10e dt: %20.10e djac: %20.10e all: %20.10e\n",c,dt,djac,c*dt*djac);
//    printf("elem_head :: %20.10e %20.10e %20.10e\n",elem_head[0],elem_head[1],elem_head[2]);
//    printf("v_wet_dry.x :: %20.10e %20.10e %20.10e\n",v_wet_dry[0].x,v_wet_dry[1].x,v_wet_dry[2].x);
//    printf("v_wet_dry.y :: %20.10e %20.10e %20.10e\n",v_wet_dry[0].y,v_wet_dry[1].y,v_wet_dry[2].y);
//    printf("grad_phi.x :: %20.10e %20.10e %20.10e\n",grad_phi[0].x,grad_phi[1].x,grad_phi[2].x);
//    printf("grad_phi.y :: %20.10e %20.10e %20.10e\n",grad_phi[0].y,grad_phi[1].y,grad_phi[2].y);
//    printf("f_wet_dry :: %20.10e %20.10e %20.10e\n",f_wet_dry[0],f_wet_dry[1],f_wet_dry[2]);
    
    integrate_triangle_gradPhi_dot_f_g_v(grad_phi, djac, dt * c, elem_head, f_wet_dry, v_wet_dry, rhs);
    for (i=0; i<NDONTRI; i++) {elem_rhs[i].c_eq += rhs[i];}
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Wet/Dry routine for calculating the wet-dry wrapped temporal terms to the FE residual.
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
 *  @param[in]  v fully not used
 *  @param[in]  v_wet_dry not used
 *  @param[in]  f not used
 *  @param[in]  f_wet_dry the wet-dry constituent concentrations (elem_c)
 *  @param[in]  djac the 2d element area
 *  @param[in]  dt not used
 *  @param[in]  c a constant
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_2d_transport_wd_temporal_triangle(SVECT2D *elem_nds, double *elem_head, SVECT2D *grad_phi, SVECT2D *v, SVECT2D *v_wet_dry, double *f, double *f_wet_dry, double djac, double *vars, DOF_3 *elem_rhs) {
    
    double dt = vars[0];
    double c = vars[1];
    
#ifdef _DEBUG
    assert(dt == 1.);
#endif
    
    int i;
    double rhs[NDONTRI]; sarray_init_dbl(rhs,NDONTRI);
    integrate_triangle_phi_f_g(djac, dt * c, elem_head, f_wet_dry, rhs);
    for (i=0; i<NDONTRI; i++) {elem_rhs[i].c_eq += rhs[i];}
}
