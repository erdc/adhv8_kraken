/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes the eddy viscosity coefficients for the SW system.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[out] ev_st (double *) the streamline eddy viscosity coefficient
 * @param[out] ev_tr (double *) the transverse eddy viscosity coefficient
 * @param[in]  eev_mode (int) a flag to determine which eddy viscosity calculation to use
 * @param[in]  eev_coef (double) a coefficient
 * @param[in]  g (double) gravity
 * @param[in]  drying_lower_limit (double) a minimum depth
 * @param[in]  h_fric (double) the depth used in friction calculations
 * @param[in]  roughness (double) friction roughness
 * @param[in]  area (double) the element area
 * @param[in]  grad_u (SVECT) the x-velocity gradient
 * @param[in]  grad_v (SVECT) the y-velocity gradient
 * @param[in]  avg_vmag (double) the elementally averaged velocity
 * @param[in]  wd_flag (int) a flag that determines if the element is wet/dry
 * 
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_sw2_get_EEVF(int eev_mode, double eev_coef, double g, double drying_lower_limit, double h_fric, double roughness, double area, SVECT2D grad_u, SVECT2D grad_v, double avg_vmag, int wd_flag, double *ev_st, double *ev_tr) {
    
    // estimated eddy viscosity
    double ev_min, vmin;
    vmin = 0.1 * sqrt(g * drying_lower_limit); /* GLB 11-13 use peclet number of 40 with min velocity to define min mixing */
    if (eev_mode == 1) {
        *ev_st = 0.0;
        *ev_tr = eev_coef * 2.0 * h_fric * sqrt(roughness) * avg_vmag;
        ev_min = (eev_coef / 0.5) * vmin * sqrt(area) / 40.;
        *ev_tr = MAX(*ev_tr, ev_min);
    }
    else if (eev_mode == 2) {
        *ev_st = 1.3 * h_fric * sqrt(roughness) * avg_vmag;
        *ev_tr = eev_coef * 0.92 * h_fric * sqrt(roughness) * avg_vmag;    /* Scaling by 10 therefore 0.092 is 0.92 */
        ev_min = (eev_coef / 0.5) * vmin * sqrt(area) / 40.;
        *ev_tr = MAX(*ev_tr, ev_min);
    }
    else if (eev_mode == 3) { /* Smag */
        *ev_st = 0.;
        *ev_tr = eev_coef * eev_coef * area * sqrt(grad_u.x * grad_u.x + grad_v.y * grad_v.y + 0.5 * (grad_u.y + grad_v.x) * (grad_u.y + grad_v.x) );
        ev_min = (eev_coef / 0.2) * vmin * sqrt(area) / 40.;
        *ev_tr = MAX(*ev_tr, ev_min);
    }
    if (wd_flag == 1 || wd_flag == 2) {
        *ev_st = 0.0;
        *ev_tr = MAX(*ev_tr, 5. * ev_min);
    }
}

