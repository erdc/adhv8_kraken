/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates the continuity advection velocities for the diffusive wave model.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \date      October, 2017
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[out] vel an array of 2D vector of velocities
 * @param[in]  nnodes the total number of nodes on the element
 * @param[in]  depth total water depth
 * @param[in]  avg_depth the elementally averaged depth
 * @param[in]  z the z-coordinate of the elemental nodes
 * @param[in]  grad_phi a pointer to the array of 2D test function gradients
 * @param[in]  roughness roughness friction factor
 * @param[in]  mannings_units_constant Mannings unit constant
 *
 * \note   CJT \:: all the velocities on the element should be the same.  I store the value in \n
 *                 all nodes in case this changes and also to try something in the wet/dry wrapper.
 * \remark
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

SVECT2D getDiffusiveWaveVelocities(int nnodes, double *depth, double avg_depth, double *z, SVECT2D *grad_phi, double roughness, double manning_units_constant, SVECT2D *vel) {
    
    int i;
    double dz_dx = 0., dz_dy = 0., elev = 0.;
    double depthpow2by3[nnodes];
    for (i=0; i<nnodes; i++) {
        elev = depth[i] + z[i];
        dz_dx += elev * grad_phi[i].x;
        dz_dy += elev * grad_phi[i].y;
        depthpow2by3[i]=pow(depth[i]*depth[i], 1. / 3.); /* depth^(2/3) */
    }
    
    double slope_check = sqrt(dz_dx * dz_dx + dz_dy * dz_dy);
    if(slope_check < SMALL) {
        dz_dx = 0.;
        dz_dy = 0.;
    } else {
        double sign_dz_dx = 1.;
        if (is_double_small(dz_dx) == NO) sign_dz_dx = fabs(dz_dx) / dz_dx;
        double sign_dz_dy = 1.;
        if (is_double_small(dz_dy) == NO) sign_dz_dy = fabs(dz_dy) / dz_dy;
        
        if (is_double_small(roughness) == NO) {
            dz_dx = -sign_dz_dx * manning_units_constant * sqrt(fabs(dz_dx)) / roughness;
            dz_dy = -sign_dz_dy * manning_units_constant * sqrt(fabs(dz_dy)) / roughness;
        } else {
            dz_dx = 0.;
            dz_dy = 0.;
        }
    }
    
    svect2d_init_array(vel,nnodes);
    SVECT2D v; v.x = 0.; v.y = 0.;
    
    if (avg_depth > 0) {
        avg_depth = pow(avg_depth, 2. / 3.);
    } else {
        avg_depth = 0.;
    }
    
    v.x = dz_dx * avg_depth;
    v.y = dz_dy * avg_depth;
    // Gajanan gkc: Gaurav mentioned in an email that velocity must have depth^(2/3)
    for (i=0; i<nnodes; i++) {
        vel[i].x = dz_dx*depthpow2by3[i]; //    vel[i].x = dz_dx*depth[i]; //    
        vel[i].y = dz_dy*depthpow2by3[i]; //    vel[i].y = dz_dy*depth[i]; //    
        if (depth[i] < 0) {
            vel[i].x = 0.;
            vel[i].y = 0.;
        }
    }
    
    return v;
}
