/* This routines calculates and returns the radial near-bed velocity due to vorticity */

#include "global_header.h"
SVECT2D sed_vorticity_velocity_components (double depth, SVECT2D velocity, SCON con, int inode) {
  
    double u_r;			/* near bed radial velocity magnitude */
    double c_2;			/* coefficient */
    double vel_mag;		/* depth averaged velocity magnitude */
    SVECT2D ave_vec;		/* the streamwise vector */
    SVECT2D e_r;			/* radial vector (normal to streamwise) */

    // return variable
    SVECT2D vor_vel;  svect2d_init(&(vor_vel));
    
    vel_mag = VECT2D_MAG(velocity);
    if (vel_mag > SMALL) {
        ave_vec.x = velocity.x / vel_mag;
        ave_vec.y = velocity.y / vel_mag;
        e_r.x = -ave_vec.y;
        e_r.y =  ave_vec.x;
        
        c_2 = sqrt (6. * con.property[1]);
        
        u_r = 6. * con.concentration[inode] * depth / c_2;
        if (u_r >  0.25 * vel_mag) u_r =  0.25 * vel_mag;
        if (u_r < -0.25 * vel_mag) u_r = -0.25 * vel_mag;
	  
        vor_vel.x = u_r * e_r.x;
        vor_vel.y = u_r * e_r.y;
    }

    return vor_vel;
}
