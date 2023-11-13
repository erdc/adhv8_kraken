/* ADH Version 2.0.0 7-11 */
/* calculates the depth average velocity and depth */
#include "global_header.h"

void tl_calculate_depthavgvel(SGRID *grid, SSW_3D *sw3d) {
    
    ID_LIST_ITEM *ptr;
    int nd=0;
    int inode=0;
    double vel_x_top=0., vel_y_top=0., ele_top=0., ele_bot=0., vel_x_bot=0., vel_y_bot=0.;
    double sum_x=0., sum_y=0., surface=0., depth=0.;
    
    if (grid->ndim != 3) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> input grid must be 3d\n");
    }
    
    /* Loop over the surface nodes */
    for (inode=0; inode<grid->nnodes_sur; inode++) {
        ptr = grid->vertical_list[inode];
        nd = ptr->id;
        /* Gajanan gkc - commeting out if condition. Bug fix, since ghost nodes weren't getting updated! */
        /* Also important for parallel XDMF output for properly plotting solution at ghost nodes.        */
        /* It is important that fe_sw3_hvel_update be called BEFORE calling this function to take in the */
        /* latest values of the solution. */
        //if (grid->node[nd].resident_pe == grid->smpi->myid) {
            /* the displacement for all the nodes in the vertical list are identical and the same as the surface */
            /* the elev_factor[nd] is what gives the correct displacement for the particular node */
            ele_top = grid->node[nd].z + sw3d->displacement[nd];  // these are all 3d variables, though displacement shouldn't be
            surface = ele_top;
            vel_x_top = sw3d->vel[nd].x;
            vel_y_top = sw3d->vel[nd].y;
            sum_x = 0.;
            sum_y = 0.;
            
            ptr=ptr->next;
            while(ptr->next != NULL) {
                nd = ptr->id;
                vel_x_bot = sw3d->vel[nd].x;
                vel_y_bot = sw3d->vel[nd].y;
                ele_bot = grid->node[nd].z + sw3d->displacement[nd];
                sum_x += 0.5*(vel_x_top + vel_x_bot)*(ele_top - ele_bot);
                sum_y += 0.5*(vel_y_top + vel_y_bot)*(ele_top - ele_bot);
                ele_top = ele_bot;
                vel_x_top = vel_x_bot;
                vel_y_top = vel_y_bot;
                ptr = ptr->next;
            }
            
            depth = surface - ele_bot;
            sum_x /= depth;
            sum_y /= depth;
            ptr = grid->vertical_list[inode];
            nd = ptr->id;
            sw3d->depth_avg_vel[inode].x = sum_x;
            sw3d->depth_avg_vel[inode].y = sum_y;
            sw3d->depth[inode] = depth;
            
        //}
    }
    
}
