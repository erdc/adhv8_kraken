/* Finds the normalized pressure at midpoint along an edge */
/* the pressure is normalized by the gravity and the reference density */
// cjt :: note :: for prism grids, subsurface midpoints can fall on quadrilateral and triangular faces,
//                since adjacent elements of different layers have a tet transition
// cjt :: note :: for tet grids, all subsurface midpoints fall on triangular faces

#include "global_header.h"


// ********************************************************************************************** //
// ********************************************************************************************** //
// Used for both Tets and Prisms
void tl_find_edge_mdpt_pressure_vertical(SGRID *grid, SSW_3D *sw3d, double density, double perturbation) {
    
    MIDPT_LIST_ITEM *myptr;
    int i;
    int top_node, bot_node;
    double pressure, press_above;
    double pressure1, press_above1; /* these are the pressures associated with the positive perturbation of the surface node 1.  This is the node 1 of the top edge. */
    double pressure2, press_above2; /* same as above but for node 2 */
    double pressure3, press_above3; /* these are the pressures associated with the negative perturbation of the surface node 1.  This is the node 1 of the top edge. */
    double pressure4, press_above4; /* same as above but for node 2 */
    double top_z = 0., bot_z = 0.;
    
    /* now calculate edge midpoints below vertices (verticles) */
    for (i=0; i<grid->num_midpts; i++)  {
        
        /* myptr is a pointer to the beginning of one of the linked lists */
        myptr = grid->midpt_list[i];
        
        /* Skip the non-vertical segments; in this loop we are dealing with the  pressure at midpoints of the vertical edges */
        if(myptr->vertical == 0) {
            continue;
        }
        
        press_above=0.;
        press_above1=0.;
        press_above3=0.;
        while(myptr->next != NULL) {
            top_node = myptr->node1;
            bot_node = myptr->node2;
            
            // cjt :: sometimes the bottom face is numbered first in the prism declaration
            //printf("before || (top_z - bot_z): %20.10f\n",(top_z - bot_z));
            //if (grid->node[top_node].z < grid->node[bot_node].z) {
            //    top_node = myptr->node2;
            //    bot_node = myptr->node1;
            //}
            //printf("after || (top_z - bot_z): %20.10f\n",(top_z - bot_z));
            
            // integrate only to verticle edge midpoint
            top_z = grid->node[top_node].z + sw3d->displacement[top_node];
            bot_z = grid->node[bot_node].z + sw3d->displacement[bot_node];
            pressure  = press_above + (top_z - bot_z) * one_8 * (3.*sw3d->density[top_node] + sw3d->density[bot_node])/density;
            
            top_z = grid->node[top_node].z + (sw3d->displacement[top_node] + sw3d->dpl_perturbation[top_node]);
            bot_z = grid->node[bot_node].z + (sw3d->displacement[bot_node] + sw3d->dpl_perturbation[bot_node]);
            pressure1 = press_above1 + (top_z - bot_z) * one_8 * (3.*sw3d->density[top_node] + sw3d->density[bot_node])/density; // from + perturb of node 1
            
            top_z = grid->node[top_node].z + (sw3d->displacement[top_node] - sw3d->dpl_perturbation[top_node]);
            bot_z = grid->node[bot_node].z + (sw3d->displacement[bot_node] - sw3d->dpl_perturbation[bot_node]);
            pressure3 = press_above3 + (top_z - bot_z) * one_8 * (3.*sw3d->density[top_node] + sw3d->density[bot_node])/density; // from - perturb of node 1
            
            /* since we are on a vertical edge, the "left" and "right" surface nodes are the same,*/
            pressure2 = pressure1; // from + perturb of node 2
            pressure4 = pressure3; // from - perturb of node 2
            
            myptr->value[0] = pressure;
            myptr->value[1] = pressure1;
            myptr->value[2] = pressure2;
            myptr->value[3] = pressure3;
            myptr->value[4] = pressure4;
            
            // integral all the way down the edge
            top_z = grid->node[top_node].z + sw3d->displacement[top_node];
            bot_z = grid->node[bot_node].z + sw3d->displacement[bot_node];
            pressure  = press_above + (top_z - bot_z) * 0.5 * (sw3d->density[top_node] + sw3d->density[bot_node])/density;
            
            top_z = grid->node[top_node].z + (sw3d->displacement[top_node] + sw3d->dpl_perturbation[top_node]);
            bot_z = grid->node[bot_node].z + (sw3d->displacement[bot_node] + sw3d->dpl_perturbation[bot_node]);
            pressure1 = press_above1 + (top_z - bot_z) * 0.5 * (sw3d->density[top_node] + sw3d->density[bot_node])/density; // from + perturb of node 1
            
            top_z = grid->node[top_node].z + (sw3d->displacement[top_node] - sw3d->dpl_perturbation[top_node]);
            bot_z = grid->node[bot_node].z + (sw3d->displacement[bot_node] - sw3d->dpl_perturbation[bot_node]);
            pressure3 = press_above3 + (top_z - bot_z) * 0.5 * (sw3d->density[top_node] + sw3d->density[bot_node])/density; // from - perturb of node 1
            
            myptr = myptr->next;
            press_above=pressure;
            press_above1=pressure1;
            press_above3=pressure3;
            
        }
    }
}

// ********************************************************************************************** //
// ********************************************************************************************** //

void tl_find_edge_mdpt_pressure(SGRID *grid, SSW_3D *sw3d, double density, double perturbation) {
    
    MIDPT_LIST_ITEM *myptr;
    int new_nd1, new_nd2;
    int old_nd1, old_nd2;
    int node_a = -1, node_b = -1, node_c = -1, node_d = -1, isTet = TRUE;  /* cjt :: initialized */
    int i;
    int node_flag;
    int old_sur_node1, old_sur_node2;
    int sur_nd1, sur_nd2;
    int top_node, bot_node;
    double pressure, press_above;
    double pressure1, press_above1; /* these are the pressures associated with the positive perturbation of the surface node 1.  This is the node 1 of the top edge. */
    double pressure2, press_above2; /* same as above but for node 2 */
    double pressure3, press_above3; /* these are the pressures associated with the negative perturbation of the surface node 1.  This is the node 1 of the top edge. */
    double pressure4, press_above4; /* same as above but for node 2 */
    double z_a,z_b,z_c,z_d,rho_a,rho_b,rho_c,rho_d;
    double z_a_t, z_b_t, z_c_t, z_d_t;
    double z_a_m, z_b_m, z_c_m, z_d_m;
    double top_z = 0., bot_z = 0.;
    double perturb_scaled = 0., dum = 0.; // the dpl scaled perturbations (cjt)
    
    /* Loop over the midpoint lists dropped from surface midpoints first (non-verticle) */
    for (i=0; i<grid->num_midpts; i++) {
        
        /* myptr is a pointer to the beginning of one of the linked lists */
        myptr = grid->midpt_list[i];
        
        /* Skip the vertical segments; they are easier and handled afterwards */
        if(myptr->vertical == 1) {
            continue;
        }
        
        old_nd1=myptr->node1;
        old_nd2=myptr->node2;
        myptr->value[0] = 0.;   /* set the midedge pressure at the surface to zero :: should actually be set to atmospheric (come back to this) */
        press_above = 0.;       /* same as note above */
        myptr->value[1] = 0.;   /* the midedge pressure at the surface when surface node 1 is perturbed */
        press_above1 = 0.;
        myptr->value[2] = 0.;
        press_above2 = 0.;
        myptr->value[3] = 0.;
        press_above3 = 0.;
        myptr->value[4] = 0.;
        press_above4 = 0.;
        old_sur_node1 = old_nd1;
        old_sur_node2 = old_nd2;
        
        
        /* The beginning of the linked list points to the surface edge, the traversed down the midpt column */
        /* This advances the pointer to fill in pressures on non-verticle subsurface edges */
        if(myptr->next != NULL) {
            myptr = myptr->next;
        }
        
        /* Now loop through the list and fill in pressures on edges */
        while(myptr->next != NULL) {
            
            /* we have the indices of the two endpoints of the edge */
            /* (note: these are sorted so that new_nd1 < new_nd2)  */
            new_nd1 = myptr->node1;
            new_nd2 = myptr->node2;
            
            /* find out which end actually drops down to the next node as we go down the edges */
            if(new_nd1 == old_nd1) {
                // old_nd1 ---- old_nd2
                //
                // a--x--b
                //  \    |
                //   \   |
                //    x  |
                //     \ |
                //       c
                node_a=old_nd1; // new_nd1
                node_b=old_nd2;
                node_c=new_nd2;
                node_flag=2;
            } else if(new_nd2 == old_nd2) {
                // old_nd1 ---- old_nd2
                //
                // b -- a
                // |   /
                // |  /
                // | /
                // c
                node_a=old_nd2; // new_nd2
                node_b=old_nd1;
                node_c=new_nd1;
                node_flag=1;
            }
            else {
                // CJT :: then this edge is from a triangle prism, since no non-vertical edges share a node
                // Because midpoints are ordered down a column, we know that old are on top and new nodes are on bottom
                // However, we must make sure the nodes on a segment are order correctly for 2D quadrilateral basis interpolation
                // Edges are define in cosntant.h. They are defined as : nd_on_QuadEdge[4][2] = { {0, 1}, {1, 2}, {2, 3}, {3, 0} };
                // This is counterclock-wise around the quad.  If the nodes on a 3d element face are given in a counter-clockwise direction,
                // Then the first node on a segment will always be new/old_nd1.
                // Since the extrusion code always numbers the 3D element top and bottom faces in counter clock-wise fashion, we're good ... I think
                // old_nd1 -------- old_nd2
                //    |                |
                //    |                |
                //    |                |
                // new_nd1 -------- new_nd2
                //
                //(-1,1)           (1,1)
                // rho_d           rho_c
                // 4(4) O----X----O 3(3) local node numbering ()
                //      |   mp1   |
                //      |         |
                //      |         |
                //      |         | mp3
                //      X   ie=1  X
                //      |         |
                //      |         |
                //      |         |
                //      |   mp2   |
                // 1(1) O----X----O 2(2) local node numbering ()
                // rho_a           rho_b
                //(-1,-1)          (1,-1)
                //
                //   0------2
                //   |\    /|
                //   | \  / |
                //   |  \/  |
                //   |  1   |
                //   3--|---5
                //    \ |  /
                //     \| /
                //      \/
                //      4
                isTet = FALSE;
                node_a = new_nd1;
                node_b = new_nd2;
                node_c = old_nd2;
                node_d = old_nd1;
                
                if (grid->node[node_a].z > grid->node[node_d].z) {
                    node_a = old_nd1;
                    node_b = old_nd2;
                    node_c = new_nd2;
                    node_d = new_nd1;
                }
                // make sure these are not on the same horizontal location
                assert((pow(grid->node[node_a].x - grid->node[node_b].x,2) + pow(grid->node[node_a].y - grid->node[node_b].y,2)) > 1e-6);
                // make sure d is above a
                assert(grid->node[node_a].z < grid->node[node_d].z);
            }
            
            // non-perturbed pressure calculation
            z_a = grid->node[node_a].z + sw3d->displacement[node_a];
            z_b = grid->node[node_b].z + sw3d->displacement[node_b];
            z_c = grid->node[node_c].z + sw3d->displacement[node_c];
            rho_a = sw3d->density[node_a] / density;
            rho_b = sw3d->density[node_b] / density;
            rho_c = sw3d->density[node_c] / density;
            
            if (isTet == TRUE) {
                pressure = 0.5*(z_b - z_c)*(0.5*rho_a + 0.25*(rho_b + rho_c)) + press_above;
                
                // midpoint below surface midpoint (note: a perturbed dpl for node a will never perturb pressure)
                if(node_flag==2) {
                    z_b_t = z_b + sw3d->dpl_perturbation[node_b];
                    z_c_t = z_c + sw3d->dpl_perturbation[node_c];
                    z_b_m = z_b - sw3d->dpl_perturbation[node_b];
                    z_c_m = z_c - sw3d->dpl_perturbation[node_c];
                    pressure2 = 0.5*(z_b_t - z_c_t)*(0.5*rho_a + 0.25*(rho_b + rho_c)) + press_above2; // from + perturb of node 2
                    pressure4 = 0.5*(z_b_m - z_c_m)*(0.5*rho_a + 0.25*(rho_b + rho_c)) + press_above4; // from - perturb of node 2
                    z_b_t = z_b;
                    z_b_m = z_b;
                    z_c_t = z_c;
                    z_c_m = z_c;
                    pressure1 = 0.5*(z_b_t - z_c_t)*(0.5*rho_a + 0.25*(rho_b + rho_c)) + press_above1; // from + perturb of node 1
                    pressure3 = 0.5*(z_b_m - z_c_m)*(0.5*rho_a + 0.25*(rho_b + rho_c)) + press_above3; // from - perturb of node 1
                } else {
                    z_b_t = z_b + sw3d->dpl_perturbation[node_b];
                    z_c_t = z_c + sw3d->dpl_perturbation[node_c];
                    z_b_m = z_b - sw3d->dpl_perturbation[node_b];
                    z_c_m = z_c - sw3d->dpl_perturbation[node_c];
                    pressure1 = 0.5*(z_b_t - z_c_t)*(0.5*rho_a + 0.25*(rho_b + rho_c)) + press_above1; // from + perturb of node 1
                    pressure3 = 0.5*(z_b_m - z_c_m)*(0.5*rho_a + 0.25*(rho_b + rho_c)) + press_above3; // from - perturb of node 1
                    z_b_t = z_b;
                    z_b_m = z_b;
                    z_c_t = z_c;
                    z_c_m = z_c;
                    pressure2 = 0.5*(z_b_t - z_c_t)*(0.5*rho_a + 0.25*(rho_b + rho_c)) + press_above2; // from + perturb of node 2
                    pressure4 = 0.5*(z_b_m - z_c_m)*(0.5*rho_a + 0.25*(rho_b + rho_c)) + press_above4; // from - perturb of node 2
                }
            } else {
                
                z_d = grid->node[node_d].z + sw3d->displacement[node_d];
                rho_d = sw3d->density[node_d] / density;

                pressure = -one_8*(rho_a + rho_b + rho_c + rho_d)*(z_a + z_b - z_c - z_d) + press_above;
                
                // if the column including nodes a and d is perturbed
                z_a_t = z_a + sw3d->dpl_perturbation[node_a];
                z_d_t = z_d + sw3d->dpl_perturbation[node_d];
                z_a_m = z_a - sw3d->dpl_perturbation[node_a];
                z_d_m = z_d - sw3d->dpl_perturbation[node_d];
                pressure1 = -one_8*(rho_a + rho_b + rho_c + rho_d)*(z_a_t + z_b - z_c - z_d_t) + press_above1; // from + perturb of node 1
                pressure3 = -one_8*(rho_a + rho_b + rho_c + rho_d)*(z_a_m + z_b - z_c - z_d_t) + press_above3; // from - perturb of node 1
                
                // if the column including nodes b and c is perturbed
                z_b_t = z_b + sw3d->dpl_perturbation[node_b];
                z_c_t = z_c + sw3d->dpl_perturbation[node_c];
                z_b_m = z_b - sw3d->dpl_perturbation[node_b];
                z_c_m = z_c - sw3d->dpl_perturbation[node_c];
                pressure2 = -one_8*(rho_a + rho_b + rho_c + rho_d)*(z_a + z_b_t - z_c_t - z_d) + press_above2; // from + perturb of node 2
                pressure4 = -one_8*(rho_a + rho_b + rho_c + rho_d)*(z_a + z_b_m - z_c_m - z_d) + press_above4; // from - perturb of node 2
                
            }
            
            /* Store the pressure in the entry 'value' */
            myptr->value[0] = pressure;
            myptr->value[1] = pressure1;
            myptr->value[2] = pressure2;
            myptr->value[3] = pressure3;
            myptr->value[4] = pressure4;
            
            /* BE sure to move to next entry in the linked list */
            myptr = myptr->next;
            old_nd1=new_nd1;
            old_nd2=new_nd2;
            press_above=pressure;
            press_above1=pressure1;
            press_above2=pressure2;
            press_above3=pressure3;
            press_above4=pressure4;
        }
    }
    
    // find pertubed midpoint pressures for all nodes on element vertical edges
    tl_find_edge_mdpt_pressure_vertical(grid, sw3d, density, perturbation);
}
