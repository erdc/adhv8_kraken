/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Initializes 3D Navier Stokes solutions and boundary conditions for a Newton iterate.
 * \author    Corey Trahan, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    David Smith, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[inout]  mod (SMODEL *) pointer to the model struct
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_ns3_init(SSUPER_MODEL *sm, int imod) {
 
    SMODEL *mod = &(sm->submodel[imod]);

    int DEBUG = OFF;
    
    // aliases
    SGRID *grid = mod->grid;
    SNS_3D *ns3 = mod->ns->d3;
    
    // initialize arrays
    int inode, i, j, istart, istring, ie;
    double mag = 0.;
    for(i = 0; i < grid->nnodes; i++) {
        ns3->vel[i].x = 0.;
        ns3->vel[i].y = 0.;
        ns3->vel[i].z = 0.;
        istart = 4 * i;
        mod->bc_mask[istart] = NO;      // x_eq
        mod->bc_mask[istart + 1] = NO;  // y_eq
        mod->bc_mask[istart + 2] = NO;  // z_eq
        mod->bc_mask[istart + 3] = NO;  // prs_eq
    }
    
    // assign boundary condition mask and temporarily use vel to calculate tangent vector
    svect2d_init_array(ns3->tanvec, grid->nnodes);
    if (ROTATE == ON && mod->flag.FLUX_WEIGHTED_NORMALS == OFF) {
        for (ie = 0; ie < grid->nelems2d; ie++) {
            if (grid->elem2d[ie].bflag == 2) { // only for sidewalls
                
                istring = grid->elem2d[ie].string;
                if(mod->str_values[istring].flow.bc_flag == BCT_VEL_NEU || mod->str_values[istring].flow.bc_flag == BCT_DIS_NEU) {
                    for (j = 0; j < grid->elem2d[ie].nnodes; j++) {
                        /* using vel to temporarily store the projected area for each node */
                        if (grid->elem2d[ie].nnodes == NDONTRI) {
                            ns3->vel[grid->elem2d[ie].nodes[j]].x += grid->elem2d[ie].djac3d * grid->elem2d[ie].nrml.x;
                            ns3->vel[grid->elem2d[ie].nodes[j]].y += grid->elem2d[ie].djac3d * grid->elem2d[ie].nrml.y;
                            ns3->vel[grid->elem2d[ie].nodes[j]].z += grid->elem2d[ie].djac3d * grid->elem2d[ie].nrml.z;
                        } else { // CJT :: can't use djac to weight for prisms, use areas, do this later
                            // get 2d quad element area to weight (should be calculated at t=0 and stored in elem2d.djac3d...)
                            // use cross product to get area of quadrilateral (cannot use this for djac)
                            SVECT nd[grid->elem2d[ie].nnodes];
                            for (i=0; i<grid->elem2d[ie].nnodes; i++) {
                                nd[i].x = grid->node[grid->elem2d[ie].nodes[i]].x;
                                nd[i].y = grid->node[grid->elem2d[ie].nodes[i]].y;
                                nd[i].z = grid->node[grid->elem2d[ie].nodes[i]].z;
                            }
                            SVECT v1; svect_subtract_array2(&v1, &nd[0], &nd[1], 1);
                            SVECT v2; svect_subtract_array2(&v2, &nd[0], &nd[3], 1);
                            SVECT w1; svect_subtract_array2(&w1, &nd[2], &nd[1], 1);
                            SVECT w2; svect_subtract_array2(&w2, &nd[2], &nd[3], 1);
                            double quad_area = one_2 * (svect_mag(svect_cross(v1,v2)) + svect_mag(svect_cross(w1,w2)));
                            assert(quad_area > 1e-6);
                            
                            
                            ns3->vel[grid->elem2d[ie].nodes[j]].x += quad_area * grid->elem2d[ie].nrml.x;
                            ns3->vel[grid->elem2d[ie].nodes[j]].y += quad_area * grid->elem2d[ie].nrml.y;
                            ns3->vel[grid->elem2d[ie].nodes[j]].z += quad_area * grid->elem2d[ie].nrml.z;
                        }
                        mod->bc_mask[4*grid->elem2d[ie].nodes[j]] = 2; // x_eq
                    }
                }

            }
        }
        
        for (i=0; i<grid->nnodes; i++) {
            mag = VECT2D_MAG(ns3->vel[i]);
            if(mag>SMALL) {
                ns3->tanvec[i].x =  ns3->vel[i].y/mag;
                ns3->tanvec[i].y = -ns3->vel[i].x/mag;
            } else {
                if (mod->bc_mask[4*i] == 2) {
                    tl_error(">> tanvec magnitude is too small.");
                }
            }
        }
    }
    else if (ROTATE == ON && mod->flag.FLUX_WEIGHTED_NORMALS == ON) {
        for (ie = 0; ie < grid->nelems2d; ie++) {
            istring = grid->elem2d[ie].string;
            if(mod->str_values[istring].flow.bc_flag == BCT_VEL_NEU || mod->str_values[istring].flow.bc_flag == BCT_DIS_NEU) {
                for (j = 0; j < grid->elem2d[ie].nnodes; j++) {
                    mod->bc_mask[4*grid->elem2d[ie].nodes[j]] = 2; // x_eq
                }
            }
        }
        fe_weighted_normals(mod);
    }
    
    // initialize velocities and pressure
    for(i = 0; i < grid->nnodes; i++) {
        ns3->vel[i].x = ns3->old_vel[i].x;
        ns3->vel[i].y = ns3->old_vel[i].y;
        ns3->vel[i].z = ns3->old_vel[i].z;
        ns3->prs[i]   = ns3->old_prs[i];
        ns3->displacement[i]    = ns3->old_displacement[i];
    }

    /*
    if (mod->flag.MG == ON) {
        for (i = 0; i < grid->nnodes; i++) {
            if (mod->grid->node[i].bflag == 0) {
                ns3->displacement[i]    = ns3->old_displacement[i];
            }
        }
    }
    */
    
    // enforce the Dirichlet boundary conditions
    int isers = 0, ndof = 4;
    for(i = 0; i < grid->nnodes; i++) {
        istart = ndof * i;
        istring = grid->node[i].string;
        if(istring > NORMAL) {
            if(mod->str_values[istring].flow.bc_flag == BCT_VEL_DIR) {
                mod->bc_mask[istart] = YES;
                mod->bc_mask[istart + 1] = YES;
                mod->bc_mask[istart + 2] = YES;
                isers = mod->str_values[istring].flow.ivx;
                ns3->vel[i].x = sseries_get_value(isers, mod->series_head, 0);
                isers = mod->str_values[istring].flow.ivy;
                ns3->vel[i].y = sseries_get_value(isers, mod->series_head, 0);
                isers = mod->str_values[istring].flow.ivz;
                ns3->vel[i].z = sseries_get_value(isers, mod->series_head, 0);
            }
            if(mod->str_values[istring].pressure.bc_flag == BCT_PRS_DIR) {
                mod->bc_mask[istart + 3] = YES;
                isers = mod->str_values[istring].pressure.iu_0;
                ns3->prs[i] = sseries_get_value(isers, mod->series_head, 0); // * (1./mod->density);
            }
            else if(mod->str_values[istring].pressure.bc_flag == BCT_HSP_DIR) {
                mod->bc_mask[istart + 3] = YES;
                double height = mod->str_values[istring].ref_height - grid->node[i].z;
                ns3->prs[i] = mod->str_values[istring].ref_pressure + mod->density * mod->gravity * height;
            }
        }

        //printf("u[%d]: %20.10f\n",i,ns3->vel[i].x);
    }
    //exit(-1);
}
