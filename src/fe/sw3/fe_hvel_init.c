/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Initializes SW 3D HVEL solutions and boundary conditions for a Newton iterate.
 * \author    Charlie Berger, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Gary Brown, Ph.D.
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

void fe_hvel_init(SSUPER_MODEL *sm, int imod) {
    
    int DEBUG = OFF;
    
    // aliases
    SMODEL *mod = &(sm->submodel[imod]);
    SGRID *grid = mod->grid;
    SSW_3D *sw3 = mod->sw->d3;
    
    // initialize arrays
    int i, j, istart, istring, ie;
    double mag = 0.;
    for(i = 0; i < grid->nnodes; i++) {
        sw3->vel[i].x = 0.;
        sw3->vel[i].y = 0.;
        sw3->vel[i].z = 0.;
        istart = 3 * i;
        mod->bc_mask[istart] = NO;      // x_eq
        mod->bc_mask[istart + 1] = NO;  // y_eq
        mod->bc_mask[istart + 2] = NO;  // c_eq
        
    }
    
    // assign boundary condition mask and temporarily use vel to calculate tangent vector
    svect2d_init_array(sw3->tanvec, grid->nnodes);
    if (ROTATE == ON && mod->flag.FLUX_WEIGHTED_NORMALS == OFF) {
        for (ie = 0; ie < grid->nelems2d; ie++) {
            istring = grid->elem2d[ie].string;
            if(mod->str_values[istring].flow.bc_flag == BCT_VEL_NEU || mod->str_values[istring].flow.bc_flag == BCT_DIS_NEU) {
                for (j = 0; j < grid->elem2d[ie].nnodes; j++) {
                    /* using vel to temporarily store the projected area for each node */
                    if (grid->elem2d[ie].nnodes == NDONTRI) {
                        sw3->vel[grid->elem2d[ie].nodes[j]].x += grid->elem2d[ie].djac3d * grid->elem2d[ie].nrml.x;
                        sw3->vel[grid->elem2d[ie].nodes[j]].y += grid->elem2d[ie].djac3d * grid->elem2d[ie].nrml.y;
                        sw3->vel[grid->elem2d[ie].nodes[j]].z += grid->elem2d[ie].djac3d * grid->elem2d[ie].nrml.z;
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
                        
                        
                        sw3->vel[grid->elem2d[ie].nodes[j]].x += quad_area * grid->elem2d[ie].nrml.x;
                        sw3->vel[grid->elem2d[ie].nodes[j]].y += quad_area * grid->elem2d[ie].nrml.y;
                        sw3->vel[grid->elem2d[ie].nodes[j]].z += quad_area * grid->elem2d[ie].nrml.z;
                    }
                    mod->bc_mask[3*grid->elem2d[ie].nodes[j]] = 2; // x_eq
                }
            }
        }
        
        for (i=0; i<grid->nnodes; i++) {
            mag = VECT2D_MAG(sw3->vel[i]);
            if(mag>SMALL) {
                sw3->tanvec[i].x =  sw3->vel[i].y/mag;
                sw3->tanvec[i].y = -sw3->vel[i].x/mag;
            } else {
                if (mod->bc_mask[3*i] == 2) {
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
                    mod->bc_mask[3*grid->elem2d[ie].nodes[j]] = 2; // x_eq
                }
            }
        }
        fe_weighted_normals(mod);
    }
    
    // initialize velocities and surface displacement
    for (i = 0; i < grid->nnodes; i++) {
        sw3->vel[i].x = sw3->old_vel[i].x;
        sw3->vel[i].y = sw3->old_vel[i].y;
        sw3->vel[i].z = sw3->old_vel[i].z;
        sw3->displacement[i] = sw3->old_displacement[i];
    }
    
    // enforce Dirichlet boundary conditions
    int isers = 0, ndof = 3;
    for (i = 0; i < grid->nnodes; i++) {
        istart = ndof * i;
        istring = grid->node[i].string;
        if(istring > NORMAL) {
            if(mod->str_values[istring].flow.bc_flag == BCT_VEL_DIR) {
                mod->bc_mask[istart] = YES;     // x_eq
                mod->bc_mask[istart + 1] = YES; // y_eq
                isers = mod->str_values[istring].flow.ivx;
                sw3->vel[i].x = sseries_get_value(isers, mod->series_head, 0);
                isers = mod->str_values[istring].flow.ivy;
                sw3->vel[i].y = sseries_get_value(isers, mod->series_head, 1);
            }
        }
    }
    
    // update subsurface displacement
    tl_vertical_adapt(mod->grid, sw3->displacement);
    //tl_vertical_adapt_ALE(mod->grid, sw3);
    
    // update displacement perturbations since surface dpl has changed
    tl_get_dpl_perturbation(grid, mod->perturbation, sw3->displacement, sw3->dpl_perturbation);
    
    // update depth and depth averaged velocities
    tl_calculate_depthavgvel(grid, sw3);
    
    // update pressures at vertices
    tl_calculate_pressure(mod->grid, sw3, mod->density, mod->perturbation);

#ifdef _MESSG
    comm_update_double(sw3->prs, 1, grid->smpi);
    comm_update_VECT2D(sw3->tanvec, grid->smpi);
    comm_update_VECT(sw3->old_vel, grid->smpi);
    //comm_update_double(sw3->ol_head, 1, grid->smpi);
#endif
    // update pressures at edge midpoints
    tl_find_edge_mdpt_pressure(mod->grid, sw3, mod->density, mod->perturbation);
#ifdef _MESSG
    comm_update_double(sw3->prs, 1, grid->smpi);
#endif
    
    for (i=0; i<grid->nnodes; i++) {
        if (sw3->prs[i] < 0 || sw3->prs_plus[i] < 0 || sw3->prs_minus[i] < 0) DEBUG = ON;
    }
    for (i=0; i<grid->num_midpts; i++) {
        for (j=0; j<5; j++) {
            if (grid->midpt_list[i]->value[j] < 0) DEBUG=ON;
        }
    }

#ifdef _DEBUG
    if (DEBUG) {
        svect2d_printScreen_array("hvel_init :: tangent vector", sw3->tanvec, grid->nnodes, __LINE__, __FILE__);   printf("\n");
        printScreen_dble_array("hvel_init :: displacement", sw3->displacement, grid->nnodes, __LINE__, __FILE__);   printf("\n");
        printScreen_dble_array("hvel_init :: displacement perturbation", sw3->dpl_perturbation, grid->nnodes, __LINE__, __FILE__); printf("\n");
        svect_printScreen_array("hvel_init",sw3->vel,"velocity",grid->nnodes,__LINE__,__FILE__); printf("\n");
        printScreen_dble_array("hvel_init :: node pressure", sw3->prs, grid->nnodes, __LINE__, __FILE__); printf("\n");
        printScreen_dble_array("hvel_init :: node pressure +", sw3->prs_plus, grid->nnodes, __LINE__, __FILE__); printf("\n");
        printScreen_dble_array("hvel_init :: node pressure -", sw3->prs_minus, grid->nnodes, __LINE__, __FILE__); printf("\n");
        
        printf("hvel_init :: mdpt pressures\n");
        for (i=0; i<grid->num_midpts; i++) {
            for (j=0; j<5; j++) {
                printf("%20.10e",grid->midpt_list[i]->value[j]);
            }
            printf("\n");
        }
        tl_check_all_pickets(__FILE__, __LINE__);
        tl_error("STOP\n");
    }
#endif
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Calculates flux-weighted boundary normals for the SW 3D HVEL stage.
 * \author    Gajanan Choudhary, M.S.E.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[inout]  mod (SMODEL *) pointer to the model struct
 *
 * \note Gajanan gkc - There still seems to be something tiny missing in the normal calculation that I am not able to pin down.
 * The problem can be noticed by testing tetrahedron circular slosh test with 2 levels of adaption and SRT = 10 
 * Run 2 tests with weighted normals turned off and turned on. There is 1 spot on the boundary where a inexplicable velocity spike occurs.
 * Please check that before using weighted normals
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_weighted_normals(SMODEL *mod) {
    
    int DEBUG = OFF;
    
    // aliases
    SGRID *grid = mod->grid;
    SSW_3D *sw3 = mod->sw->d3;
    ID_LIST_ITEM *ptr;
    
    // initialize arrays
    int i, j, istring, ie, icol;
    double mag = 0.;
    
    // assign boundary condition mask and temporarily use vel to calculate tangent vector
    double *columnfluxbc = tl_alloc(sizeof(double), grid->num_midpts);
    SVECT2D *normvec = tl_alloc(sizeof(SVECT2D), grid->num_midpts);
    int *nodetoedge[2]; /* Each boundary node column connected to exactly two 2d element columns */
    nodetoedge[0] = tl_alloc(sizeof(int), grid->nnodes);
    nodetoedge[1] = tl_alloc(sizeof(int), grid->nnodes);

    sarray_init_value_dbl(columnfluxbc, grid->num_midpts, 0.0);
    sarray_init_value_int(nodetoedge[0], grid->nnodes, UNSET_INT);
    sarray_init_value_int(nodetoedge[1], grid->nnodes, UNSET_INT);
    svect2d_init_array(normvec, grid->num_midpts);
    svect2d_init_array(sw3->tanvec, grid->nnodes);

    for (icol=0; icol<grid->num_midpts; icol++){
        ptr = mod->grid->sidewall_list[icol];
        while (ptr->next != NULL) {
            ie = ptr->id; /* on of the 2D elements in the sidewall. */
            istring = grid->elem2d[ie].string;
            if(mod->str_values[istring].flow.bc_flag == BCT_VEL_NEU || mod->str_values[istring].flow.bc_flag == BCT_DIS_NEU) {
                // read normal velocity/discharge
                int isers = mod->str_values[istring].flow.isigma;
                double hydro_flux_normal = sseries_get_value(isers, mod->series_head, 0); // Note: NOT taking (-) here unlike in hvel_boundary_resid.
                if (mod->str_values[istring].flow.bc_flag == BCT_DIS_NEU) {
                    hydro_flux_normal /= mod->str_values[istring].total_area;
                }
                columnfluxbc[icol] += grid->elem2d[ie].djac3d * hydro_flux_normal;
            } /* end if bc flag. */
            ptr = ptr->next;
        } /* end while ptr->next !=NULL */
    } /* end for icol= (0 to grid->num_midpts) */
    for (icol=0; icol<grid->num_midpts; icol++){
        ptr = mod->grid->sidewall_list[icol];
        while (ptr->next != NULL) {
            ie = ptr->id; /* on of the 2D elements in the sidewall. */
            normvec[icol].x = grid->elem2d[ie].nrml.x;
            normvec[icol].y = grid->elem2d[ie].nrml.y;
            for(i=0; i<NDONTRI; i++){
                if (nodetoedge[0][grid->elem2d[ie].nodes[i]] == UNSET_INT || nodetoedge[0][grid->elem2d[ie].nodes[i]] == icol){ /* First pass */
                    nodetoedge[0][grid->elem2d[ie].nodes[i]] = icol;
                }
                else if (nodetoedge[1][grid->elem2d[ie].nodes[i]] == UNSET_INT || nodetoedge[1][grid->elem2d[ie].nodes[i]] == icol){ /* First pass */
                    nodetoedge[1][grid->elem2d[ie].nodes[i]] = icol;
                }
                else {
                    tl_error("\n3 boundaries connected to 1 node. Something is wrong.\n");
                }
            }
            ptr = ptr->next;
        } /* Now nodetoedge[row][surnode] contains which surface node is connected to which 2 sidewall element columns */
    }


    for (i=0; i<grid->nnodes; i++) {
        int icol1 = nodetoedge[0][i];
        int icol2 = nodetoedge[1][i];
        if (icol1 != UNSET_INT){
            if (icol2 != UNSET_INT) {
                double q1 = columnfluxbc[icol1];
                double q2 = columnfluxbc[icol2];
                //printf("\nNode %5i : col ids (%5i, %5i) : fluxes (%15.6e, %15.6e)", i, icol1, icol2, q1, q2);
                if (fabs(q1 - q2)<SMALL){
                    sw3->vel[i].x = normvec[icol1].x + normvec[icol2].x;
                    sw3->vel[i].y = normvec[icol1].y + normvec[icol2].y;
                }
                else if (fabs(q1) > SMALL && fabs(q1 + q2)<SMALL){
                    tl_error("\nSum of non-zero fluxes close to zero; Flux weighted normals not yet supported for this case. ");
                }
                else{
                    sw3->vel[i].x = q1*normvec[icol1].x + q2*normvec[icol2].x;
                    sw3->vel[i].y = q1*normvec[icol1].y + q2*normvec[icol2].y;
                }
            }
            else {
                /* This happens in parallel runs. So we have to set the normal somehow. */
                //tl_error("\nBoundary node column connected to only 1 boundary edge. Impossible.\n");
                sw3->vel[i].x = normvec[icol1].x;
                sw3->vel[i].y = normvec[icol1].y;
            }
        }
    }
    
    for (i=0; i<grid->nnodes; i++) {
        mag = VECT2D_MAG(sw3->vel[i]);
        if(mag>SMALL) {
            sw3->tanvec[i].x =  sw3->vel[i].y/mag;
            sw3->tanvec[i].y = -sw3->vel[i].x/mag;
        } else {
            if (mod->bc_mask[3*i] == 2) {
                tl_error(">> tanvec magnitude is too small.");
            }
        }
    }
    
    normvec = tl_free(sizeof(SVECT2D), grid->num_midpts, normvec);
    columnfluxbc = tl_free(sizeof(double), grid->num_midpts, columnfluxbc);
    nodetoedge[0] = tl_free(sizeof(int), grid->nnodes, nodetoedge[0]);
    nodetoedge[1] = tl_free(sizeof(int), grid->nnodes, nodetoedge[1]);


    //svect2d_printScreen_array("hvel_init :: tangent vector", sw3->tanvec, grid->nnodes, __LINE__, __FILE__);   printf("\n");
    //exit(-2);


#ifdef _DEBUG
    if (DEBUG) {
        svect2d_printScreen_array("hvel_init :: tangent vector", sw3->tanvec, grid->nnodes, __LINE__, __FILE__);   printf("\n");
    }
#endif
}
