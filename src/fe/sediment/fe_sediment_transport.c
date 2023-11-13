/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_sediment_transport.c This file collections functions responsible for computing the next time-step */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Sets the SedLib nodal call flag.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] sedlib_call_flag (int *)    determines whether sedlib is called for a given node
 *  @param[in]     nnodes           (int)      the total number of nodes on the grid
 *  @param[in]     depth            (double *) the nodal water depths
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void get_sedlib_call_flag(int nnodes, double *depth, int *sedlib_call_flag) {
    int inode = 0;
    for (inode=0; inode < nnodes; inode++) {
        sedlib_call_flag[inode] = 1;
        if (depth[inode] < 0.001) sedlib_call_flag[inode] = 0;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     The FE engine for computing the next time-step of the sediment transport model.
 *             Returns NO if the nonlinear step fails, YES if successful.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in,out]  mod  (SMODEL *) a pointer to an AdH model where depth should be updated
 * @param[in]      ndim (int)      the dimensionality of the grid
 *
 * \details   Solves the following weak, discrete diffusive wave equation: \n
 * \f{eqnarray*}{
 * \weakTrns{2d}{e}{h}{\phidd{i}}{\ctbjh}{i}{bl}{j} \\
 * \weakTrns{2d}{e}{h}{\phidd{i}}{\cjh \, \depth{h}}{i}{sl}{j} \\
 * \weakTrns{3d}{e}{h}{\phiddd{i}}{\cjh}{i}{sl}{j}
 * \f}
 *
 * \note CJT\:: handles 2D bed load and 2D/3D suspended load transport
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int fe_sediment_transport(SMODEL *mod, int ndim) {

#ifdef _SEDLIB

#ifdef _DEBUG
    assert(ndim == 2 || ndim == 3);
    if (ndim == 2) assert(mod->sw->d2);
    if (ndim == 3) assert(mod->sw->d3);
#endif
    
    int i, j, inode=0, ie=0, ised=0, ilayer=0;
    
    // aliases
    SSED *sed = mod->sed;
    SGRID *grid = mod->grid;
    SSW_2D *sw2d = mod->sw->d2;
    SSW_3D *sw3d = mod->sw->d3;
    int nnodes = grid->nnodes;
    int nnodes_bed = grid->nnodes_bed;
    int nnodes_bed_prev = grid->nnodes_bed_prev;
#ifdef _MESG
    SMPI *smpi = mod->grid->smpi;
#endif
    
    // sedlib call flag
    int sedlib_call_flag[grid->nnodes_bed];
    if (ndim == 2) {
        get_sedlib_call_flag(nnodes_bed, sw2d->head, sed->call_flag);
    } else if (ndim == 3) {
        get_sedlib_call_flag(nnodes_bed, sw3d->depth, sed->call_flag);
    }
#ifdef _MESSG  /* also provides barrier */
    comm_update_int_surf(sed->call_flag, 1, mod->grid->smpi);
#endif

    // initialize variables over all nodes
    for (inode=0; inode<nnodes; inode++) {
        for (ised=0; ised < sed->nsed; ised++) {
            sed->susload[ised].older_c[inode] = sed->susload[ised].old_c[inode];
            sed->susload[ised].old_c[inode] = sed->susload[ised].c[inode];
            sed->susload[ised].source[inode] = 0.0;
        }
    }

    // initialize variables over bed nodes
    for (inode=0; inode<nnodes_bed; inode++) {
        for (ised=0; ised < sed->nsed; ised++) {
            sed->bedload[ised].older_c[inode] = sed->bedload[ised].old_c[inode];
            sed->bedload[ised].old_c[inode] = sed->bedload[ised].c[inode];
            sed->bedload[ised].older_thick[inode] = sed->bedload[ised].old_thick[inode];
            sed->bedload[ised].old_thick[inode] = sed->bedload[ised].thick[inode];
            sed->bedload[ised].source[inode] = 0.0;
            sed->susload[ised].rouse_coef[inode] = 1.0;
            sed->old_bed_displacement[inode] = sed->bed_displacement[inode];
            sed->old_as_ceiling[inode] = sed->as_ceiling[inode];
        }
    }

    // if running in 2d initialize transport corrections
    if (ndim == 2) {
        for (inode=0; inode<nnodes_bed; inode++) {
            for (ised=0; ised < sed->nsed; ised++) {
                sed->susload[ised].mfcf[inode] = 1.0;
                svect2d_init(&(sed->susload[ised].vcf[inode]));
                svect2d_init(&(sed->susload[ised].vor_vel[inode]));
            }
        }
    } 

    // get bathymetry gradients
    if (ndim == 2) {
        tl_nodal_grad_avg_2d(grid, mod->str_values, sw2d->bed_elevation, sed->bed_gradient);
    } else if (ndim == 3) {
        tl_nodal_grad_avg_3d(grid, mod->str_values, sw3d->bed_elevation, sed->bed_gradient);
    }

    // initialize active and bed layers
    sbed_layer_copy(sed->old_active_layer, sed->active_layer, 1, sed->nsed, nnodes_bed);
    sbed_layer_copy(sed->old_bed_layer, sed->bed_layer, sed->nlayers, sed->nsed, nnodes_bed);

    // call sediment bed source/sink wrapper
    for (inode=0; inode<nnodes_bed; inode++) {
        if (sed->call_flag[inode]) {
        
          /* SEDLIB BED LINK *******************************************************************************/
          sedlib_link_bed(mod, inode, ndim, 0, sed->bed_gradient[inode]);
          /*************************************************************************************************/

        }
    }

    // are there any material types with the bed displacement turned off (have to think about this for 3d)
        for (ie=0; ie<grid->nelems2d; ie++) {
            if (mod->mat[grid->elem2d[ie].mat].sw->bed_disp_flag == NO) {
                for (inode=0; inode<grid->elem2d[ie].nnodes; inode++) {
                    if (ndim == 3) {
                        sed->no_displacement_flag[grid->nodeID_3d_to_2d_bed[grid->elem2d[ie].nodes[inode]]] = 1;
                    } else {
                        sed->no_displacement_flag[grid->elem2d[ie].nodes[inode]] = 1;
                    }
                }
            }
        }
        for (i=0; i<nnodes_bed; i++) {
            if (sed->no_displacement_flag[i]) {
                sed->as_ceiling[i] += sed->old_bed_displacement[i] - sed->bed_displacement[i];
                sed->bed_displacement[i] = sed->old_bed_displacement[i];
            }
        }
#ifdef _MESSG
       comm_update_int_surf(sed->no_displacement_flag, 1, mod->grid->smpi);
#endif


    /* calculate transport factors */
    for (inode=0; inode<nnodes_bed; inode++) {
        if (sed->call_flag[inode]) {

            /* SEDLIB TRANSPORT LINK *******************************************************************************/
            if (ndim == 2 && mod->flag.VORTICITY) {
                //sed_vorticity_velocity_components(inode, &tvel);
            }
            sedlib_link_transport(mod, inode, ndim, 0, sed->bed_gradient[inode]);
            /*******************************************************************************************************/  
        } else {
            for (ised=0; ised < sed->nsed; ised++) {
              sed->bedload[ised].thick[inode] = 1.0;
              sed->susload[ised].mfcf[inode] = 1.0;
              svect2d_init(&(sed->bedload[ised].v[inode]));
              svect2d_init(&(sed->bedload[ised].flux[inode]));
            }
        }
    }
    //sed_sediment_diversion_adjustment();

    /***********************************************************************************/
    /***********************************************************************************/
    /* newton solve ********************************************************************/
    mod->nsys = 1;
    mod->nsys_sq = 1;
    mod->solver_info.refresh = YES;
    mod->solver_info.LINEAR_PROBLEM = YES;
   
    for (mod->ised=0; mod->ised < mod->sed->nsed; mod->ised++) {
        
        // suspended transport solve
        mod->solver_info.PRN_NEWTON = SLT; 
        if (ndim==2) {
            if (fe_newton(mod, mod->grid, fe_sw2_transport_init, fe_sw2_transport_update, fe_sw2_transport_resid, fe_sw2_transport_load, fe_sw2_transport_inc) == NO) {
                return (NO);
            }
        } else {
            if (fe_newton(mod, mod->grid, fe_sw3_transport_init, fe_sw3_transport_update, fe_sw3_transport_resid, fe_sw3_transport_load, fe_sw3_transport_inc) == NO) {
                return (NO);
            }
        }

        // bedload transport solve
        if (BEDLOAD && mod->sed->grain[mod->ised].type == SND) {
            mod->solver_info.PRN_NEWTON = BLT;
            if (fe_newton(mod, mod->grid, fe_sw2_bedload_init, fe_sw2_bedload_update, fe_sw2_bedload_resid, fe_sw2_bedload_load, fe_sw2_bedload_inc) == NO) {
                return (NO);
            }
        }

    }

    /***********************************************************************************/
    /***********************************************************************************/

    /* calculate sediment bed fluxes */
    svect2d_init_array(sed->bedload_vector, mod->grid->nnodes_bed);
    for (inode=0; inode <nnodes_bed; inode++) {
        for (ised=0; ised < sed->nsed; ised++) {
            sed->bedload_vector[inode].x += sed->bedload[ised].v[inode].x * sed->bedload[ised].c[inode] * sed->grain[ised].reference_c * sed->bedload[ised].thick[inode] * mod->density;
            sed->bedload_vector[inode].y += sed->bedload[ised].v[inode].y * sed->bedload[ised].c[inode] * sed->grain[ised].reference_c * sed->bedload[ised].thick[inode] * mod->density;
        }
    }

    /* calculate sediment sus fluxes */
    svect2d_init_array(sed->susload_vector, nnodes);
    if (ndim == 2) {
        for (inode=0; inode < nnodes; inode++) {
            for (ised=0; ised < sed->nsed; ised++) {
                sed->susload_vector[inode].x += ((sw2d->vel[inode].x + sed->susload[ised].vcf[inode].x) * sed->susload[ised].mfcf[inode]) * sed->susload[ised].c[inode] * sed->grain[ised].reference_c * sw2d->head[inode] * mod->density;
                sed->susload_vector[inode].y += ((sw2d->vel[inode].y + sed->susload[ised].vcf[inode].y) * sed->susload[ised].mfcf[inode]) * sed->susload[ised].c[inode] * sed->grain[ised].reference_c * sw2d->head[inode] * mod->density;
            }
        }
    } else if (ndim == 3) { // really this should be a 3d vector I think ...
        for (inode=0; inode < nnodes; inode++) {
            for (ised=0; ised < sed->nsed; ised++) {
                sed->susload_vector[inode].x += sw3d->vel[inode].x * sed->susload[ised].c[inode] * sed->grain[ised].reference_c * sw3d->depth[inode] * mod->density;
                sed->susload_vector[inode].y += sw3d->vel[inode].y * sed->susload[ised].c[inode] * sed->grain[ised].reference_c * sw3d->depth[inode] * mod->density;
            }
        }
    }

    // update shallow water bed displacement and elevation for two way coupling
    if (ndim == 2) {
        for (inode=0; inode < nnodes_bed; inode++) {
            sw2d->bed_displacement[inode] = sed->bed_displacement[inode];
            sw2d->bed_elevation[inode] = grid->node[inode].z + sed->bed_displacement[inode];
        }
        ssediment_calculate_2d_error(mod->sed, mod->sw->d2, mod->grid, mod->mat, mod->dt);
    } else if (ndim == 3) {

        sarray_copy_dbl(sw3d->older_bed_displacement, sw3d->old_bed_displacement, nnodes);
        sarray_copy_dbl(sw3d->old_bed_displacement, sw3d->bed_displacement, nnodes);

        int gnode = UNSET_INT;
        for (inode=0; inode < nnodes_bed; inode++) {
            gnode = grid->nodeID_2d_to_3d_bed[inode];
            sw3d->bed_displacement[gnode] = sed->bed_displacement[inode];       // bed displacement is a 3d variable
            sw3d->bed_elevation[inode] = grid->node[gnode].z + sed->bed_displacement[inode];  // bed elevation is a 2d variable
        }
        ssediment_calculate_3d_error(mod->sed, mod->sw->d3, mod->grid, mod->mat, mod->dt);
    }

    /* reset to defaults*/
    mod->solver_info.PRN_NEWTON = OFF;
#endif

    return (YES);
}

