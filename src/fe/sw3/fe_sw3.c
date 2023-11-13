/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     The FE engine for computing the next time-step of the 3D SW model. Returns
 *            NO if the nonlinear step fails, YES if successful.
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[in,out] mod (SMODEL *) a pointer to a 3D SW wave model struct
 *
 * \details   Solves the following weak, discrete SW 3D equations:
 * \nSTEP 1: - (h,u,v): \n
 * \f{eqnarray*}{
 *   \weakSWDaContReducedKinematic{i} \\
 *   \weakSWMxDDD{e}{i}{h} \\
 *   \weakSWMxDDD{e}{i}{h}
 * \f}
 * \nSTEP2 - Either of the following (w): \n
 * \f{eqnarray*}{
 *    \weakSWContSurDDD{e}{i}{h} \\
 *    \weakSWContBedDDD{e}{i}{h}
 * \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int fe_sw3(SMODEL *mod) {
    
    assert(mod->sw->d3);
    assert(mod->grid);
    
    // alias
    SSW_3D *sw3d = mod->sw->d3;
    int nnodes = mod->grid->nnodes;
    
    // prep solutions
    mod->old_dt = mod->dt;
    int i, j, inode = 0;
    for (inode = 0; inode < nnodes; inode++) {
        sw3d->older_displacement[inode]  = sw3d->old_displacement[inode];
        sw3d->older_vel[inode].x = sw3d->old_vel[inode].x;
        sw3d->older_vel[inode].y = sw3d->old_vel[inode].y;
        sw3d->older_vel[inode].z = sw3d->old_vel[inode].z;
        sw3d->old_displacement[inode]  = sw3d->displacement[inode];
        sw3d->old_vel[inode].x = sw3d->vel[inode].x;
        sw3d->old_vel[inode].y = sw3d->vel[inode].y;
        sw3d->old_vel[inode].z = sw3d->vel[inode].z;
        sw3d->density[inode] = mod->density;
        sw3d->old_grid_speed[inode] = sw3d->grid_speed[inode];
    }
    
    // reallocate stored elemental residuals
    if ((mod->flag.ADAPTED_THE_GRID == YES) && (sw3d->elem_rhs_realloc == 0)) {
        for (i=0; i<MAX_NNODES_ON_ELEM3D; i++) {
            sw3d->elem_rhs_supg_dacont[i] = (double *) tl_realloc(sizeof(double), mod->grid->nelems3d, mod->grid->nelems3d_old,sw3d->elem_rhs_supg_dacont[i]);
            sw3d->elem_rhs_supg_cont[i] = (double *) tl_realloc(sizeof(double), mod->grid->nelems3d, mod->grid->nelems3d_old, sw3d->elem_rhs_supg_cont[i]);
            for (j=0; j<mod->grid->nelems3d; j++) {
                sw3d->elem_rhs_supg_dacont[i][j] = 0.;
                sw3d->elem_rhs_supg_cont[i][j] = 0.;
            }
        }
        sw3d->elem_rhs_realloc = 1;
    } else {
        for (i=0; i<MAX_NNODES_ON_ELEM3D; i++) {
            for (j=0; j<mod->grid->nelems3d; j++) {
                sw3d->elem_rhs_supg_dacont[i][j] = 0.;
                sw3d->elem_rhs_supg_cont[i][j] = 0.;
            }
        }
    }
    
    // baroclinic SW 3D density
    switch(mod->flag.BAROCLINIC) {
        case 1:
            tl_density_calculator_metric(mod->density, NULL, 1., mod->con[mod->salinity_id].concentration, mod->con[mod->salinity_id].property[0], nnodes, sw3d->density, mod->flag.EOS);
            break;
        case 10:
            tl_density_calculator_metric(mod->density, mod->con[mod->temperature_id].concentration, 1., NULL, mod->con[mod->temperature_id].property[0], nnodes, sw3d->density, mod->flag.EOS);
            break;
        case 11:
            tl_density_calculator_metric(mod->density, mod->con[mod->temperature_id].concentration, 1., mod->con[mod->salinity_id].concentration, mod->con[mod->salinity_id].property[0], nnodes, sw3d->density, mod->flag.EOS);
            break;
    }
    
    // CJT :: calculate nodal bathymetry gradients for 2D/3D DB coupling
    // CJT :: this is only an approximation, so w_bed isn't 100% accurate like w_surface
    int nd1, nd2, nd3, ie;
    double elem_z[NDONTRI];
    SVECT2D grad_z;
    double total_weight[mod->grid->nnodes];
    svect2d_init_array(mod->sw->d3->grad_bed,mod->grid->nnodes);
    sarray_init_dbl(total_weight,mod->grid->nnodes);
    for (ie=0; ie<mod->grid->nelems2d; ie++) {
        if (mod->grid->elem2d[ie].bflag == 1) {
            nd1 = mod->grid->elem2d[ie].nodes[0];
            nd2 = mod->grid->elem2d[ie].nodes[1];
            nd3 = mod->grid->elem2d[ie].nodes[2];
            elem_z[0] = mod->grid->node[nd1].z;
            elem_z[1] = mod->grid->node[nd2].z;
            elem_z[2] = mod->grid->node[nd3].z;
            GRAD2D_FUNC(mod->grid->elem2d[ie].grad_shp, elem_z, grad_z);
            total_weight[nd1] += fabs(mod->grid->elem2d[ie].djac);
            total_weight[nd2] += fabs(mod->grid->elem2d[ie].djac);
            total_weight[nd3] += fabs(mod->grid->elem2d[ie].djac);
            mod->sw->d3->grad_bed[nd1].x += grad_z.x * fabs(mod->grid->elem2d[ie].djac);
            mod->sw->d3->grad_bed[nd2].x += grad_z.x * fabs(mod->grid->elem2d[ie].djac);
            mod->sw->d3->grad_bed[nd3].x += grad_z.x * fabs(mod->grid->elem2d[ie].djac);
            mod->sw->d3->grad_bed[nd1].y += grad_z.y * fabs(mod->grid->elem2d[ie].djac);
            mod->sw->d3->grad_bed[nd2].y += grad_z.y * fabs(mod->grid->elem2d[ie].djac);
            mod->sw->d3->grad_bed[nd3].y += grad_z.y * fabs(mod->grid->elem2d[ie].djac);
        }
    }
    for (i=0; i<mod->grid->nnodes; i++) {
        if (total_weight[i] > SMALL) {
            mod->sw->d3->grad_bed[i].x /= fabs(total_weight[i]);
            mod->sw->d3->grad_bed[i].y /= fabs(total_weight[i]);
        }
    }
    
    mod->nsys = 3;
    mod->nsys_sq = 9;
    
    return (YES);
}
