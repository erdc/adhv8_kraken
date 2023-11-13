/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     The FE engine for computing the next time-step of the 3D NS model. Returns
 *            NO if the nonlinear step fails, YES if successful.
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

int fe_ns3(SMODEL *mod) {

    assert(mod->ns->d3);
    assert(mod->grid);
    
    // alias
    SNS_3D *ns3d = mod->ns->d3;
    int nnodes = mod->grid->nnodes;
    
    // prep solutions
    mod->old_dt = mod->dt;
    int i, j, inode = 0;
    for (inode = 0; inode < nnodes; inode++) {
        ns3d->older_prs[inode]    = ns3d->old_prs[inode];
        ns3d->older_vel[inode].x  = ns3d->old_vel[inode].x;
        ns3d->older_vel[inode].y  = ns3d->old_vel[inode].y;
        ns3d->older_vel[inode].z  = ns3d->old_vel[inode].z;
        ns3d->old_prs[inode]   = ns3d->prs[inode];
        ns3d->old_vel[inode].x = ns3d->vel[inode].x;
        ns3d->old_vel[inode].y = ns3d->vel[inode].y;
        ns3d->old_vel[inode].z = ns3d->vel[inode].z;
        ns3d->density[inode] = mod->density;
        ns3d->older_displacement[inode]  = ns3d->old_displacement[inode];
        ns3d->old_displacement[inode]    = ns3d->displacement[inode];
    }

    /*
    if (mod->flag.MG == ON) {
        for (inode = 0; inode < nnodes; inode++) {
            if (mod->grid->node[inode].bflag == 0) {
                ns3d->older_displacement[inode]  = ns3d->old_displacement[inode];
                ns3d->old_displacement[inode]    = ns3d->displacement[inode];
            }
        }
    }
    */
    
    // reallocate stored elemental residuals
    /*
    if ((mod->flag.ADAPTED_THE_GRID == YES) && (ns3d->elem_rhs_realloc ==0)) {
        for (i=0; i<MAX_NNODES_ON_ELEM3D; i++) {
            ns3d->elem_rhs_supg_dacont[i] = (double *) tl_realloc(sizeof(double), mod->grid->nelems3d, mod->grid->nelems3d_old,ns3d->elem_rhs_supg_dacont[i]);
            ns3d->elem_rhs_supg_cont[i] = (double *) tl_realloc(sizeof(double), mod->grid->nelems3d, mod->grid->nelems3d_old, ns3d->elem_rhs_supg_cont[i]);
            for (j=0; j<mod->grid->nelems3d; j++) {
                ns3d->elem_rhs_supg_dacont[i][j] = 0.;
                ns3d->elem_rhs_supg_cont[i][j] = 0.;
            }
        }
        ns3d->elem_rhs_realloc = 1;
    } else {
        for (i=0; i<MAX_NNODES_ON_ELEM3D; i++) {
            for (j=0; j<mod->grid->nelems3d; j++) {
                ns3d->elem_rhs_supg_dacont[i][j] = 0.;
                ns3d->elem_rhs_supg_cont[i][j] = 0.;
            }
        }
    }
    */
    
    return (YES);
}
