/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Calculates the 2D and 3D transport residuals
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[in,out] mod (SMODEL *) a pointer to the diffusive wave model struct
 *
 *  \details   Solves the following weak, discrete 2D or 3D transport equation: \n
 *    \f{eqnarray*}{
 *    \weakTrnsDA{2d}{e}{h}{\phidd{i}}{\cjhDA \, \depth{h}}{i}{t}{j} \\
 *    \weakTrns{3d}{e}{h}{\phiddd{i}}{\cjh}{i}{t}{j} \\
 *  \f}
 *  \n where \f$\cjh\f$ is constituent j's concentration, \f$\cjhDA\f$ is depth-averaged constituent j's concentration \n
 *  \f$ \depth{h} \f$ is the water depth and \f$ \velc{h}{t} \f$ and \f$ \velcDA{h}{t} \f$ are water velocities and  depth-averaged velocities.
 *
 *
 * \note cjt \:: assumes the no 1D elements are used on 3D grids for now
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_transport_resid(SSUPER_MODEL *sm, int imod) {

    SMODEL *mod = &(sm->submodel[imod]);
#ifdef _DEBUG
    assert(mod->flag.TRANSPORT == ON);
    assert(mod->nsys == 1);
#endif
    
    int DEBUG = OFF;
    int DEBUG_PICKETS = OFF;
#ifdef _DEBUG
    if (mod->t_prev > DEBUG_TIME) debug.residual = ON;
    if (debug.residual == ON) {
        DEBUG = ON;
    }
#endif
    
    int i = 0, ie = UNSET_INT;
    double elem_rhs[NDONPRISM]; // allocated to max elemental nodes (triangular prisms)

#ifdef _PETSC
    PetscInt row_index, ierr;
    //PetscScalar value;
    PetscScalar values[3];
    // Set these to 0 in case solving with max_nsys=3
    values[1]=0; values[2]=0;

    int GnodeID, ifmap; // For readability when determining place in PETSc matrix
#endif
    
    // alias
    SGRID *grid = mod->grid;
    
    /* initializes the residual */
    if (mod->amICoupled == NO) {
#ifdef _PETSC
        // Zero residual
        ierr = VecZeroEntries(sm->residual);
#else
        sarray_init_dbl(sm->residual, sm->nnodes * sm->max_nsys);
#endif
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 3D ELEMENT CONTRIBUTIONS
    for (ie = 0; ie < grid->nelems3d; ie++) {
        sarray_init_dbl(elem_rhs,NDONPRISM);
        fe_3d_transport_body_resid(mod, elem_rhs, ie, 0., UNSET_INT, PERTURB_NONE, 1, DEBUG);
#ifdef _DEBUG
        if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
#endif
        for (i = 0; i < grid->elem3d[ie].nnodes; i++) {
            // add to residual if not a Dirichlet BC
            if (mod->bc_mask[grid->elem3d[ie].nodes[i]] == NO) {
#ifdef _PETSC
                //value = -elem_rhs[i];
                values[0] = -elem_rhs[i];
                GnodeID = grid->elem3d[ie].nodes[i];
                ifmap = mod->fmap[GnodeID];
                if(grid->node[GnodeID].resident_pe == grid->smpi->myid){
                    //VecSetValuesLocal(sm->residual,1,&(ifmap),&(value),ADD_VALUES);
                    VecSetValuesBlockedLocal(sm->residual,1,&(ifmap),values,ADD_VALUES);
                }
#else
                sm->residual[mod->fmap[grid->elem3d[ie].nodes[i]]] -= elem_rhs[i];
#endif
            }
        }
    }

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 2D ELEMENT CONTRIBUTIONS
    if (grid->ndim == 2) { // 2D TRANSPORT BODY CONTRIBUTIONS
        for (ie = 0; ie < grid->nelems2d; ie++) {
            if (grid->wd_flag[ie] == 0 || grid->wd_flag[ie] == 1) {
                //if (mod->str_values[grid->elem2d[ie].string].phys_flag == SW2_FLAG) {
                sarray_init_dbl(elem_rhs,NDONPRISM);
                fe_2d_transport_body_resid(mod, elem_rhs, ie, 0., UNSET_INT, PERTURB_NONE, 1, DEBUG);
#ifdef _DEBUG
                if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
#endif
                for (i = 0; i < grid->elem2d[ie].nnodes; i++) {
                    // add to residual if not a Dirichlet BC (otherwise = 0)
                    if (mod->bc_mask[grid->elem2d[ie].nodes[i]] == NO) {
#ifdef _PETSC
                        //value = -elem_rhs[i];
                        values[0] = -elem_rhs[i];
                        GnodeID = grid->elem2d[ie].nodes[i];
                        ifmap = mod->fmap[GnodeID];
                        if(grid->node[GnodeID].resident_pe == grid->smpi->myid){
                            //VecSetValuesLocal(sm->residual,1,&(ifmap),&(value),ADD_VALUES);
                            VecSetValuesBlockedLocal(sm->residual,1,&(ifmap),values,ADD_VALUES);
                        }
#else
                        sm->residual[mod->fmap[grid->elem2d[ie].nodes[i]]] -= elem_rhs[i];
#endif
                    }
                }
                //}
            }
        }
    } else { // 3D TRANSPORT BOUNDARY CONTRIBUTIONS
        for (ie = 0; ie < grid->nelems2d; ie++) {
            sarray_init_dbl(elem_rhs,NDONPRISM);
            fe_3d_transport_boundary_resid(mod, elem_rhs, ie, 0., UNSET_INT, PERTURB_NONE, 1, DEBUG);
#ifdef _DEBUG
            if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
#endif
            for (i = 0; i < grid->elem2d[ie].nnodes; i++) {
                // add to residual if not a Dirichlet BC (otherwise = 0)
                if (mod->bc_mask[grid->elem2d[ie].nodes[i]] == NO) {
#ifdef _PETSC
                    //value = -elem_rhs[i];
                    values[0] = -elem_rhs[i];
                    GnodeID = grid->elem2d[ie].nodes[i];
                    ifmap = mod->fmap[GnodeID];
                    if(grid->node[GnodeID].resident_pe == grid->smpi->myid){
                        //VecSetValuesLocal(sm->residual,1,&(ifmap),&(value),ADD_VALUES);
                        VecSetValuesBlockedLocal(sm->residual,1,&(ifmap),values,ADD_VALUES);
                    }
#else
                    sm->residual[mod->fmap[grid->elem2d[ie].nodes[i]]] -= elem_rhs[i];
#endif
                }
            }
        }
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 1D ELEMENT CONTRIBUTIONS
    for (ie = 0; ie < grid->nelems1d; ie++) {
        //if (mod->str_values[grid->elem1d[ie].string].phys_flag == SW2_FLAG) {
        sarray_init_dbl(elem_rhs,NDONPRISM);
        fe_2d_transport_boundary_resid(mod, elem_rhs, ie, 0., UNSET_INT, PERTURB_NONE, 1, DEBUG);
#ifdef _DEBUG
        if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
#endif
        for (i = 0; i < grid->elem1d[ie].nnodes; i++) {
            // add to residual if not a Dirichlet BC (otherwise = 0)
            if (mod->bc_mask[grid->elem1d[ie].nodes[i]] == NO) {
#ifdef _PETSC
                //value = -elem_rhs[i];
                values[0] = -elem_rhs[i];
                GnodeID = grid->elem1d[ie].nodes[i];
                ifmap = mod->fmap[GnodeID];
                row_index = ifmap+sm->Istart;
                if(grid->node[GnodeID].resident_pe == grid->smpi->myid){
                    //VecSetValuesLocal(sm->residual,1,&(ifmap),&(value),ADD_VALUES);
                    VecSetValuesBlockedLocal(sm->residual,1,&(ifmap),values,ADD_VALUES);
                }
#else
                sm->residual[mod->fmap[grid->elem1d[ie].nodes[i]]] -= elem_rhs[i];
#endif
            }
        }
        //}
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#ifdef _PETSC
    // Assemble the PETSc residual vector
    ierr = VecAssemblyBegin(sm->residual);
    ierr = VecAssemblyEnd(sm->residual);
#else
    if (mod->amICoupled == NO) {
        
        // CHECK RESIDUAL FOR SMALL ENTRIES
        for (i = 0; i < grid->nnodes; i++) {
            if (fabs(sm->residual[i]) < SMALL) {sm->residual[i] = 0.0;}
        }
        
#ifdef _DEBUG
        if (DEBUG == ON
#ifdef _MESSG
            && grid->smpi->myid == 0
#endif
            ) {
            printScreen_resid("transport residual", sm->residual, grid->nnodes, mod->nsys, __LINE__, __FILE__);
            //exit(-1);
        }
#endif
    }
#endif
}

