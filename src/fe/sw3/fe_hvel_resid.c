/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Calculates the 3D depth-averaged and momentum shallow water residuals.
 * \author    Charlie Berger, Ph.D.
 * \author    Gaurav Savant, Ph.D.
 * \author    Corey Trahan, Ph.D.
 * \bug       none
 * \warning   none
 * \copyright AdH
 *
 * @param[in,out] mod (SMODEL *) a pointer to an AdH model struct
 *
 * \details   Solves the following weak, discrete equation: \n
 * \f{eqnarray*}{
 *   \weakSWDaContReducedKinematic{i} \\
 *   \weakSWMxDDD{e}{i}{h} \\
 *   \weakSWMxDDD{e}{i}{h}
 * \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_hvel_resid(SSUPER_MODEL *sm, int imod) {
    
    SMODEL *mod = &(sm->submodel[imod]);
#ifdef _DEBUG
    assert(mod->flag.SW3_FLOW == ON);
    assert(mod->nsys == 3);
#endif
    
    int DEBUG = OFF;
    int DEBUG_PICKETS = OFF;
#ifdef _DEBUG
    if (mod->t_prev > DEBUG_TIME) debug.residual = ON;
    if (debug.residual == ON) {
        DEBUG = ON;
    }
    time_t time1;  time(&time1);
#endif
    
    int i, j, k, icol, ie2d, ie3d, iseg, surface_nodeID;
    int new_column_flag = YES;
    ID_LIST_ITEM *ptr;

#ifdef _PETSC
    // Declare PETSc-specific variables
    PetscInt ierr;
    PetscScalar values[3];
    int GlobalNodeID;
#endif
    
    // alias
    SSW_3D *sw3 = mod->sw->d3;
    SGRID *grid = mod->grid;
    
    DOF_3 elem_rhs[NDONPRISM];
    
    /* initialize the residual */
    if (mod->amICoupled == NO) {
#ifdef _PETSC
        ierr = VecZeroEntries(sm->residual);
#else
        sarray_init_dbl(sm->residual, grid->nnodes * mod->max_nsys);
#endif
    }
    
    /* initialize the supg storage arrays */
    for (i=0; i<MAX_NNODES_ON_ELEM3D; i++) {
        sarray_init_dbl(sw3->elem_rhs_supg_dacont[i], grid->nelems3d);
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 2D ELEMENT CONTRIBUTIONS
    
    for (ie2d=0; ie2d<grid->nelems2d; ie2d++) {
        
        // intialize and copy element variables
        dof3_init_array(elem_rhs, NDONPRISM);
        fe_hvel_boundary_resid(mod, elem_rhs, ie2d, 0, UNSET_INT, PERTURB_NONE, 0, NULL, DEBUG);
        
        /* sets the Dirichlet boundary conditions */
        for (i=0; i<grid->elem2d[ie2d].nnodes; i++) {
            j = mod->nsys * grid->elem2d[ie2d].nodes[i];
            if (mod->bc_mask[j] == YES) {
                elem_rhs[i].x_eq = 0.;
            }
            if (mod->bc_mask[j + 1] == YES) {
                elem_rhs[i].y_eq = 0.;
            }
            if (mod->str_values[grid->elem2d[ie2d].string].flow.bc_flag == BCT_FRS) {
                if (mod->bc_mask[j + 2] == YES) {
                    elem_rhs[i].c_eq = 0.;
                }
            }
        }
        
#ifdef _DEBUG
        if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
        if (DEBUG) {
            printf("2d element: %d \n",ie2d);
            dof3_printScreen_array("elem_rhs",elem_rhs, grid->elem2d[ie2d].nnodes,__LINE__,__FILE__);
        }
#endif
        // Avengers Assemble
        for (i = 0; i < grid->elem2d[ie2d].nnodes; i++) {
#ifdef _PETSC
            GlobalNodeID = grid->elem2d[ie2d].nodes[i];
            if(grid->node[GlobalNodeID].resident_pe == grid->smpi->myid){
                j = mod->fmap[GlobalNodeID];
                values[0] = -elem_rhs[i].x_eq;
                values[1] = -elem_rhs[i].y_eq;
                values[2] = 0;
                ierr = VecSetValuesBlockedLocal(sm->residual, 1, &j, values,
                        ADD_VALUES);

                // Set continuity eq for surface node separately
                iseg = find_vertical_segment(grid, grid->elem2d[ie2d].nodes[i],
                        grid->vertical_hash);
                surface_nodeID = grid->vertical_list[iseg]->id;
                j = mod->fmap[surface_nodeID] * mod->nsys + 2;
                values[0] = -elem_rhs[i].c_eq;
                ierr = VecSetValuesLocal(sm->residual, 1, &j, &(values[0]),
                        ADD_VALUES);
            }
#else
            j = mod->nsys * mod->fmap[grid->elem2d[ie2d].nodes[i]];
            sm->residual[j] -= elem_rhs[i].x_eq;
            sm->residual[j + 1] -= elem_rhs[i].y_eq;
            
            // The depth integrated continuity equation makes contributions to the surface node only
            iseg = find_vertical_segment(grid, grid->elem2d[ie2d].nodes[i], grid->vertical_hash);
            surface_nodeID = grid->vertical_list[iseg]->id;
            sm->residual[mod->fmap[surface_nodeID] * mod->nsys + 2] -= elem_rhs[i].c_eq;
#endif
        }
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 3D ELEMENT CONTRIBUTIONS
    
    for(icol=0; icol<grid->ncolumns; icol++)  {
        new_column_flag = YES;
        
        ptr = grid->column_list[icol];
        while(ptr->next != NULL) {
            ie3d = ptr->id;
            
            dof3_init_array(elem_rhs, NDONPRISM);
            fe_hvel_body_resid(mod, elem_rhs, ie3d, 0., UNSET_INT, PERTURB_NONE, 0, NULL, DEBUG);
            
            /* sets the Dirichlet boundary conditions */
            for (i = 0; i < grid->elem3d[ie3d].nnodes; i++) {
                j = mod->nsys * grid->elem3d[ie3d].nodes[i];
                if (mod->bc_mask[j] == YES) {
                    elem_rhs[i].x_eq = 0.;
                }
                if (mod->bc_mask[j] == 2) {
                    elem_rhs[i].y_eq = elem_rhs[i].x_eq * sw3->tanvec[grid->elem3d[ie3d].nodes[i]].x + elem_rhs[i].y_eq * sw3->tanvec[grid->elem3d[ie3d].nodes[i]].y;
                    elem_rhs[i].x_eq = 0.;
                }
                if (mod->bc_mask[j + 1] == YES) {
                    elem_rhs[i].y_eq = 0.;
                }
            }
            
#ifdef _DEBUG
            if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
            if (DEBUG) {
                printf("3d element: %d \n",ie3d);
                dof3_printScreen_array("elem_rhs",elem_rhs, grid->elem3d[ie3d].nnodes,__LINE__,__FILE__);
            }
#endif
            // Avengers Assemble
            for (i = 0; i < grid->elem3d[ie3d].nnodes; i++) {
#ifdef _PETSC
                GlobalNodeID = grid->elem3d[ie3d].nodes[i];
                if(grid->node[GlobalNodeID].resident_pe == grid->smpi->myid){
                    values[0] = -elem_rhs[i].x_eq;
                    values[1] = -elem_rhs[i].y_eq;
                    values[2] = 0;
                    j = mod->fmap[GlobalNodeID];
                    ierr = VecSetValuesBlockedLocal(sm->residual, 1, &j, values,
                            ADD_VALUES);

                    // Set continuity eq for surface node separately
                    iseg = find_vertical_segment(grid, grid->elem3d[ie3d].nodes[i],
                            grid->vertical_hash);
                    surface_nodeID = grid->vertical_list[iseg]->id;
                    j = mod->fmap[surface_nodeID] * mod->nsys + 2;
                    values[0] = -elem_rhs[i].c_eq;
                    ierr = VecSetValuesLocal(sm->residual, 1, &j, &(values[0]),
                        ADD_VALUES);
                }
#else
                j = mod->nsys * mod->fmap[grid->elem3d[ie3d].nodes[i]];
                sm->residual[j]     -= elem_rhs[i].x_eq;
                sm->residual[j + 1] -= elem_rhs[i].y_eq;
                
                // The depth integrated continuity equation makes contributions to the surface node only
                iseg = find_vertical_segment(grid, grid->elem3d[ie3d].nodes[i], grid->vertical_hash);
                surface_nodeID = grid->vertical_list[iseg]->id;
                sm->residual[mod->fmap[surface_nodeID] * mod->nsys + 2] -= elem_rhs[i].c_eq;
#endif
            }
            new_column_flag = NO;
            ptr = ptr->next;
        }
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
#ifdef _PETSC
    ierr = VecAssemblyBegin(sm->residual);
    ierr = VecAssemblyEnd(sm->residual);
#else
    if (mod->amICoupled == NO) {
#ifdef _DEBUG
        if (DEBUG == ON
#ifdef _MESSG
            && grid->smpi->myid == 0
#endif
            ) {
            //printScreen_resid("hvel residual", sm->residual, grid->nnodes, nsys, __LINE__, __FILE__);
            printf("\n");
            printf("printing residual: fe_hvel_resid @ line %d:%s \n",__LINE__, __FILE__);
            for (i = 0; i < mod->grid->nnodes; i++) {
                j = i * mod->nsys;
                printf("i=%d \t {%20.10f,%20.10f,%20.10f} \t final residual = ", (i + 1),mod->grid->node[i].x,mod->grid->node[i].y,mod->grid->node[i].z);
                for (k = 0; k < mod->nsys; k++)
                    printf(" %30.20e ", sm->residual[j + k]);
                printf("\n");
            }
        }
#endif
    }
    
#ifdef _DEBUG
    if (DEBUG == ON) {
        time_t time2;  time(&time2);
        TIME_IN_HVEL_RESID += difftime(time2,time1);
    }
#endif
    
#endif
}
