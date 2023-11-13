/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 * \brief     Calculates the 3D continuity shallow water residual.
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
 *    \weakSWContSurDDD{e}{i}{h} \\
 *    \weakSWContBedDDD{e}{i}{h}
 * \f}
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_wvel_resid(SSUPER_MODEL *sm, int imod) {
    
    SMODEL *mod = &(sm->submodel[imod]);
    SSW_3D *sw3 = mod->sw->d3;
    SGRID *grid = mod->grid;
    
#ifdef _DEBUG
    assert(mod->flag.SW3_FLOW == ON);
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
    
    int i, j, k, ie, ie2d, ie3d, icol, idof;
    ID_LIST_ITEM *ptr;
    int new_column_flag = YES;
    
    double elem_rhs[NDONPRISM];            // the full elemental residual for each equation
    double elem_rhs_noBed[NDONPRISM];      // the elemental residual w/o bed contributions
    sarray_init_dbl(elem_rhs, NDONPRISM);
    sarray_init_dbl(elem_rhs_noBed, NDONPRISM);

#ifdef _PETSC
    // Declare PETSc-specific variables
    PetscInt ierr;
    int ifmap;
    PetscScalar values[3];
#endif
    
    /* Store elemental residual through surface/bed into a nodal array for monitoring */
    sarray_init_dbl(sw3->vertical_node_flux, grid->nnodes);
    
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
        sarray_init_dbl(sw3->elem_rhs_supg_cont[i], grid->nelems3d);
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 2D ELEMENT CONTRIBUTIONS
    
    for (ie2d=0; ie2d<grid->nelems2d; ie2d++) {
        
        // intialize and copy element variables
        fe_wvel_boundary_resid(mod, elem_rhs, elem_rhs_noBed, ie2d, 0, UNSET_INT, PERTURB_NONE, 0, DEBUG);
        
        /* sets the Dirichlet boundary conditions */
        for(i = 0; i < grid->elem2d[ie2d].nnodes; i++) {
            if(mod->bc_mask[grid->elem2d[ie2d].nodes[i]]==YES) {
                elem_rhs[i]=0.;
                elem_rhs_noBed[i]=0.;
            }
        }
        
#ifdef _DEBUG
        if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
        if (DEBUG) {
            printf("2d element: %d \n",ie2d);
            dof1_printScreen_array("elem_rhs",elem_rhs, grid->elem2d[ie2d].nnodes,__LINE__,__FILE__);
        }
#endif
        // Avengers Assemble
        for (i = 0; i < grid->elem2d[ie2d].nnodes; i++) {
#ifdef _PETSC
            j = grid->elem2d[ie2d].nodes[i];
            sw3->vertical_node_flux[j] -= elem_rhs_noBed[i];
            if(grid->node[j].resident_pe == grid->smpi->myid){
                //sw3->vertical_node_flux[j] -= elem_rhs_noBed[i];
                ifmap = mod->fmap_wvel[j];
                values[0] = -elem_rhs[i];
                values[1] = 0;
                values[2] = 0;
                ierr = VecSetValuesBlockedLocal(sm->residual, 1, &(ifmap), &(values),
                            ADD_VALUES);
            }
#else
            j = grid->elem2d[ie2d].nodes[i];
            sw3->vertical_node_flux[j] -= elem_rhs_noBed[i];
            sm->residual[mod->fmap_wvel[j]] -= elem_rhs[i];
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
            fe_wvel_body_resid(mod, elem_rhs, ie3d, 0., UNSET_INT, PERTURB_NONE, 0, DEBUG);
            
            /* sets the Dirichlet boundary conditions */
            for(i = 0; i < grid->elem3d[ie3d].nnodes; i++) {
                if(mod->bc_mask[grid->elem3d[ie3d].nodes[i]] == YES) {
                    elem_rhs[i]=0.;
                }
            }
            
#ifdef _DEBUG
            if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
            if (DEBUG) {
                printf("3d element: %d \n",ie3d);
                dof1_printScreen_array("elem_rhs",elem_rhs, grid->elem3d[ie3d].nnodes,__LINE__,__FILE__);
            }
#endif
            // Avengers Assemble
            for (i = 0; i < grid->elem3d[ie3d].nnodes; i++) {
#ifdef _PETSC
                j = grid->elem3d[ie3d].nodes[i];
                sw3->vertical_node_flux[j] -= elem_rhs[i];
                if(grid->node[j].resident_pe == grid->smpi->myid){
                    //sw3->vertical_node_flux[j] -= elem_rhs[i];
                    ifmap = mod->fmap_wvel[j];
                    values[0] = -elem_rhs[i];
                    values[1] = 0;
                    values[2] = 0;
                    ierr = VecSetValuesBlockedLocal(sm->residual, 1, &(ifmap), &(values),
                                ADD_VALUES);
                }
#else
                j = grid->elem3d[ie3d].nodes[i];
                sw3->vertical_node_flux[j] -= elem_rhs[i];
                sm->residual[mod->fmap_wvel[j]] -= elem_rhs[i];
#endif
            }
            new_column_flag = NO;
            ptr = ptr->next;
        }
    }
    
    //if (mod->t_prev > 1400) exit(-1);

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
            printf("printing residual: fe_wvel_resid @ line %d:%s \n",__LINE__, __FILE__);
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
#endif
}
