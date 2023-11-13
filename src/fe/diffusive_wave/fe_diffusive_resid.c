/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Calculates the diffusive wave residual for overland flow.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out]  model a pointer to an AdH model where Newton residual should be updated
 *
 * \details   Solves the following weak, discrete diffusive wave equation: \n
 * \f{eqnarray*}{
 *  \resid{i}{dw}{} =\bodyTime{\,2d}{e}{\phidd{i}}{\depth{h}} \,-
 *           \bodyConv{\,2d}{e}{\phidd{i}}{ (\depth{h} \, k \, \grad \elev{h}) } \,-
 *           \bodySource{\,2d}{e}{\phidd{i}}{\srcGen{dw}{h}{}} \,+
 *           \bcConv{}{e}{\phidd{i}}{(\depth{h} \, k \, \grad \elev{h})} \,+
 *           \bodySUPG{\supg{i}{dw}{e}}
 * \f}
 * which can also be written as \n
 * \f{eqnarray*}{
 *  \resid{i}{dw}{} =\bodyTime{\,2d}{e}{\phidd{i}}{\depth{h}} \,-
 *           \bodyConv{\,2d}{e}{\phidd{i}}{ (\depth{h} \, \velb{h}} \,-
 *           \bodySource{\,2d}{e}{\phidd{i}}{\srcGen{dw}{h}{}} \,+
 *           \bcConv{}{e}{\phidd{i}}{(\velb{h} \,\depth{h})} \,+
 *           \bodySUPG{\supg{i}{dw}{e}}
 * \f}
 * where \n
 * \f{eqnarray*}{
 * u &= k \, \deriv{\elev{h}}{x} \\
 * v &= k \, \deriv{\elev{h}}{y} \\
 * k &= \frac{c^2 H^{4/3}}{n^2 \, \|\velh\|}
 * \f}
 * where c is the conversion factor and n is mannings coefficient.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_diffusive_resid(SSUPER_MODEL *sm, int imod) {
    
    SMODEL *mod = &(sm->submodel[imod]);
#ifdef _DEBUG
    // should only call this return for stand-alone diffusive wave model
    assert(mod->flag.DIFFUSIVE_WAVE == ON);
    assert(mod->nsys == 1);
#endif
    
    int DEBUG = OFF;
    int DEBUG_PICKETS = OFF;
    
    int ie, i;
    double elem_rhs[NDONTRI];     /* the mod->residual for each equation */
   
#ifdef _PETSC
    PetscInt row_index, ierr;
    PetscScalar values[3];
    values[1] = 0; values[2] = 0;
    int GnodeID, ifmap;
#endif 
    
    // aliases
    SSW_2D *sw_diff = mod->sw->d2;
    SGRID *grid = mod->grid;
    
    int myid = 0;
#ifdef _MESSG
    myid = grid->smpi->myid;
#endif
    
    /* initializes the residual */
    if (mod->amICoupled == NO) {
#ifdef _PETSC
        // Zero residual
        ierr = VecZeroEntries(sm->residual);
#else
        sarray_init_dbl(sm->residual, grid->nnodes * mod->nsys);
#endif
    }
    
#ifdef _DEBUG
    if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__, __LINE__);
#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /* loops over the valid elements and loads the element contributions */
    for (ie = 0; ie < grid->nelems2d; ie++) {
        if (mod->str_values[grid->elem2d[ie].string].phys_flag == OL_FLAG
#ifdef _ADH_GROUNDWATER
#ifdef _DWGW_COUPLING
                && ( grid->elem2d[ie].bflag == 0 /* Surface element of a 3D GW model*/
                    || grid->elem2d[ie].bflag == UNSET_INT /* Usual DW model*/)
#endif
#endif
                ) {
            
            // forms the 2D element residual
            fe_diffusive_body_resid(mod, elem_rhs, ie, 0, UNSET_INT, UNSET_INT, 1, DEBUG);
#ifdef _DEBUG
            if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__, __LINE__);
#endif
            
            // assembles the elemental residual contributions in global
            for (i = 0; i < grid->elem2d[ie].nnodes; i++) {
#ifdef _PETSC
                values[0] = -elem_rhs[i];
                GnodeID = grid->elem2d[ie].nodes[i];
                if(grid->node[GnodeID].resident_pe == myid){
                    VecSetValuesBlockedLocal(sm->residual,1,&(GnodeID),values,ADD_VALUES);
                }
#else
                sm->residual[grid->elem2d[ie].nodes[i]] -= elem_rhs[i];
#endif
                
            }
        }
    }
    
    /* loops over the valid elements and loads the element contributions*/
    for (ie = 0; ie < grid->nelems1d; ie++) {
        if (mod->str_values[ grid->elem1d[ie].string ].phys_flag == OL_FLAG) {

            fe_diffusive_boundary_resid(mod, elem_rhs, ie, 0, UNSET_INT, UNSET_INT, 1, DEBUG);
#ifdef _DEBUG
            if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__, __LINE__);
#endif
            /* assembles the element contributions */
            for (i = 0; i < grid->elem1d[ie].nnodes; i++) {
#ifdef _PETSC
                values[0] = -elem_rhs[i];
                GnodeID = grid->elem1d[ie].nodes[i];
                if(grid->node[GnodeID].resident_pe == myid){
                    VecSetValuesBlockedLocal(sm->residual,1,&(GnodeID),values,ADD_VALUES);
                }
#else
                sm->residual[grid->elem1d[ie].nodes[i]] -= elem_rhs[i];
#endif
                
            }
        }
    }
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
#ifdef _PETSC
    ierr = VecAssemblyBegin(sm->residual);
    ierr = VecAssemblyEnd(sm->residual);

    // TODO: I don't think we need to check for small
    // entries with PETSc. Just in case, I've sketched out
    // code to check. This code has NOT been tested.
    // There may be a better/cleaner/more efficient way
    // to impelement small value checking for the residual
    // than what I've sketched here. - SCE
    //
    ///* checks the residual for very small entries */
    //PetscScalar const *read_values;

    //if (mod->amICoupled == NO) {
    //    ierr = VecGetArrayRead(sm->residual,&(read_values));
    //    values[0] = 0; values[1] = 0; values[2] = 0;
    //    for (i = 0; i < grid->nnodes; i++) {
    //        if (fabs(read_values[i*sm->max_nsys]) < SMALL) {
    //            ierr = VecSetValuesBlockedLocal(sm->residual,1,&i,values,INSERT_VALUES);
    //        }
    //    }
    //    ierr = VecRestoreArrayRead(sm->residual,&(read_values));

    //    // Assemble new residual without small values
    //    ierr = VecAssemblyBegin(sm->residual);
    //    ierr = VecAssemblyEnd(sm->residual);
    //}
#else
    /* checks the residual for very small entries */
    if (mod->amICoupled == NO) {
        for (i = 0; i < grid->nnodes; i++) {
            if (fabs(sm->residual[i]) < SMALL) {
                sm->residual[i] = 0.0;
            }
        }
    }
 
#ifdef _DEBUG
    if (DEBUG_PICKETS == ON) tl_check_all_pickets(__FILE__, __LINE__);
#endif
    

#endif
}

