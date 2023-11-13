
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_diffusive_load.c This file collections functions responsible for loading the DW matrix */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void perturb_dw_head(SMODEL *mod, int nnodes, int ie, int dim, double *elem_mat, int DEBUG);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Forms the SW2D diffusive wave Newton Jacobi matrix.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out]  mod a pointer to an AdH model where depth should be updated
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
 * where c is the conversion factor and n is mannings coefficient.*
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_diffusive_load(SSUPER_MODEL *sm, int imod) {
    
    SMODEL *mod = &(sm->submodel[imod]);
#ifdef _DEBUG
    assert(mod->flag.DIFFUSIVE_WAVE == ON);
    assert(mod->nsys == 1);
#endif
    
    int DEBUG = OFF;
    int DEBUG_MATRIX = OFF;
    int DEBUG_PICKETS = OFF;
#ifdef _DEBUG
    if (mod->t_prev > DEBUG_TIME) debug.load = ON;
    if (debug.load == ON) {
        DEBUG = ON;
    }
    if (debug.matrix == ON) {
        DEBUG_MATRIX = ON;
    }
#endif
    
    // aliases
    SSW_2D *sw2 = mod->sw->d2;
    SGRID *grid = mod->grid;
    
    // local variable declarations
    int i = 0, j = 0, ie = UNSET_INT;
    int nmatrix = NDONQUAD * NDONQUAD; // define using max 2D element size (quadrilateral)
    double elem_mat[nmatrix];
#ifdef _PETSC
    int ierr;
#endif
    
    // initialize matrix
    if (mod->amICoupled == NO) { // if this is a monolothic run, do this in the super load routine
#ifdef _PETSC
        ierr = MatZeroEntries(sm->A);
#else
        init_adh_matrix(grid->nnodes, mod->max_nsys_sq, sm->matrix, sm->diagonal);
#endif
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 2D ELEMENT CONTRIBUTIONS
    
    for (ie = 0; ie < grid->nelems2d; ie++) {
        if (mod->str_values[grid->elem2d[ie].string].phys_flag == OL_FLAG
#ifdef _ADH_GROUNDWATER
#ifdef _DWGW_COUPLING
                && ( grid->elem2d[ie].bflag == 0 /* Surface element of a 3D GW model*/
                    || grid->elem2d[ie].bflag == UNSET_INT /* Usual DW model*/)
#endif
#endif
                ) {
            
            sarray_init_dbl(elem_mat, nmatrix);
            perturb_dw_head(mod, grid->elem2d[ie].nnodes, ie, 2, elem_mat, DEBUG);
            
#ifdef _DEBUG
            if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
            if (DEBUG)  {
                printf("2d element: %d\n",ie);
                sarray_printScreen_dbl(elem_mat, nmatrix, "elem_mat");
            }
#endif
            // assembles the element contributions
            fe_global_matrix_assemble_scalar(grid, mod->nsys, grid->elem2d[ie].nnodes, grid->elem2d[ie].nodes, mod->fmap, elem_mat//, sm->diagonal, sm->matrix
#ifdef _PETSC
                    , sm->A, sm->Istart, sm->Iend, sm->ownership_range
#else
        , sm->diagonal, sm->matrix
#endif
                    );
        }
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 1D ELEMENT CONTRIBUTIONS
    
    for (ie = 0; ie < grid->nelems1d; ie++) {
        if (grid->elem1d[ie].string != UNSET_INT) {
            if (mod->str_values[grid->elem1d[ie].string].phys_flag == OL_FLAG) {
                
                sarray_init_dbl(elem_mat, nmatrix);
                perturb_dw_head(mod, grid->elem1d[ie].nnodes, ie, 1, elem_mat, DEBUG);
                
#ifdef _DEBUG
                if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
                if (DEBUG)  {
                    printf("1d element: %d\n",ie);
                    sarray_printScreen_dbl(elem_mat, nmatrix, "elem_mat");
                }
#endif
                // assembles the element contributions
                fe_global_matrix_assemble_scalar(grid, mod->nsys, grid->elem1d[ie].nnodes, grid->elem1d[ie].nodes, mod->fmap, elem_mat//, sm->diagonal, sm->matrix
#ifdef _PETSC
                    , sm->A, sm->Istart, sm->Iend, sm->ownership_range
#else
        , sm->diagonal, sm->matrix
#endif
                    );
                
            }
        }
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
#ifdef _PETSC
    // Assemble the PETSc Mat containing the Jacobian
    ierr = MatAssemblyBegin(sm->A, MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(sm->A, MAT_FINAL_ASSEMBLY);
#else
    // checks the diagonal for nonexistent entries
    if (mod->amICoupled == NO) { // if this is a monolothic run, do this in the super load routine
    //    check_matrix_diagonal_for_nonexistant_entries(grid->nnodes, mod->max_nsys_sq, mod->nsys, sm->diagonal);
    
#ifdef _MESSG
        comm_update_double(sm->diagonal, 1, mod->grid->smpi); // CJT :: might should be sm->supersmpi
#endif
    
#ifdef _DEBUG
        if (DEBUG_MATRIX == ON
#ifdef _MESSG
            && grid->smpi->myid == 0
#endif
            ) {
            printScreen_matrix("diffusive wave matrix", sm->diagonal, sm->matrix, grid->nnodes, mod->max_nsys_sq, __LINE__, __FILE__);
            exit(-1);
        }
#endif
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Adds perturbation to depth for each node in element for Newton Jacobian calculation.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] elem_mat stores the Jacobian, elemental matrix
 *  @param[in]  mod a pointer to an AdH model where depth should be updated
 *  @param[in] ie the element
 *  @param[in] dim the dimension of the resid
 *  @param[in] DEBUG a debug option
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void perturb_dw_head(SMODEL *mod, int nodes_on_element, int ie, int dim, double *elem_mat, int DEBUG) {
    
#ifdef _DEBUG
    assert(dim == 1 || dim == 2);
    assert(nodes_on_element == 2 || nodes_on_element == 3);
#endif
    
    int i,j;
    double epsilon = 0., epsilon2 = 0.;
    double perturbation = sqrt(SMALL);
    double elem_rhs_h_P[nodes_on_element];    // +head perturbation
    double elem_rhs_h_M[nodes_on_element];    // -head perturbation
    
    // get perturbation scale over element
    int GlobalNodeID = UNSET_INT;
    double max_h = 0.;
    for (i=0; i<nodes_on_element; i++) {
        if (dim == 2) {
            GlobalNodeID = mod->grid->elem2d[ie].nodes[i];
        } else {
            GlobalNodeID = mod->grid->elem1d[ie].nodes[i];
        }
        if (fabs(mod->sw->d2->head[GlobalNodeID]) > max_h) max_h = fabs(mod->sw->d2->head[GlobalNodeID]);
    }
    
    // *****************************************************************************
    // perturbation of h ***********************************************************
    NUM_DIFF_EPSILON_GENERAL(epsilon, epsilon2, max_h, perturbation);    // calculates epsilon and 2*epsilon
    
#ifdef _DEBUG
    if (DEBUG) printf("\ndiffusive wave perturbation size: %30.20f\n",epsilon);
#endif
    
    for (i=0; i<nodes_on_element; i++) {
        
        
        if (dim == 1) {
            
#ifdef _DEBUG
            if (DEBUG) { printf("\n\nepsilon: %20.10e \t perturbration: %20.10e nodes_on_elem: %d\n",epsilon,perturbation,nodes_on_element);
                printf("\nload :: pertubing +h for %dd element %d\n",dim,ie);
            }
#endif
            fe_diffusive_boundary_resid(mod,elem_rhs_h_P,ie,perturbation,i,PERTURB_H,+1,DEBUG);
            
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -h for %dd element %d\n",dim,ie);
#endif
            fe_diffusive_boundary_resid(mod,elem_rhs_h_M,ie,perturbation,i,PERTURB_H,-1,DEBUG);
            
        } else {
            
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +h for %dd element %d\n",dim,ie);
#endif
            fe_diffusive_body_resid(mod,elem_rhs_h_P,ie,perturbation,i,PERTURB_H,+1,DEBUG);
            
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -h for %dd element %d\n",dim,ie);
#endif
            fe_diffusive_body_resid(mod,elem_rhs_h_M,ie,perturbation,i,PERTURB_H,-1,DEBUG);
            
        }
        
        elem_matrix_deriv_1dof(i, nodes_on_element, elem_rhs_h_P, elem_rhs_h_M, elem_mat, epsilon2); // gradient
    }
}



