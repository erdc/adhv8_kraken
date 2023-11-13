/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_bedload_load.c This file collections functions responsible for loading the SW 2D bed load matrix
 /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void perturb_bl_c(SMODEL *mod, int nnodes, int ie, int dim, double *elem_mat, int DEBUG);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Forms the 2D bed load Newton Jacobi matrix.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out]  mod a pointer to an AdH model where depth should be updated
 *
 * \details   Solves the following weak, discrete diffusive wave equation: \n
 * \f{eqnarray*}{ \weakBedTransport \f}
 * \n where \f$ c^{\,j,\,h}_{bl} \f$ is the bed concentration of constituent j and \f$ \delta^{\,h\,} \f$ is the bed load thickness.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_bedload_load(SMODEL *mod) {
    
#ifdef _DEBUG
    assert(mod->flag.SEDIMENT == ON);
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
    
    // initialize matrix
    if (mod->amICoupled == NO) {
        init_adh_matrix(grid->nnodes, mod->max_nsys_sq, mod->matrix, mod->diagonal);
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 2D ELEMENT CONTRIBUTIONS
    
    for (ie = 0; ie < grid->nelems2d; ie++) {
        if (grid->elem2d[ie].string != UNSET_INT) {
            
            if (mod->str_values[grid->elem2d[ie].string].phys_flag == SW3_FLAG && mod->str_values[grid->elem2d[ie].string].flow.bc_flag != BCT_BED) continue;
            
            sarray_init_dbl(elem_mat, nmatrix);
            perturb_bl_c(mod, grid->elem2d[ie].nnodes, ie, 2, elem_mat, DEBUG);
            
#ifdef _DEBUG
            if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
            if (DEBUG) {
                printf("2d element: %d\n",ie);
                sarray_dbl_printScreen(elem_mat, nmatrix, "elem_mat");
            }
#endif
            // assembles the element contributions
            fe_global_matrix_assemble_scalar(grid, mod->nsys, grid->elem2d[ie].nnodes, grid->elem2d[ie].nodes, mod->fmap, elem_mat, mod->diagonal, mod->matrix);
        }
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 1D ELEMENT CONTRIBUTIONS
    
    for (ie = 0; ie < grid->nelems1d; ie++) {
        if (grid->elem1d[ie].string != UNSET_INT) {
            
            sarray_init_dbl(elem_mat, nmatrix);
            perturb_bl_c(mod, grid->elem1d[ie].nnodes, ie, 1, elem_mat, DEBUG);
            
#ifdef _DEBUG
            if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
            if (DEBUG) {
                printf("1d element: %d\n",ie);
                sarray_dbl_printScreen(elem_mat, nmatrix, "elem_mat");
            }
#endif
            // assembles the element contributions
            fe_global_matrix_assemble_scalar(grid, mod->nsys, grid->elem1d[ie].nnodes, grid->elem1d[ie].nodes, mod->fmap, elem_mat, mod->diagonal, mod->matrix);
            
        }
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    if (mod->amICoupled == NO) {
        
        // checks the diagonal for nonexistent entries
        check_matrix_diagonal_for_nonexistant_entries(grid->nnodes, mod->max_nsys_sq, mod->nsys, mod->diagonal);
#ifdef _MESSG
        comm_update_double(mod->diagonal, 1, mod->grid->smpi);
#endif
        
#ifdef _DEBUG
        if (DEBUG_matrix == ON
#ifdef _MESSG
            && grid->smpi->myid == 0
#endif
            ) {
            printScreen_matrix("2d bedload matrix", mod->diagonal, mod->matrix, grid->nnodes, mod->max_nsys_sq, __LINE__, __FILE__);
            exit(-1);
        }
#endif
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Adds perturbation to bed load concentration for each node in element for Newton Jacobian calculation.
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
void perturb_bl_c(SMODEL *mod, int nodes_on_element, int ie, int dim, double *elem_mat, int DEBUG) {
    
#ifdef _DEBUG
    assert(dim == 1 || dim == 2);
    assert(nodes_on_element == 2 || nodes_on_element == 3);
#endif
    
    int i;
    double epsilon = 0., epsilon2 = 0.;
    double perturbation = 0.1;
    double elem_rhs_P[nodes_on_element];    // +c perturbation
    double elem_rhs_M[nodes_on_element];    // -c perturbation
    
    int ised = mod->ised;
    double c_eps = 0.;
    
    for (i=0; i<nodes_on_element; i++) {
        
        if (dim == 2) {
            GlobalNodeID = mod->grid->elem2d[ie].nodes[i];
        } else {
            GlobalNodeID = mod->grid->elem1d[ie].nodes[i];
        }
        
        c_eps = mod->sed->bedload[ised].c[GlobalNodeID];
        if (c_eps < 1.0) {
            c_eps = 1.0;
        }
        if (c_eps > mod->sed->grain[ised].reference_c) {
            c_eps = mod->sed->grain[ised].reference_c;
        }
        NUM_DIFF_EPSILON_GENERAL(epsilon, epsilon2, c_eps, perturbation);
        
        if (dim == 1) {
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +c for %dd element\n",dim,ie);
#endif
            fe_bedload_1d_elem_resid(mod,elem_rhs_P,ie,pertubation,i,PERTURB_C,+1,DEBUG);
            
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -c for %dd element\n",dim,ie);
#endif
            fe_bedload_1d_elem_resid(mod,elem_rhs_M,ie,pertubation,i,PERTURB_C,-1,DEBUG);
            
        } else {
            
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +c for %dd element\n",dim,ie);
#endif
            fe_bedload_2d_elem_resid(mod,elem_rhs_P,ie,pertubation,i,PERTURB_C,+1,DEBUG);
            
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -c for %dd element\n",dim,ie);
#endif
            fe_bedload_2d_elem_resid(mod,elem_rhs_M,ie,pertubation,i,PERTURB_C,-1,DEBUG);
            
        } else {
            
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> dimension (dim) must be either 1 or 2.");
            
        }
        
        elem_matrix_deriv_1dof(i, nodes_on_element, elem_rhs_h_P, elem_rhs_h_M, elem_mat, epsilon2); // gradient
    }
}
