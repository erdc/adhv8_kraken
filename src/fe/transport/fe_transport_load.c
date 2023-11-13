/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_transport_load.c This file collections functions for loading the 2D and 3D transport matrix */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void perturb_c(SMODEL *mod, int nnodes, int ie, int dim, double *elem_mat, int DEBUG);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Forms the 2D or 3D transport Newton Jacobi matrix.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *   @param[in,out]  mod a pointer to an AdH model where depth should be updated
 *
 *  \details   Solves the following weak, discrete 2D or 3D transport equation: \n
 *  \f{eqnarray*}{
 *    \weakTrnsDA{2d}{e}{h}{\phidd{i}}{\cjhDA \, \depth{h}}{i}{t}{j} \\
 *    \weakTrns{3d}{e}{h}{\phiddd{i}}{\cjh}{i}{t}{j} \\
 *  \f}
 *  \n where \f$\cjh\f$ is constituent j's concentration, \f$\cjhDA\f$ is depth-averaged constituent j's concentration \n
 *  \f$ \depth{h} \f$ is the water depth and \f$ \velc{h}{t} \f$ and \f$ \velcDA{h}{t} \f$ are water velocities and  depth-averaged velocities.
 *
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_transport_load(SSUPER_MODEL *sm, int imod) {

    SMODEL *mod = &(sm->submodel[imod]);
#ifdef _DEBUG
    assert(mod->flag.TRANSPORT == ON);
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
    SGRID *grid = mod->grid;
    
    // local variable declarations
    int i = 0, j = 0, ie = UNSET_INT;
    int nmatrix = NDONPRISM * NDONPRISM; // define using max AdH element size (Triangular Prism)
    double elem_mat[nmatrix];
#ifdef _PETSC
    PetscInt ierr;
#endif
    
    
    // initialize matrix
    if (mod->amICoupled == NO) {
#ifdef _PETSC
        // Zero the PETSc matrix
        ierr = MatZeroEntries(sm->A);   // TODO: DO WE NEED THIS - SAM
#else
        init_adh_matrix(grid->nnodes, mod->max_nsys_sq, sm->matrix, sm->diagonal);
#endif
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 3D ELEMENT CONTRIBUTIONS

    for (ie = 0; ie < grid->nelems3d; ie++) {
        
        sarray_init_dbl(elem_mat, nmatrix);
        perturb_c(mod, grid->elem3d[ie].nnodes, ie, 3, elem_mat, DEBUG);
        
#ifdef _DEBUG
        if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
        if (DEBUG) {
            printf("3d element: %d\n",ie);
            sarray_printScreen_dbl(elem_mat, nmatrix, "elem_mat");
        }
#endif
        // assembles the element contributions
        fe_global_matrix_assemble_scalar(grid, mod->nsys, grid->elem3d[ie].nnodes, grid->elem3d[ie].nodes, mod->fmap, elem_mat//, sm->diagonal, sm->matrix
#ifdef _PETSC
                , sm->A, sm->Istart, sm->Iend, sm->ownership_range
#else
        , sm->diagonal, sm->matrix
#endif
                );
    }
    //printScreen_matrix("transport matrix after 3d additions", sm->diagonal, sm->matrix, grid->nnodes, mod->max_nsys_sq, __LINE__, __FILE__);
    //exit(-1);
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 2D ELEMENT CONTRIBUTIONS

    for (ie = 0; ie < grid->nelems2d; ie++) {
        //if (grid->wd_flag[ie] == 0 || grid->wd_flag[ie] == 1 || grid->wd_flag[ie] == 2) {
            if (grid->elem2d[ie].string != UNSET_INT) {
                
                sarray_init_dbl(elem_mat, nmatrix);
                perturb_c(mod, grid->elem2d[ie].nnodes, ie, 2, elem_mat, DEBUG);
                
#ifdef _DEBUG
                if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
                if (DEBUG) {
                    printf("2d element: %d\n",ie);
                    sarray_printScreen_dbl(elem_mat, grid->elem2d[ie].nnodes * grid->elem2d[ie].nnodes, "elem_mat");
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
        //}
    }
    //printScreen_matrix("transport matrix after 2d additions", sm->diagonal, sm->matrix, grid->nnodes, mod->max_nsys_sq, __LINE__, __FILE__);
    //exit(-1);
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 1D ELEMENT CONTRIBUTIONS
    
    for (ie = 0; ie < grid->nelems1d; ie++) {
        if (grid->elem1d[ie].string != UNSET_INT) {
            //if (mod->str_values[grid->elem1d[ie].string].phys_flag == SW2_FLAG) {
            
            sarray_init_dbl(elem_mat, nmatrix);
            perturb_c(mod, grid->elem1d[ie].nnodes, ie, 1, elem_mat, DEBUG);
            
#ifdef _DEBUG
            if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
            if (DEBUG) {
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
            //}
        }
    }
    
    //printScreen_matrix("transport matrix after 1d additions", sm->diagonal, sm->matrix, grid->nnodes, mod->max_nsys_sq, __LINE__, __FILE__);
    //exit(-1);

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
#ifdef _PETSC
    // Assemble the PETSc Mat containing the Jacobian
    ierr = MatAssemblyBegin(sm->A, MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(sm->A, MAT_FINAL_ASSEMBLY);
#else
    // checks the diagonal for nonexistent entries
    if (mod->amICoupled == NO) {
        check_matrix_diagonal_for_nonexistant_entries(grid->nnodes, mod->max_nsys_sq, sm->nsys, sm->diagonal);
        
#ifdef _MESSG
        comm_update_double(sm->diagonal, 1, mod->grid->smpi);
#endif
        
#ifdef _DEBUG
        if (DEBUG_MATRIX == ON
#ifdef _MESSG
            && grid->smpi->myid == 0
#endif
            ) {
            printScreen_matrix("transport matrix", sm->diagonal, sm->matrix, grid->nnodes, mod->max_nsys_sq, __LINE__, __FILE__);
            exit(-1);
        }
#endif
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Adds perturbation to transport concentration for each node in element for Newton Jacobian calculation.
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
 *  \note CJT \:: assumes there are no 1D element in a 3D run for now
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void perturb_c(SMODEL *mod, int nodes_on_element, int ie, int dim, double *elem_mat, int DEBUG) {
    
#ifdef _DEBUG
    assert(dim == 1 || dim == 2 || dim == 3);
    assert(nodes_on_element > 0 && nodes_on_element < NDONPRISM);
    assert(ie > -1);
#endif
    
    int i, index, GlobalNodeID = UNSET_INT;
    double perturbation = 0.1;
    double epsilon = 0., epsilon2 = 0., c_eps = 0.;
    double elem_rhs_P[nodes_on_element];    // +c perturbation
    double elem_rhs_M[nodes_on_element];    // -c perturbation
    
    sarray_init_dbl(elem_rhs_P,nodes_on_element);
    sarray_init_dbl(elem_rhs_M,nodes_on_element);
    
    double reference_conc = 0;
    double *c = NULL;
    if (mod->is_sediment_running) {
#ifdef _SEDIMENT
        index = mod->ised;
        c = mod->sed->susload[index].c;
        reference_conc = mod->sed->grain[index].reference_c;
#endif
    } else {
        index = mod->itrns;
        c = mod->con[index].concentration;
        reference_conc = mod->con[index].property[0];
    }
    
    double max_c_value = 0.;
    if (dim == 3) {
        for (i = 0; i<nodes_on_element; i++) {
            if (fabs(c[i]) > max_c_value) {
                max_c_value = fabs(c[i]);
            }
        }
    }
    
    for (i=0; i<nodes_on_element; i++) {
        
        /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        /* calculate perturbation  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        
        if (dim == 3) {
            GlobalNodeID = mod->grid->elem3d[ie].nodes[i];
        } else if (dim == 2) {
            GlobalNodeID = mod->grid->elem2d[ie].nodes[i];
        } else {
            GlobalNodeID = mod->grid->elem1d[ie].nodes[i];
        }
        
        if (dim == 1 || dim == 2) {
            c_eps = c[GlobalNodeID];
            if (c_eps < 1.0) {
                c_eps = 1.0;
            }
            if (c_eps > reference_conc) {
                c_eps = reference_conc;
            }
        } else {
            c_eps = max_c_value; // in 3d
        }
        
        
        NUM_DIFF_EPSILON_GENERAL(epsilon, epsilon2, c_eps, perturbation);
        //epsilon = 100;
        //epsilon2 = 2 * epsilon;
        
#ifdef _DEBUG
        if (DEBUG) printf("\nc perturbation size: %30.20f\n",epsilon);
#endif
        
        /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        /* calculate perturbed residuals ++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        
        if (dim == 1) {
            
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +c for %dd element %d\n",dim,ie);
#endif
            fe_2d_transport_boundary_resid(mod,elem_rhs_P,ie,epsilon,i,PERTURB_C,+1,DEBUG);
            
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -c for %dd element %d\n",dim,ie);
#endif
            fe_2d_transport_boundary_resid(mod,elem_rhs_M,ie,epsilon,i,PERTURB_C,-1,DEBUG);
            
        } else if (dim == 2) {
            
            
            if (mod->flag.SW2_TRANSPORT == ON || mod->flag.NS2_TRANSPORT) { // 2D TRANSPORT BOUNDARIES
#ifdef _DEBUG
                if (DEBUG) printf("\nload :: pertubing +c for %dd element %d\n",dim,ie);
#endif
                fe_2d_transport_body_resid(mod,elem_rhs_P,ie,epsilon,i,PERTURB_C,+1,DEBUG);
#ifdef _DEBUG
                if (DEBUG) printf("\nload :: pertubing -c for %dd element %d\n",dim,ie);
#endif
                fe_2d_transport_body_resid(mod,elem_rhs_M,ie,epsilon,i,PERTURB_C,-1,DEBUG);
                
            } else { // 3D TRANSPORT BOUNDARIES
#ifdef _DEBUG
                if (DEBUG) printf("\nload :: pertubing +c for %dd element %d\n",dim,ie);
#endif
                fe_3d_transport_boundary_resid(mod,elem_rhs_P,ie,epsilon,i,PERTURB_C,+1,DEBUG);
#ifdef _DEBUG
                if (DEBUG) printf("\nload :: pertubing -c for %dd element %d\n",dim,ie);
#endif
                fe_3d_transport_boundary_resid(mod,elem_rhs_M,ie,epsilon,i,PERTURB_C,-1,DEBUG);
            }
            
        } else {
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +c for %dd element %d\n",dim,ie);
#endif
            fe_3d_transport_body_resid(mod,elem_rhs_P,ie,epsilon,i,PERTURB_C,+1,DEBUG);
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -c for %dd element %d\n",dim,ie);
#endif
            fe_3d_transport_body_resid(mod,elem_rhs_M,ie,epsilon,i,PERTURB_C,-1,DEBUG);
            
        }
        
        /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        /* calculate residual gradient ++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        
        elem_matrix_deriv_1dof(i, nodes_on_element, elem_rhs_P, elem_rhs_M, elem_mat, epsilon2);
        
    }
    
}
