/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_sw2_load.c This file collections functions responsible for loading the SW 2D matrix */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"
static int matrix_increment = 0; // prints each matrix to a separate file for comparing

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void perturb_gw_head(SMODEL *mod, int nnodes, int ie, int dim, double *elem_mat, int DEBUG);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Forms the GW Newton Jacobi matrix.
 *  \author    Matthew Farthing, Ph.D.
 *  \author    Gajanan Choudhary, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out]  mod a pointer to an AdH model where depth should be updated
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#ifndef _DEBUG
//#define _DEBUG
#define _DEBUG_LOCAL
#endif
void fe_gw_load(SSUPER_MODEL *sm, int imod) {
#ifndef _PETSC

    SMODEL *mod = &(sm->submodel[imod]);

#ifdef _DEBUG
    assert(mod->flag.GW_FLOW == ON);
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
    time_t time1;  time(&time1);
#endif

    // aliases
    SGW *gw = mod->sgw;
    SGRID *grid = mod->grid;
    
    // local variable declarations
    int i = 0, j = 0, ie = UNSET_INT;
    int nmatrix = MAX_NNODES_ON_ELEM3D * MAX_NNODES_ON_ELEM3D; // define using max element size 
    double elem_mat[nmatrix];
    int nodes_per_elem = 0;

    // initialize matrix
    if (mod->amICoupled == NO) { // if this is a monolothic run, do this in the super load routine
        init_adh_matrix(grid->nnodes, mod->max_nsys_sq, sm->matrix, sm->diagonal);
    }

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 3D ELEMENT CONTRIBUTIONS
    for (ie = 0; ie < grid->nelems3d; ie++) {
        sarray_init_dbl(elem_mat, nmatrix);
        perturb_gw_head(mod, grid->elem3d[ie].nnodes, ie, 3, elem_mat, DEBUG);

        //enforce Dirichlet boundary conditions
        nodes_per_elem = grid->elem3d[ie].nnodes;
        for (i=0; i < nodes_per_elem; i++) {
            int GlobalNodeID = grid->elem3d[ie].nodes[i];
            int istr = mod->grid->node[GlobalNodeID].string;
            if (istr > NORMAL) {
                int index = 0;
                if (mod->str_values[istr].flow.bc_flag == BCT_PRS_DIR ||
                        mod->str_values[istr].flow.bc_flag == BCT_DIR) {
                    /* Zero the row */
                    for (j=0; j < nodes_per_elem; j++) { //zero row
                        index = i * nodes_per_elem + j;

                        /*mwf debug
                        printf("zeroing elem %d node %d istr %d j= %d index= %d was%g \n",ie,i,istr,j,index,elem_mat[index]);
                        */
                        elem_mat[index] = 0.;
                    }
                    /* Zero the column - Gajanan gkc adding based on fe_wvel_load.c */
                    for(j=i; j<nodes_per_elem*nodes_per_elem; j+=nodes_per_elem) {
                        elem_mat[j] = 0.;
                    }

                    //diagonal
                    index = i * nodes_per_elem + i;
                    elem_mat[index] = 1.0;
                    /*mwf debug
                    printf("setting elem %d node %d diagonal istr %d index= %d val %g \n",ie,i,istr,index,elem_mat[index]);
                    */
                }
            }
        } //end Dirichlet enforcement
        /*mwf debug 
        printf("3d element: %d djac= %g \n",ie,grid->elem3d[ie].djac);
        sarray_printScreen_dbl(elem_mat, nmatrix, "elem_mat");
        */

#ifdef _DEBUG
        if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
        if (DEBUG)  {
            printf("3d element: %d\n",ie);
            sarray_printScreen_dbl(elem_mat, nmatrix, "elem_mat");
        }
#endif
        // assembles the element contributions
        fe_global_matrix_assemble_scalar(grid, mod->nsys, grid->elem3d[ie].nnodes, grid->elem3d[ie].nodes, mod->fmap, elem_mat, sm->diagonal, sm->matrix
#ifdef _PETSC
                , sm->A
#endif
                );
    }

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 2D ELEMENT CONTRIBUTIONS
    for (ie = 0; ie < grid->nelems2d; ie++) {
        sarray_init_dbl(elem_mat, nmatrix);
        perturb_gw_head(mod, grid->elem2d[ie].nnodes, ie, 2, elem_mat, DEBUG);

        nodes_per_elem = grid->elem2d[ie].nnodes;

#ifdef _DEBUG
        if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
        if (DEBUG)  {
            printf("2d element: %d\n",ie);
            sarray_printScreen_dbl(elem_mat, nmatrix, "elem_mat");
        }
#endif
        // assembles the element contributions
        fe_global_matrix_assemble_scalar(grid, mod->nsys, grid->elem2d[ie].nnodes, grid->elem2d[ie].nodes, mod->fmap, elem_mat, sm->diagonal, sm->matrix
#ifdef _PETSC
                , sm->A
#endif
                );
    }

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    if (mod->amICoupled == NO) { // if this is a monolothic run, do this in the super load routine
        /* Gajanan gkc commenting out the following lines. This is because GW diagonal seems to have
         * many entries smaller than 1e-7, which get replaced to NOT_QUITE_SMALL due to the following
         * function call. This causes the shock speed in the celia test case to differ from 
         * that obtained in an earlier AdH repo - fawkes v5.0.4.
         * */
        // //checks the diagonal for nonexistent entries
        // check_matrix_diagonal_for_nonexistant_entries(grid->nnodes, mod->max_nsys_sq, mod->nsys, sm->diagonal);
#ifdef _MESSG
        comm_update_double(sm->diagonal, 1, mod->grid->smpi);
#endif

#ifdef _DEBUG
        if (DEBUG_MATRIX == ON
#ifdef _MESSG
           && grid->smpi->myid == 0
#endif
           ) {
            printScreen_matrix("groundwater matrix", sm->diagonal, sm->matrix, grid->nnodes, mod->max_nsys_sq, __LINE__, __FILE__);
            //exit(-1);
        }
#endif
    }
    
    
#ifdef _DEBUG // timings
    time_t time2;  time(&time2);
    TIME_IN_GW_LOAD += difftime(time2,time1);
#endif

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
void perturb_gw_head(SMODEL *mod, int nodes_on_element, int ie, int dim, double *elem_mat, int DEBUG) {
          
    assert(dim == 2 || dim == 3);
    assert(nodes_on_element == 3 || nodes_on_element == 4);
          
    int i,j;
    double epsilon = 0., epsilon2 = 0.;
    double perturbation = sqrt(SMALL);

    double elem_rhs_h_P[nodes_on_element];    // +head perturbation
    double elem_rhs_h_M[nodes_on_element];    // -head perturbation
          
    // get perturbation scale over element
    int GlobalNodeID = UNSET_INT; int node_string = UNSET_INT; 
    int dirichlet_node = NO;
    double max_h = 0.;
    for (i=0; i<nodes_on_element; i++) {
        if (dim == 2) {
            GlobalNodeID = mod->grid->elem2d[ie].nodes[i];
        } else {
            GlobalNodeID = mod->grid->elem3d[ie].nodes[i];
        }
        // *****************************************************************************
        // perturbation of psi ***********************************************************
        NUM_DIFF_EPSILON(epsilon, epsilon2, mod->sgw->gw_phead[GlobalNodeID], perturbation);    // calculates epsilon and 2*epsilon
        /*mwf match AdH fawkes version for now? */
        epsilon = NOT_QUITE_SMALL;
        epsilon2= 2.0*epsilon;


        node_string = mod->grid->node[GlobalNodeID].string;
        dirichlet_node = NO;
        if (node_string > NORMAL) {
            if (mod->str_values[node_string].flow.bc_flag == BCT_PRS_DIR || mod->str_values[node_string].flow.bc_flag == BCT_DIR) {
                dirichlet_node = YES;
            }
        }
        if (dim == 2) {

#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +psi for %dd element %d\n",dim,ie);
#endif
            fe_gw_boundary_resid(mod,elem_rhs_h_P,ie,epsilon,i,PERTURB_H,+1,DEBUG);

#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -psi for %dd element %d\n",dim,ie);
#endif
            fe_gw_boundary_resid(mod,elem_rhs_h_M,ie,epsilon,i,PERTURB_H,-1,DEBUG);

            elem_matrix_deriv_1dof(i, nodes_on_element, elem_rhs_h_P, elem_rhs_h_M, elem_mat, epsilon2); // gradient

        } else if (dim == 3) { /* 3D */

#ifdef _DEBUG
            if (DEBUG) printf("\ngroundwater perturbation size: %30.20f\n",epsilon);
#endif
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +psi for %dd element %d\n",dim,ie);
#endif
            fe_gw_body_resid(mod,elem_rhs_h_P,ie,epsilon,i,PERTURB_H,+1,DEBUG);

#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -psi for %dd element %d\n",dim,ie);
#endif
            fe_gw_body_resid(mod,elem_rhs_h_M,ie,epsilon,i,PERTURB_H,-1,DEBUG);

            /*mwf debug
            printf("\ngroundwater perturbation= %30.20f  double perturbation size: %30.20f\n",epsilon,epsilon2);
            printf("\nload :: i= %d pertubing +psi for %d element %d\n",i,dim,ie);
            for (j=0; j < nodes_on_element; j++) {
                printf("elem_rhs_P[%d]= %30.20f ",j,elem_rhs_h_P[j]);
            }
            printf("\nload :: i= %d pertubing -psi for %d element %d\n",i,dim,ie);
            for (j=0; j < nodes_on_element; j++) {
                printf("elem_rhs_M[%d]= %30.20f ",j,elem_rhs_h_M[j]);
            }
            */

            elem_matrix_deriv_1dof(i, nodes_on_element, elem_rhs_h_P, elem_rhs_h_M, elem_mat, epsilon2); // gradient
        } /* end switch on dirichlet nodes */
    }/* end first node loop */

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
void perturb_gw_head_skip_dirichlet(SMODEL *mod, int nodes_on_element, int ie, int dim, double *elem_mat, int DEBUG) {

    assert(dim == 2 || dim == 3);
    assert(nodes_on_element == 3 || nodes_on_element == 4);

    int i,j;
    double epsilon = 0., epsilon2 = 0.;
    double perturbation = sqrt(SMALL);

    double elem_rhs_h_P[nodes_on_element];    // +head perturbation
    double elem_rhs_h_M[nodes_on_element];    // -head perturbation

    // get perturbation scale over element
    int GlobalNodeID = UNSET_INT; int node_string = UNSET_INT; 
    int dirichlet_node = NO;
    double max_h = 0.;
    for (i=0; i<nodes_on_element; i++) {
        if (dim == 2) {
            GlobalNodeID = mod->grid->elem2d[ie].nodes[i];
        } else {
            GlobalNodeID = mod->grid->elem3d[ie].nodes[i];
        }
        // *****************************************************************************
        // perturbation of psi ***********************************************************
        NUM_DIFF_EPSILON(epsilon, epsilon2, mod->sgw->gw_phead[GlobalNodeID], perturbation);    // calculates epsilon and 2*epsilon
        /*mwf match AdH fawkes version for now? */
        epsilon = NOT_QUITE_SMALL;
        epsilon2= 2.0*epsilon;


        node_string = mod->grid->node[GlobalNodeID].string;
        dirichlet_node = NO;
        if (node_string > NORMAL) {
            if (mod->str_values[node_string].flow.bc_flag == BCT_PRS_DIR || mod->str_values[node_string].flow.bc_flag == BCT_DIR) {
                dirichlet_node = YES;
            }
        }
        if (dim == 2 && dirichlet_node == NO) {

#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +psi for %dd element %d\n",dim,ie);
#endif
            fe_gw_boundary_resid(mod,elem_rhs_h_P,ie,epsilon,i,PERTURB_H,+1,DEBUG);

#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -psi for %dd element %d\n",dim,ie);
#endif
            fe_gw_boundary_resid(mod,elem_rhs_h_M,ie,epsilon,i,PERTURB_H,-1,DEBUG);

            elem_matrix_deriv_1dof(i, nodes_on_element, elem_rhs_h_P, elem_rhs_h_M, elem_mat, epsilon2); // gradient

        } else if (dim == 3 && dirichlet_node == NO) { /* 3D */

#ifdef _DEBUG
            if (DEBUG) printf("\ngroundwater perturbation size: %30.20f\n",epsilon);
#endif
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +psi for %dd element %d\n",dim,ie);
#endif
            fe_gw_body_resid(mod,elem_rhs_h_P,ie,epsilon,i,PERTURB_H,+1,DEBUG);

#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -psi for %dd element %d\n",dim,ie);
#endif
            fe_gw_body_resid(mod,elem_rhs_h_M,ie,epsilon,i,PERTURB_H,-1,DEBUG);

            /*mwf debug
            printf("\ngroundwater perturbation= %30.20f  double perturbation size: %30.20f\n",epsilon,epsilon2);
            printf("\nload :: i= %d pertubing +psi for %d element %d\n",i,dim,ie);
            for (j=0; j < nodes_on_element; j++) {
                printf("elem_rhs_P[%d]= %30.20f ",j,elem_rhs_h_P[j]);
            }
            printf("\nload :: i= %d pertubing -psi for %d element %d\n",i,dim,ie);
            for (j=0; j < nodes_on_element; j++) {
                printf("elem_rhs_M[%d]= %30.20f ",j,elem_rhs_h_M[j]);
            }
            */

            elem_matrix_deriv_1dof(i, nodes_on_element, elem_rhs_h_P, elem_rhs_h_M, elem_mat, epsilon2); // gradient
        } /* end switch on dirichlet nodes */
    }/* end first node loop */

    /* make sure element matrix has row and col zeroed */
    for (i=0; i<nodes_on_element; i++) {
        if (dim == 2) {
            GlobalNodeID = mod->grid->elem2d[ie].nodes[i];
        } else {
            GlobalNodeID = mod->grid->elem3d[ie].nodes[i];
        }
        node_string = mod->grid->node[GlobalNodeID].string;
        dirichlet_node = NO;
        if (node_string > NORMAL) {
            if (mod->str_values[node_string].flow.bc_flag == BCT_PRS_DIR || mod->str_values[node_string].flow.bc_flag == BCT_DIR) {
                dirichlet_node = YES;
            }
        }
        if (dirichlet_node == YES) {
            /* zero row */
            for (j = 0; j < nodes_on_element; j++) {
                elem_mat[i * nodes_on_element + j] = 0.;
            }
            /*zero column -- increment for Dirchlet nodes should always be zero -- */
            for (j=i; j < nodes_on_element*nodes_on_element; j+= nodes_on_element) {
                elem_mat[j] = 0.0;
            }
            elem_mat[i + i * nodes_on_element] = 1.;
        }
    } /* end modfications for dirichlet boundary */
}


#ifdef _DEBUG_LOCAL
#undef _DEBUG
#endif

