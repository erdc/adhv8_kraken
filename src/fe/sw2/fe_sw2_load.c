/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_sw2_load.c This file collections functions responsible for loading the SW 2D matrix */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"
static int matrix_increment = 0; // prints each matrix to a separate file for comparing

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void perturb_sw2_head(SMODEL *mod, int nnodes, int ie, int dim, DOF_3 *elem_mat_h, int DEBUG);
void perturb_sw2_vel( SMODEL *mod, int nnodes, int ie, int dim, DOF_3 *elem_mat_u, DOF_3 *elem_mat_v, int DEBUG);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Forms the SW2D Newton Jacobi matrix.
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
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_sw2_load(SSUPER_MODEL *sm, int imod) {
    
    SMODEL *mod = &(sm->submodel[imod]);
    
#ifdef _DEBUG
    assert(mod->flag.SW2_FLOW == ON);
    assert(mod->nsys == 3);
    assert(mod->nsys_sq == 9);
    assert(sm->diagonal != NULL);
    assert(sm->matrix != NULL);
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
    SSW_2D *sw2 = mod->sw->d2;
    SGRID *grid = mod->grid;
    
    // local variable declarations
    int i = 0, j = 0, ie = UNSET_INT, GlobalNodeID = UNSET_INT, node_string = UNSET_INT, nnodes = 0;
    int nmatrix = NDONQUAD * NDONQUAD; // define using max 2D element size (quadrilateral)
    DOF_3 elem_mat_u[nmatrix], elem_mat_v[nmatrix], elem_mat_h[nmatrix];
#ifdef _PETSC
    int ierr = 0;
#endif
    
    // initialize matrix
    if (mod->amICoupled == NO) { // if this is a monolothic run, do this in the super load routine
#ifdef _PETSC
        // Zero the PETSc matrix
        ierr = MatZeroEntries(sm->A);   // TODO: DO WE NEED THIS - SAM
#else
        init_adh_matrix(grid->nnodes, mod->max_nsys_sq, sm->matrix, sm->diagonal);
#endif
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 2D ELEMENT CONTRIBUTIONS
    for (ie = 0; ie < grid->nelems2d; ie++) {
        if (grid->wd_flag[ie] == 0 || grid->wd_flag[ie] == 1 || grid->wd_flag[ie] == 2) {
            
            if (grid->elem2d[ie].string != UNSET_INT) {
                if (mod->str_values[grid->elem2d[ie].string].phys_flag == SW2_FLAG) {
                    nnodes = grid->elem2d[ie].nnodes;
                    
                    // perturb solutions
                    dof3_init_array(elem_mat_h, nmatrix);
                    dof3_init_array(elem_mat_u, nmatrix);
                    dof3_init_array(elem_mat_v, nmatrix);
                    perturb_sw2_head(mod, nnodes, ie, 2, elem_mat_h, DEBUG);
                    perturb_sw2_vel( mod, nnodes, ie, 2, elem_mat_u, elem_mat_v, DEBUG);
                    
                    // enforce Dirichlet boundary conditions
                    for (i=0; i<nnodes; i++) {
                        
                        GlobalNodeID = grid->elem2d[ie].nodes[i];
                        node_string = grid->node[i].string;
                        
                        // enforce Dirchlet condition for nodes with depth less than zero
                        if (mod->sw->d2->old_head[GlobalNodeID] <= 0.) {
                            fe_assign_mom_db_dof3(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_h);
                        }
                        
                        if (node_string > NORMAL) {
                            if (mod->str_values[node_string].ol_flow.bc_flag == BCT_PRS_DIR) {
                                fe_assign_dacont_db_dof3(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_h);
                            } else if (mod->str_values[node_string].ol_flow.bc_flag == BCT_VEL_DIR) {
                                fe_assign_mom_db_dof3(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_h);
                            } else if (mod->str_values[node_string].ol_flow.bc_flag == BCT_VEL_PRS_DIR) {
                                fe_assign_dacont_db_dof3(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_h);
                                fe_assign_mom_db_dof3(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_h);
                            }
                        }
                    }
#ifdef _DEBUG
                    if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
                    if (DEBUG) {
                        printf("2d element: %d\n", ie);
                        dof3_printScreen_array("elem_mat_h",elem_mat_h, NDONTRI * NDONTRI,__LINE__,__FILE__);
                        dof3_printScreen_array("elem_mat_u",elem_mat_u, NDONTRI * NDONTRI,__LINE__,__FILE__);
                        dof3_printScreen_array("elem_mat_v",elem_mat_v, NDONTRI * NDONTRI,__LINE__,__FILE__);
                    }
#endif
                    // assemble global Jacobi matrix
                    fe_global_matrix_assemble_sw2(grid, mod->nsys, nnodes, grid->elem2d[ie].nodes, mod->fmap,
                                                  elem_mat_u, elem_mat_v, elem_mat_h,
#ifdef _PETSC
                                                  sm->A, sm->Istart, sm->Iend, sm->ownership_range
#else
                                                  sm->diagonal, sm->matrix
#endif
                                                  );
                }
            }
        }
    }
    
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 1D ELEMENT CONTRIBUTIONS
    
    for (ie = 0; ie < grid->nelems1d; ie++) {
        if (grid->elem1d[ie].string != UNSET_INT) {
            if (mod->str_values[grid->elem1d[ie].string].phys_flag == SW2_FLAG) {
                nnodes = grid->elem1d[ie].nnodes;
                
                // perturb solutions
                dof3_init_array(elem_mat_h, nmatrix);
                dof3_init_array(elem_mat_u, nmatrix);
                dof3_init_array(elem_mat_v, nmatrix);
                perturb_sw2_head(mod, nnodes, ie, 1, elem_mat_h, DEBUG);
                perturb_sw2_vel( mod, nnodes, ie, 1, elem_mat_u, elem_mat_v, DEBUG);
                
                // enforce Dirichlet boundary conditions
                for (i=0; i<nnodes; i++) {
                    
                    GlobalNodeID = grid->elem1d[ie].nodes[i];
                    node_string = grid->node[i].string;
                    
                    // enforce Dirchlet condition for nodes with depth less than zero
                    if (mod->sw->d2->old_head[GlobalNodeID] <= 0.) { // CJT :: for elem2d this is <=
                        fe_assign_mom_db_dof3(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_h);
                    }
                    
                    if (node_string > NORMAL) {
                        if (mod->str_values[node_string].ol_flow.bc_flag == BCT_PRS_DIR) {
                            fe_assign_dacont_db_dof3(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_h);
                        } else if (mod->str_values[node_string].ol_flow.bc_flag == BCT_VEL_DIR) {
                            fe_assign_mom_db_dof3(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_h);
                        } else if (mod->str_values[node_string].ol_flow.bc_flag == BCT_VEL_PRS_DIR) {
                            fe_assign_dacont_db_dof3(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_h);
                            fe_assign_mom_db_dof3(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_h);
                        }
                    }
                }
#ifdef _DEBUG
                if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
                if (DEBUG) {
                    printf("1d element: %d\n",ie);
                    dof3_printScreen_array("elem_mat_h",elem_mat_h, NDONTRI * NDONTRI,__LINE__,__FILE__);
                    dof3_printScreen_array("elem_mat_u",elem_mat_u, NDONTRI * NDONTRI,__LINE__,__FILE__);
                    dof3_printScreen_array("elem_mat_v",elem_mat_v, NDONTRI * NDONTRI,__LINE__,__FILE__);
                }
#endif
                // assemble global Jacobi matrix
                fe_global_matrix_assemble_sw2(grid, mod->nsys, nnodes, grid->elem1d[ie].nodes, mod->fmap,elem_mat_u, elem_mat_v, elem_mat_h,
#ifdef _PETSC
                                                  sm->A, sm->Istart, sm->Iend, sm->ownership_range
#else
                                                  sm->diagonal, sm->matrix
#endif
                                                  );
            }
        }
    }
    
    //printScreen_matrix("sw2d matrix", sm->diagonal, sm->matrix, grid->nnodes, mod->max_nsys_sq, __LINE__, __FILE__);
    //exit(-1);
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
#ifdef _PETSC
    // Assemble the PETSc matrix containing the Jacobian
    ierr = MatAssemblyBegin(sm->A, MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(sm->A, MAT_FINAL_ASSEMBLY);
#else
    /* checks the diagonal for nonexistent entries */
    /* SEH Could also be moved to sw2_resid */
    /* velocity at dry nodes fixed - dss */
    if (mod->amICoupled == NO) { // if this is a monolothic run, do this in the super load routine
        for (i = 0; i < grid->nnodes; i++) {
            j = i * mod->max_nsys_sq;
            if (fabs(sm->diagonal[j]) < NOT_QUITE_SMALL) {
                sm->diagonal[j] = NOT_QUITE_SMALL;
                if (fabs(sm->residual[3 * i]) < NOT_QUITE_SMALL)
                    sm->residual[3 * i] = -NOT_QUITE_SMALL * sw2->vel[i].x;
            }
            j = j + 4;
            if (fabs(sm->diagonal[j]) < NOT_QUITE_SMALL) {
                sm->diagonal[j] = NOT_QUITE_SMALL;
                if (fabs(sm->residual[3 * i + 1]) < NOT_QUITE_SMALL)
                    sm->residual[3 * i + 1] = -NOT_QUITE_SMALL * sw2->vel[i].y;
            }
            j = j + 4;
            if (fabs(sm->diagonal[j]) < NOT_QUITE_SMALL)
                sm->diagonal[j] = NOT_QUITE_SMALL;
            
        }
#ifdef _MESSG
        comm_update_double(sm->diagonal, mod->max_nsys_sq, mod->grid->smpi); // CJT :: questionable whether needed
#endif
#ifdef _DEBUG // print matrix to either screen or file
        if (DEBUG_MATRIX == ON
#ifdef _MESSG
            && grid->smpi->myid == 0
#endif
            ) {
            
            // print matrix to file
            char buffer[32];
            FILE *ffile;
            matrix_increment++;
            snprintf(buffer, sizeof(char) * 32, "matrix%i.txt", matrix_increment);
            ffile = fopen(buffer, "wb");
            printFile_matrix(ffile,"sw2d matrix", sm->diagonal, sm->matrix, grid->nnodes, mod->max_nsys_sq, __LINE__, __FILE__);
            fclose(ffile);
            if (matrix_increment == 2) exit(-1);
            
            // print matrix to screen
            //printScreen_matrix("sw2d matrix", sm->diagonal, sm->matrix, grid->nnodes, mod->max_nsys_sq, __LINE__, __FILE__);
            //exit(-1);
        }
#endif
    }
    
#ifdef _DEBUG // timings
    time_t time2;  time(&time2);
    TIME_IN_2D_SW_LOAD += difftime(time2,time1);
#endif
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
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
 *  @param[in,out] elem_mat_h stores the Jacobian, elemental matrix for the head perturbation
 *  @param[in]  mod a pointer to an AdH model where depth should be updated
 *  @param[in] nodes_on_element the number of nodes on the element
 *  @param[in] ie the element
 *  @param[in] dim the dimension of the resid
 *  @param[in] DEBUG a debug option
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void perturb_sw2_head(SMODEL *mod, int nodes_on_element, int ie, int dim, DOF_3 *elem_mat_h, int DEBUG) {
    
#ifdef _DEBUG
    assert(dim == 1 || dim == 2);
#endif
    
    int i,j, GlobalNodeID = UNSET_INT;
    double epsilon = 0., epsilon2 = 0.;
    double perturbation = sqrt(SMALL);
    DOF_3 elem_rhs_h_P[nodes_on_element];    // +head perturbation initialized by elem_resid
    DOF_3 elem_rhs_h_M[nodes_on_element];    // -head perturbation initialized by elem_resid
    
    //DEBUG = ON;
    
    for (i=0; i<nodes_on_element; i++) {
        if (dim == 1) {
            GlobalNodeID = mod->grid->elem1d[ie].nodes[i];
            NUM_DIFF_EPSILON(epsilon, epsilon2, mod->sw->d2->head[GlobalNodeID], perturbation);    // calculates epsilon and 2*epsilon
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) boundary perturbation of depth ++++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: boundary pertubing +h for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon);
#endif
            fe_sw2_boundary_resid(mod,elem_rhs_h_P,ie,epsilon,i,PERTURB_H,+1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) boundary perturbation of depth ++++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: boundary pertubing -h for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon);
#endif
            fe_sw2_boundary_resid(mod,elem_rhs_h_M,ie,epsilon,i,PERTURB_H,-1,DEBUG);
            
        } else if (dim == 2) {
            GlobalNodeID = mod->grid->elem2d[ie].nodes[i];
            NUM_DIFF_EPSILON(epsilon, epsilon2, mod->sw->d2->head[GlobalNodeID], perturbation);    // calculates epsilon and 2*epsilon
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) body perturbation of depth ++++++++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: body pertubing +h for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon);
#endif
            fe_sw2_body_resid(mod,elem_rhs_h_P,ie,epsilon,i,PERTURB_H,+1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) body perturbation of depth ++++++++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: body pertubing -h for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon);
#endif
            fe_sw2_body_resid(mod,elem_rhs_h_M,ie,epsilon,i,PERTURB_H,-1,DEBUG);
            
        } else {
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> dimension (dim) must be either 1 or 2.");
        }
        
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // calculate Jacobian ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        elem_matrix_deriv(i, nodes_on_element, elem_rhs_h_P, elem_rhs_h_M, elem_mat_h, epsilon2); // gradient
    }
    
    //DEBUG =  OFF;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Adds perturbation to depth-averaged velocities for each node in element for Newton Jacobian calculation.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] elem_mat_u stores the Jacobian, elemental matrix for the x-velocity perturbation
 *  @param[in,out] elem_mat_v stores the Jacobian, elemental matrix for the y-velocity perturbation
 *  @param[in]  mod a pointer to an AdH model where depth should be updated
 *  @param[in] ie the element
 *  @param[in] dim the dimension of the resid
 *  @param[in] DEBUG a debug option
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void perturb_sw2_vel(SMODEL *mod, int nodes_on_element, int ie, int dim, DOF_3 *elem_mat_u, DOF_3 *elem_mat_v, int DEBUG) {
    
#ifdef _DEBUG
    assert(dim == 1 || dim == 2);
#endif
    
    int i, GnodeID = UNSET_INT;
    double epsilon_u = 0., epsilon_v = 0., epsilon_2u = 0., epsilon_2v = 0.;
    double perturbation = sqrt(SMALL);
    DOF_3 elem_rhs_u_P[nodes_on_element];
    DOF_3 elem_rhs_u_M[nodes_on_element];
    DOF_3 elem_rhs_v_P[nodes_on_element];
    DOF_3 elem_rhs_v_M[nodes_on_element];
    
    for (i=0; i<nodes_on_element; i++) {
        if (dim == 1) {
#ifdef _DEBUG
            assert(i<NDONSEG);
#endif
            GnodeID = mod->grid->elem1d[ie].nodes[i];
            NUM_DIFF_EPSILON(epsilon_u, epsilon_2u, mod->sw->d2->vel[GnodeID].x, perturbation);    // calculates epsilon and 2*epsilon
            NUM_DIFF_EPSILON(epsilon_v, epsilon_2v, mod->sw->d2->vel[GnodeID].y, perturbation);    // calculates epsilon and 2*epsilon
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) boundary perturbation of u-velocity +++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +u for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_u);
#endif
            fe_sw2_boundary_resid(mod,elem_rhs_u_P,ie,epsilon_u,i,PERTURB_U,+1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) boundary perturbation of u-velocity +++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -u for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_u);
#endif
            fe_sw2_boundary_resid(mod,elem_rhs_u_M,ie,epsilon_u,i,PERTURB_U,-1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) boundary perturbation of v-velocity +++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +v for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_v);
#endif
            fe_sw2_boundary_resid(mod,elem_rhs_v_P,ie,epsilon_v,i,PERTURB_V,+1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) boundary perturbation of v-velocity +++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -v for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_v);
#endif
            fe_sw2_boundary_resid(mod,elem_rhs_v_M,ie,epsilon_v,i,PERTURB_V,-1,DEBUG);
            
        } else if (dim == 2) {
            
            GnodeID = mod->grid->elem2d[ie].nodes[i];
            NUM_DIFF_EPSILON(epsilon_u, epsilon_2u, mod->sw->d2->vel[GnodeID].x, perturbation);    // calculates epsilon and 2*epsilon
            NUM_DIFF_EPSILON(epsilon_v, epsilon_2v, mod->sw->d2->vel[GnodeID].y, perturbation);    // calculates epsilon and 2*epsilon
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) body perturbation of u-velocity +++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +u for %dd element %d :: perturbation: %20.10e :: u: %20.10f\n",dim,ie,epsilon_u,mod->sw->d2->vel[GnodeID].x);
#endif
            fe_sw2_body_resid(mod,elem_rhs_u_P,ie,epsilon_u,i,PERTURB_U,+1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) body perturbation of u-velocity +++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -u for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_u);
#endif
            fe_sw2_body_resid(mod,elem_rhs_u_M,ie,epsilon_u,i,PERTURB_U,-1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) body perturbation of v-velocity +++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +v for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_v);
#endif
            fe_sw2_body_resid(mod,elem_rhs_v_P,ie,epsilon_v,i,PERTURB_V,+1,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) body perturbation of v-velocity +++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -v for %dd element %d :: perturbation: %20.10e :: v: %20.10f\n",dim,ie,epsilon_v,mod->sw->d2->vel[GnodeID].y);
#endif
            fe_sw2_body_resid(mod,elem_rhs_v_M,ie,epsilon_v,i,PERTURB_V,-1,DEBUG);
            
        } else {
            
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> dimension (dim) must be either 1 or 2.");
            
        }
        
        
        if (DEBUG) {
            //printf("\n");
            //printf("ndim: %d \t eps :: u: %20.10e \t v: %20.20f\n",dim,epsilon_2u,epsilon_2v);
            //dof3_printScreen_array("elem_rhs_v_P",elem_rhs_v_P, nodes_on_element,__LINE__,__FILE__);
            //dof3_printScreen_array("elem_rhs_v_M",elem_rhs_v_M, nodes_on_element,__LINE__,__FILE__);
            //exit(-2);
        }
        
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // calculate Jacobian ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        elem_matrix_deriv(i, nodes_on_element, elem_rhs_u_P, elem_rhs_u_M, elem_mat_u, epsilon_2u);
        elem_matrix_deriv(i, nodes_on_element, elem_rhs_v_P, elem_rhs_v_M, elem_mat_v, epsilon_2v);
    }
}

