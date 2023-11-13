/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_hvel_load.c This file collections functions responsible for loading the SW 3D matrix */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void perturb_dpl(SMODEL *mod, int nodes_on_elem, int ie, int dim, DOF_3 *elem_mat_dpl, int DEBUG);
void perturb_vx_vy(SMODEL *mod, int nodes_on_elem, int ie, int dim, DOF_3 *elem_mat_u, DOF_3 *elem_mat_v, int DEBUG);
void print_hvel_mat_to_file(SSUPER_MODEL *sm, int imod);

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
 *  \details Elemental matrix arrays - the "_u", for example, indicates the derivative with respect to u
 *   since these are DOF_3's there is one of these for each equation type.
 *   A 4-node element example the y-velocity derivative of the depth-averaged continuity matrix is: \n
 *  \f{eqnarray*}{ elem\_mat\_v[i].c\_eq =
 *    \begin{bmatrix}
 *    \deriv{R_{c}^0}{v^0} & \deriv{R_{c}^0}{v^1} & \deriv{R_{c}^0}{v^2} & \deriv{R_{c}^0}{v^3} \\
 *    \deriv{R_{c}^1}{v^0} & \deriv{R_{c}^1}{v^1} & \deriv{R_{c}^1}{v^2} & \deriv{R_{c}^1}{v^3} \\
 *    \deriv{R_{c}^2}{v^0} & \deriv{R_{c}^2}{v^1} & \deriv{R_{c}^2}{v^2} & \deriv{R_{c}^2}{v^3} \\
 *    \deriv{R_{c}^3}{v^0} & \deriv{R_{c}^3}{v^1} & \deriv{R_{c}^3}{v^2} & \deriv{R_{c}^3}{v^3}
 *    \end{bmatrix}
 *  \f}
 * \n where i=0 corresponds to \f$ \deriv{R_{c}^0}{v^0} \f$, i=1 to \f$ \deriv{R_{c}^0}{v^1} \f$, etc.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_hvel_load(SSUPER_MODEL *sm, int imod) {
    
    SMODEL *mod = &(sm->submodel[imod]);
#ifdef _DEBUG
    assert(mod->flag.SW3_FLOW == ON);
    assert(mod->nsys == 3);
    assert(mod->nsys_sq == 9);
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
    SSW_3D *sw3 = mod->sw->d3;
    SGRID *grid = mod->grid;
    
    // local variables
    int i, j, k, icol, iend, inode, ie, ie2d = UNSET_INT, ie2d_sur = UNSET_INT, ie2d_bed = UNSET_INT, ie3d = UNSET_INT, GnodeID = UNSET_INT, nnodes = 0;
    int new_column_flag = YES;
    ID_LIST_ITEM *ptr;
    int nmatrix = NDONPRISM * NDONPRISM; // define using max 3D element size
    
    // elemental matrix stored as arrays
    DOF_3 elem_mat_dpl[nmatrix], elem_mat_u[nmatrix], elem_mat_v[nmatrix];

#ifdef _PETSC
    int ierr = 0;
#endif
    
    // initialize matrix
    if (mod->amICoupled == NO) { // if this is a monolothic run, do this in the super load routine
        assert(sm->nsubmodels == 1);
#ifdef _PETSC
        // Zero the PETSc matrix
        ierr = MatZeroEntries(sm->A);
#else
        init_adh_matrix(grid->nnodes, mod->max_nsys_sq, sm->matrix, sm->diagonal);
#endif
    }
    
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 3D ELEMENT CONTRIBUTIONS
    
    for(icol=0; icol<grid->ncolumns; icol++)  {
        new_column_flag = YES;
        
        // loop over the 3d elements in the column
        ptr = grid->column_list[icol];
        while(ptr->next != NULL) {
            ie3d = ptr->id;
            nnodes = grid->elem3d[ie3d].nnodes;
            
#ifdef _DEBUG
            assert(icol == grid->elem3d[ie3d].icol);
#endif
            
            dof3_init_array(elem_mat_dpl, nmatrix);
            dof3_init_array(elem_mat_u,   nmatrix);
            dof3_init_array(elem_mat_v,   nmatrix);
            
            perturb_dpl(mod, grid->elem3d[ie3d].nnodes, ie3d, 3, elem_mat_dpl, DEBUG);
            perturb_vx_vy(mod, grid->elem3d[ie3d].nnodes, ie3d, 3, elem_mat_u, elem_mat_v,DEBUG);
            
            // matrix dirichlet bcs
            for(i=0; i<nnodes; i++) {
                GnodeID = grid->elem3d[ie3d].nodes[i];
                j = mod->nsys * GnodeID;
                if(mod->bc_mask[j] == YES) {
                    fe_assign_mom_db_dof3(i, nnodes, elem_mat_u, elem_mat_v, elem_mat_dpl);
                }
            }
#ifdef _DEBUG
            if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
            if (DEBUG) {
                printf("3d element: %d \n",ie3d);
                dof3_printScreen_array("elem_mat_dpl",elem_mat_dpl, nmatrix,__LINE__,__FILE__);
                dof3_printScreen_array("elem_mat_u",elem_mat_u, nmatrix,__LINE__,__FILE__);
                dof3_printScreen_array("elem_mat_v",elem_mat_v, nmatrix,__LINE__,__FILE__);
            }
#endif
            
            // assemble global Jacobi Matrix
            fe_global_matrix_assemble_sw3(3, grid, mod->nsys, nnodes, ie3d, mod->fmap, elem_mat_u, elem_mat_v, elem_mat_dpl,
#ifdef _PETSC
                    sm->A
#else
                    sm->diagonal, sm->matrix
#endif
                    );
            
            // advance do next element
            new_column_flag = NO;
            ptr = ptr->next;
        }
    }
    
#ifdef _DEBUG
    if (DEBUG) {
        printScreen_matrix("FINAL sw3 matrix 3d element addtion", sm->diagonal, sm->matrix, grid->nnodes, mod->max_nsys_sq, __LINE__, __FILE__);
    }
#endif
    
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 2D ELEMENT CONTRIBUTIONS
    
    for(ie2d = 0; ie2d < grid->nelems2d; ie2d++) {
        nnodes = grid->elem2d[ie2d].nnodes;
        
        dof3_init_array(elem_mat_dpl, nmatrix);
        dof3_init_array(elem_mat_u, nmatrix);
        dof3_init_array(elem_mat_v, nmatrix);
        
        perturb_dpl(mod, grid->elem2d[ie2d].nnodes, ie2d, 2, elem_mat_dpl, DEBUG);
        perturb_vx_vy(mod, grid->elem2d[ie2d].nnodes, ie2d, 2, elem_mat_u, elem_mat_v, DEBUG);
        
        // matrix dirichlet bcs
        for(i = 0; i<grid->elem2d[ie2d].nnodes; i++) {
            GnodeID = grid->elem2d[ie2d].nodes[i];
            j = mod->nsys * GnodeID;
            if(mod->bc_mask[j]==YES) {
                fe_assign_mom_db_dof3(i, grid->elem2d[ie2d].nnodes, elem_mat_u, elem_mat_v, elem_mat_dpl);
            }
        }
#ifdef _DEBUG
        if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
        if (DEBUG) {
            //if (grid->node[0].id == 0 || grid->node[1].id == 0 || grid->node[2].id == 0) {
            printf("2d element: %d\n",ie2d);
            dof3_printScreen_array("elem_mat_dpl",elem_mat_dpl, NELEMAT,__LINE__,__FILE__);
            dof3_printScreen_array("elem_mat_u",elem_mat_u, NELEMAT,__LINE__,__FILE__);
            dof3_printScreen_array("elem_mat_v",elem_mat_v, NELEMAT,__LINE__,__FILE__);
        }
#endif
        fe_global_matrix_assemble_sw3(2, grid, mod->nsys, nnodes, ie2d, mod->fmap, elem_mat_u, elem_mat_v, elem_mat_dpl,
#ifdef _PETSC
                    sm->A
#else
                    sm->diagonal, sm->matrix
#endif
                    );
    }
    
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
   
#ifdef _PETSC
    //ierr = MatAssemblyBegin(sm->A, MAT_FINAL_ASSEMBLY);
    //ierr = MatAssemblyEnd(sm->A, MAT_FINAL_ASSEMBLY);
    
    // Below is for adding small value to zeros on diagonal
    if (mod->amICoupled == NO) { // if this is a monolothic run, do this in the super load routine
        assert(sm->nsubmodels == 1);
        // Check that all diagonal elements assigned (may be zero though)
        //PetscBool missing_diagonal;
        //PetscInt dd;
        //ierr = MatMissingDiagonal(sm->A, &(missing_diagonal), &(dd));
        //if(missing_diagonal){
        //    printf("Missing a diagonal entry on proc %i\n",mod->grid->smpi->myid);
        //}
        int num_diag = grid->my_nnodes*mod->max_nsys;
        //PetscScalar diag_entries[num_diag];
        double small_number = NOT_QUITE_SMALL;
        // Read in the diagonal entries of the matrix
        // TODO: FIX THIS
        // MatGetValues isn't properly indexed for the global matrix
        //printf("Before the get\n");
        //for(i = 0; i < num_diag; i++){
        //    ierr = MatGetValues(sm->A,1,&i,1,&i,&(diag_entries[i]));
        //    printf("d[%i]=%f\n",i,diag_entries[i]);
        //    //ierr = MatGetValuesLocal(sm->A,1,&i,1,&i,&(diag_entries[i]));
        //}
        //printf("Check diagonal correctly read:\n");
        //for(i=0;i<num_diag;i++){
        //    printf("%i %f\n",i,diag_entries[i]);
        //}
        for(i = 0; i < num_diag; i++){
            ierr = MatSetValuesLocal(sm->A,1,&i,1,&i,&small_number,ADD_VALUES);
            //if(fabs(diag_entries[i]) < small_number){
            //    ierr = MatSetValuesLocal(sm->A,1,&i,1,&i,&small_number,ADD_VALUES);
            //}
        }
    }

    // Above is for adding small value to zeros on diagonal
    ierr = MatAssemblyBegin(sm->A, MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(sm->A, MAT_FINAL_ASSEMBLY);
#else
    // checks the diagonal for nonexistent entries
    // SEH Could also be moved to sw2_resid
    // velocity at dry nodes fixed - dss
    if (mod->amICoupled == NO) { // if this is a monolothic run, do this in the super load routine
        assert(sm->nsubmodels == 1);
        double small_number = NOT_QUITE_SMALL;
        for (i = 0; i < grid->nnodes; i++) {
            j = i * mod->max_nsys_sq;
            if (fabs(sm->diagonal[j]) < small_number) {
                sm->diagonal[j] = small_number;
            }
            j = j + 4;
            if (fabs(sm->diagonal[j]) < small_number) {
                sm->diagonal[j] = small_number;
            }
            j = j + 4;
            if (fabs(sm->diagonal[j]) < small_number)
                sm->diagonal[j] = small_number;
        }
#ifdef _MESSG
        comm_update_double(sm->diagonal, mod->max_nsys_sq, mod->grid->smpi); // CJT :: questionable whether needed
#endif
        
#ifdef _DEBUG
        if (DEBUG_MATRIX == ON) {
            printScreen_matrix("sw3d matrix", sm->diagonal, sm->matrix, grid->nnodes, mod->max_nsys_sq, __LINE__, __FILE__);
            tl_error("matrix printing from sw3 hvel load");
        }
#endif
    }
    
#ifdef _DEBUG
    if (DEBUG) {
        time_t time2;  time(&time2);
        TIME_IN_HVEL_LOAD += difftime(time2,time1);
    }
#endif
#endif
    
    
    //printScreen_matrix("sw3d matrix", sm->diagonal, sm->matrix, grid->nnodes, mod->max_nsys_sq, __LINE__, __FILE__);
    //exit(-1);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Adds perturbation to elevation displacement for each node in element for Newton Jacobian calculation.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] elem_mat_dpl stores the Jacobian, elemental matrix for the head perturbation
 *  @param[in] mod a pointer to an AdH model where depth should be updated
 *  @param[in] nodes_on_elem the number of nodes on the element
 *  @param[in] ie the element
 *  @param[in] dim the dimension of the resid
 *  @param[in] DEBUG a debug option
 *
 *  note CJT \:: displacement has already been perturbed and stored in an array
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void perturb_dpl(SMODEL *mod, int nodes_on_elem, int ie, int dim, DOF_3 *elem_mat_dpl, int DEBUG) {
    
#ifdef _DEBUG
    assert(dim == 2 || dim == 3);
#endif
    
    int i, j, node_column_count = 0, GnodeID, GsurfnodeID_2, GsurfnodeID, iseg, id = UNSET_INT;
    
    DOF_3 elem_rhs_dpl_P[nodes_on_elem];    // +dpl perturbation initialized by elem_resid
    DOF_3 elem_rhs_dpl_M[nodes_on_elem];    // -dpl perturbation initialized by elem_resid
    
    double p2 = 2. * mod->perturbation;
    
    int node_in_column_flag[nodes_on_elem];
    int node_already_perturbed_flag[nodes_on_elem];
    sarray_init_value_int(node_already_perturbed_flag, nodes_on_elem, NO);
    
    for(i=0; i<nodes_on_elem; i++) {
        
        if (dim == 3) {GnodeID = mod->grid->elem3d[ie].nodes[i];}
        else {         GnodeID = mod->grid->elem2d[ie].nodes[i];}
        iseg = find_vertical_segment(mod->grid, GnodeID, mod->grid->vertical_hash);
        GsurfnodeID = mod->grid->vertical_list[iseg]->id;
        
        if (node_already_perturbed_flag[i] == NO) {
            
            // flag all nodes underneath (in same column) for perturbation
            sarray_init_value_int(node_in_column_flag, nodes_on_elem, FALSE);
            node_column_count = 0;
            for(j=0; j<nodes_on_elem; j++) {
                if (dim == 3) {
                    id = mod->grid->elem3d[ie].nodes[j];
                } else {
                    id = mod->grid->elem2d[ie].nodes[j];
                }
                iseg = find_vertical_segment(mod->grid, id, mod->grid->vertical_hash);
                GsurfnodeID_2 = mod->grid->vertical_list[iseg]->id;
                
                if(GsurfnodeID == GsurfnodeID_2) {
                    node_in_column_flag[j] = YES;            // flag all nodes in the same column for perturbation
                    node_already_perturbed_flag[j] = YES;    // so we don't increment this node again (do we need both?)
                    node_column_count++;
                }
            }
#ifdef _DEBUG
            // There should be two nodes in a column on a prism element
            if (nodes_on_elem == NDONPRISM && node_column_count != 2) {
                printf("ERROR: count = %d\n",node_column_count);
                tl_error("nodes in element column count are wrong!\n");
            }
#endif
            
            if (dim == 2) { // BOUNDARY
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // (+) boundary perturbation of displacement +++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
                if (DEBUG) printf("\nload :: 2D perturbing +dpl for %dd element: %d\n",dim,ie);
#endif
                fe_hvel_boundary_resid(mod, elem_rhs_dpl_P, ie, 0, i, PERTURB_DPL, +1, node_in_column_flag, DEBUG);
                
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // (-) boundary perturbation of displacement +++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
                if (DEBUG) printf("\nload :: 2D perturbing -dpl for %dd element: %d\n",dim,ie);
#endif
                fe_hvel_boundary_resid(mod, elem_rhs_dpl_M, ie, 0, i, PERTURB_DPL, -1, node_in_column_flag, DEBUG);
                
            } else { // BODY
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // (+) body perturbation of displacement +++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
                if (DEBUG) printf("\nload :: 3D perturbing +dpl for %dd element: %d\n",dim,ie);
#endif
                fe_hvel_body_resid(mod, elem_rhs_dpl_P, ie, 0, i, PERTURB_DPL, +1, node_in_column_flag, DEBUG);
                
                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // (-) body perturbation of displacement +++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
                if (DEBUG) printf("\nload :: 3D perturbing -dpl for %dd element: %d\n",dim,ie);
#endif
                fe_hvel_body_resid(mod, elem_rhs_dpl_M, ie, 0, i, PERTURB_DPL, -1, node_in_column_flag, DEBUG);
            }
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // calculate Jacobian ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            elem_matrix_deriv(i, nodes_on_elem, elem_rhs_dpl_P, elem_rhs_dpl_M, elem_mat_dpl, p2);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // Vel BC residual replacement +++++++++++++++++++++++++++++++++++++++++++++++++
            // CJT DOING THIS IN BODY RESID NOW
            ///*
            if (dim == 3) {
                int irow = UNSET_INT;
                for (j=0; j<nodes_on_elem; j++) {
                    GnodeID = mod->grid->elem3d[ie].nodes[j];
                    if (mod->bc_mask[GnodeID * mod->nsys] == 2) {
                        irow = j * nodes_on_elem + i;
                        elem_mat_dpl[irow].y_eq = mod->sw->d3->tanvec[GnodeID].x*elem_mat_dpl[irow].x_eq + mod->sw->d3->tanvec[GnodeID].y * elem_mat_dpl[irow].y_eq;
                        elem_mat_dpl[irow].x_eq = 0.0;
                    }
                }
            }
            //*/
        }
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Adds perturbation to lateral 3D velocity for each node in element for Newton Jacobian calculation.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] elem_mat_u stores the Jacobian, elemental matrix for the x-velocity perturbation
 *  @param[in,out] elem_mat_v stores the Jacobian, elemental matrix for the y-velocity perturbation
 *  @param[in] mod a pointer to an AdH model where depth should be updated
 *  @param[in] nodes_on_elem the number of nodes on the element
 *  @param[in] ie the element
 *  @param[in] dim the dimension of the resid
 *  @param[in] DEBUG a debug option
 *
 *  \note stuff
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void perturb_vx_vy(SMODEL *mod, int nodes_on_elem, int ie, int dim, DOF_3 *elem_mat_u, DOF_3 *elem_mat_v, int DEBUG) {
    
#ifdef _DEBUG
    assert(dim == 2 || dim == 3);
#endif
    
    int i, j, irow, GnodeID;
    double epsilon = 0., epsilon_u = 0., epsilon_v = 0., epsilon_2u = 0., epsilon_2v = 0.;
    double perturbation = sqrt(SMALL);
    
    DOF_3 elem_rhs_u_P[nodes_on_elem]; /* the residual resulting from a positive (P) perturbation of u */
    DOF_3 elem_rhs_u_M[nodes_on_elem]; /* the residual resulting from a negative (M) perturbation of u */
    DOF_3 elem_rhs_v_P[nodes_on_elem];
    DOF_3 elem_rhs_v_M[nodes_on_elem];
    
    int node_in_column_flag[nodes_on_elem];
    sarray_init_value_int(node_in_column_flag, nodes_on_elem, FALSE);
    
    // get maximum velocity over the element
    SVECT2D max_vel; svect2d_init(&max_vel);
    max_vel.x = 1.;
    max_vel.y = 1.;
    // cjt :: I'm not sure this should be done this way, since it changes the nodes perturbation value depending on what element it is solving
    for (i=0; i<nodes_on_elem; i++) {
        if (dim==3) {
            GnodeID = mod->grid->elem3d[ie].nodes[i];
        } else {
            GnodeID = mod->grid->elem2d[ie].nodes[i];
        }
        if (fabs(mod->sw->d3->vel[GnodeID].x) > max_vel.x) max_vel.x = fabs(mod->sw->d3->vel[GnodeID].x);
        if (fabs(mod->sw->d3->vel[GnodeID].y) > max_vel.y) max_vel.y = fabs(mod->sw->d3->vel[GnodeID].y);
    }
    
    for (i=0; i<nodes_on_elem; i++) {
        if (dim == 2) { // BOUNDARY TERMS
            GnodeID = mod->grid->elem2d[ie].nodes[i];
            NUM_DIFF_EPSILON(epsilon_u, epsilon_2u, mod->sw->d3->vel[GnodeID].x, perturbation);    // calculates epsilon and 2*epsilon
            NUM_DIFF_EPSILON(epsilon_v, epsilon_2v, mod->sw->d3->vel[GnodeID].y, perturbation);    // calculates epsilon and 2*epsilon
            
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) boundary perturbation of u-velocity +++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +u for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_u);
#endif
            fe_hvel_boundary_resid(mod,elem_rhs_u_P,ie,epsilon_u,i,PERTURB_U,+1,node_in_column_flag,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) boundary perturbation of u-velocity +++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -u for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_u);
#endif
            fe_hvel_boundary_resid(mod,elem_rhs_u_M,ie,epsilon_u,i,PERTURB_U,-1,node_in_column_flag,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) boundary perturbation of v-velocity +++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +v for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_v);
#endif
            fe_hvel_boundary_resid(mod,elem_rhs_v_P,ie,epsilon_v,i,PERTURB_V,+1,node_in_column_flag,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) boundary perturbation of v-velocity +++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -v for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_v);
#endif
            fe_hvel_boundary_resid(mod,elem_rhs_v_M,ie,epsilon_v,i,PERTURB_V,-1,node_in_column_flag,DEBUG);
            
        } else if (dim == 3) { // BODY TERMS
            GnodeID = mod->grid->elem3d[ie].nodes[i];
            NUM_DIFF_EPSILON(epsilon_u, epsilon_2u, mod->sw->d3->vel[GnodeID].x, perturbation);    // calculates epsilon and 2*epsilon
            NUM_DIFF_EPSILON(epsilon_v, epsilon_2v, mod->sw->d3->vel[GnodeID].y, perturbation);    // calculates epsilon and 2*epsilon
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) body perturbation of u-velocity +++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +u for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_u);
#endif
            fe_hvel_body_resid(mod,elem_rhs_u_P,ie,epsilon_u,i,PERTURB_U,+1,node_in_column_flag,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) body perturbation of u-velocity +++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -u for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_u);
#endif
            fe_hvel_body_resid(mod,elem_rhs_u_M,ie,epsilon_u,i,PERTURB_U,-1,node_in_column_flag,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) body perturbation of v-velocity +++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +v for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_v);
#endif
            fe_hvel_body_resid(mod,elem_rhs_v_P,ie,epsilon_v,i,PERTURB_V,+1,node_in_column_flag,DEBUG);
            
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) body perturbation of v-velocity +++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -v for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_v);
#endif
            fe_hvel_body_resid(mod,elem_rhs_v_M,ie,epsilon_v,i,PERTURB_V,-1,node_in_column_flag,DEBUG);
            
        } else {
            
            printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
            tl_error(">> dimension (dim) must be either 1 or 2.");
            
        }
        
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // calculate Jacobian ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        elem_matrix_deriv(i, nodes_on_elem, elem_rhs_u_P, elem_rhs_u_M, elem_mat_u, epsilon_2u);
        elem_matrix_deriv(i, nodes_on_elem, elem_rhs_v_P, elem_rhs_v_M, elem_mat_v, epsilon_2v);
        
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Vel BC residual replacement +++++++++++++++++++++++++++++++++++++++++++++++++
        // CJT :: DOING THIS IN BODY RESID NOW
        ///*
        if (dim == 3) {
            for (j=0; j<nodes_on_elem; j++) {
                GnodeID = mod->grid->elem3d[ie].nodes[j];
                if (mod->bc_mask[ GnodeID * mod->nsys] == 2) {
                    irow = j * nodes_on_elem + i;
                    elem_mat_u[irow].y_eq = mod->sw->d3->tanvec[GnodeID].x * elem_mat_u[irow].x_eq + mod->sw->d3->tanvec[GnodeID].y * elem_mat_u[irow].y_eq;
                    elem_mat_u[irow].x_eq = 0.;
                    elem_mat_v[irow].y_eq = mod->sw->d3->tanvec[GnodeID].x * elem_mat_v[irow].x_eq + mod->sw->d3->tanvec[GnodeID].y * elem_mat_v[irow].y_eq;
                    elem_mat_v[irow].x_eq = 0.;
                }
            }
        }
        //*/
    }
}




void print_hvel_mat_to_file(SSUPER_MODEL *sm, int imod) {
#ifndef _PETSC
    
    SMODEL *mod = &(sm->submodel[imod]);

    int i = 0, j = 0, k = 0, l = 0, col = 0;
    FILE *fp = fopen("hvel_umf_matrix.txt","w");

    printf("hvel matrix:\n");
    for(i = 0; i < mod->grid->nnodes; i++){
        // Print matrix row diagonal
        for(j = 0; j < mod->nsys_sq; j++){
            fprintf(fp,"%i %i %.16e\n",i*3+(j/3)+1,i*3+(j%3)+1,sm->diagonal[i*mod->nsys_sq + j]);
        }
        for(j = 0; j < sm->matrix[i].size; j++){
            col = sm->matrix[i].index[j];
            for(k = j*mod->nsys_sq; k < (j+1)*mod->nsys_sq; k++){
                l = k - j*mod->nsys_sq;
                fprintf(fp,"%i %i %.16e\n",i*3+(l/3)+1,col*3+(l%3)+1,sm->matrix[i].value[k]);
            }
        }

    }
    fclose(fp);

#endif
}
