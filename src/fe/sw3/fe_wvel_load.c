/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_wvel_load.c This file collections functions for loading the 3D continuity matrix */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void perturb_w(SMODEL *mod, int nnodes, int ie, int dim, double *elem_mat, int DEBUG);
void print_wvel_mat_to_file(SSUPER_MODEL *sm, int imod);

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
 *  \f{eqnarray*}{ \weakSwContNoKinematic{i} \f}
 *  \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_wvel_load(SSUPER_MODEL *sm, int imod) {
    
    SMODEL *mod = &(sm->submodel[imod]);
#ifdef _DEBUG
    assert(mod->flag.SW3_FLOW == ON);
    assert(mod->nsys == 1);
    assert(mod->nsys_sq = 1);
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
    int i = 0, j = 0, icol = 0, ie = UNSET_INT, new_column_flag = YES, nnodes = 0, GNodeID = UNSET_INT;
    ID_LIST_ITEM *ptr;
    int nmatrix = NDONPRISM * NDONPRISM; // define using max AdH element size (Triangular Prism)
    double elem_mat[nmatrix];

#ifdef _PETSC
    int ierr = 0;
#endif
    
    // initialize matrix
    if (mod->amICoupled == NO) { // if this is a monolothic run, do this in the super load routine
#ifdef _PETSC
        // Zero the PETSc matrix
        ierr = MatZeroEntries(sm->A);
#else
        init_adh_matrix(grid->nnodes, mod->max_nsys_sq, sm->matrix, sm->diagonal);
#endif
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 2D ELEMENT CONTRIBUTIONS
    
    for (ie = 0; ie < grid->nelems2d; ie++) {
        nnodes = grid->elem2d[ie].nnodes;
        sarray_init_dbl(elem_mat, nmatrix);
        perturb_w(mod, grid->elem2d[ie].nnodes, ie, 2, elem_mat, DEBUG);
        
        // dirichlet bcs
        for(i=0; i<nnodes; i++) {
            GNodeID = grid->elem2d[ie].nodes[i];
            if(mod->bc_mask[GNodeID] == YES) {
                /* zero the row */
                for(j=0; j<nnodes; j++) {
                    elem_mat[j + i*nnodes] = 0.;
                }
                /* zero the column */
                for(j=i; j<nnodes*nnodes; j+=nnodes) {
                    elem_mat[j] = 0.;
                }
                /* put a 1. on the diagonal */
                elem_mat[i + i*nnodes] = 1.;
            }
        }
        
        
#ifdef _DEBUG
        if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
        if (DEBUG) {
            printf("2d element: %d\n",ie);
            sarray_printScreen_dbl(elem_mat, nmatrix, "elem_mat");
        }
#endif
        // assembles the element contributions
        fe_global_matrix_assemble_scalar(grid, mod->nsys, grid->elem2d[ie].nnodes, grid->elem2d[ie].nodes,
                                     mod->fmap_wvel, elem_mat,
#ifdef _PETSC
                                     sm->A, sm->Istart, sm->Iend, sm->ownership_range
#else
                                     sm->diagonal, sm->matrix
#endif
                                     );
    }
    
    //    printScreen_matrix("FINAL wvel matrix 3D element addtion", sm->diagonal, sm->matrix, grid->nnodes, mod->nsys_sq, __LINE__, __FILE__);
    //    exit(-1);
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 3D ELEMENT CONTRIBUTIONS
    
    for(icol=0; icol<grid->ncolumns; icol++)  {
        new_column_flag = YES;
        
        // loop over the 3d elements in the column
        ptr = grid->column_list[icol];
        while(ptr->next != NULL) {
            ie = ptr->id;
            nnodes = grid->elem3d[ie].nnodes;
            
#ifdef _DEBUG
            assert(icol == grid->elem3d[ie].icol);
#endif
            sarray_init_dbl(elem_mat, nmatrix);
            perturb_w(mod, grid->elem3d[ie].nnodes, ie, 3, elem_mat, DEBUG);
            
            // dirichlet bcs
            for(i=0; i<nnodes; i++) {
                GNodeID = grid->elem3d[ie].nodes[i];
                if(mod->bc_mask[GNodeID] == YES) {
                    /* zero the row */
                    for(j=0; j<nnodes; j++) {
                        elem_mat[j + i*nnodes] = 0.;
                    }
                    /* zero the column */
                    for(j=i; j<nnodes*nnodes; j+=nnodes) {
                        elem_mat[j] = 0.;
                    }
                    /* put a 1 on the diagonal */
                    elem_mat[i + i*nnodes] =1.;
                }
            }
            
            
#ifdef _DEBUG
            if (DEBUG_PICKETS) tl_check_all_pickets(__FILE__,__LINE__);
            if (DEBUG) {
                printf("3d element: %d\n",ie);
                sarray_printScreen_dbl(elem_mat, nmatrix, "elem_mat");
            }
#endif
            // assembles the element contributions
            fe_global_matrix_assemble_scalar(grid, mod->nsys, grid->elem3d[ie].nnodes, grid->elem3d[ie].nodes,mod->fmap_wvel, elem_mat,
#ifdef _PETSC
                                         sm->A, sm->Istart, sm->Iend, sm->ownership_range
#else
                                         sm->diagonal, sm->matrix
#endif
                                         );
            
            new_column_flag = NO;
            ptr = ptr->next;
        }
    }
    
    //#ifdef _DEBUG
    //    if (DEBUG) {
    //        printScreen_matrix("FINAL wvel matrix 3D element addtion", sm->diagonal, sm->matrix, grid->nnodes, mod->nsys_sq, __LINE__, __FILE__);
    //        exit(-1);
    //    }
    //#endif
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
#ifdef _PETSC
    ierr = MatAssemblyBegin(sm->A, MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(sm->A, MAT_FINAL_ASSEMBLY);
#else
    // checks the diagonal for nonexistent entries
    if (mod->amICoupled == NO) { // if this is a monolothic run, do this in the super load routine
        check_matrix_diagonal_for_nonexistant_entries(grid->nnodes, mod->max_nsys_sq, mod->nsys, sm->diagonal);
#ifdef _MESSG
        comm_update_double(sm->diagonal, 1, mod->grid->smpi);
#endif
        
#ifdef _DEBUG
        if (DEBUG_MATRIX == ON
#ifdef _MESSG
            && grid->smpi->myid == 0
#endif
            ) {
            printScreen_matrix("wvel matrix", sm->diagonal, sm->matrix, grid->nnodes, mod->max_nsys_sq, __LINE__, __FILE__);
            exit(-1);
        }
#endif
    }
    
#endif
//#ifdef _DEBUG
//    print_wvel_mat_to_file(sm,imod);
//#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Adds velocity perturbation to the 3D continuity eqn for each node in element for Newton Jacobian calculation.
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
 *  \note CJT \:: elem_rhs_bed is collected separately for monitoring. This term should be 0 if discretely
 *  consistent.  It represents flux through the opposite side of the kinematic bc (since only one
 *  side can be applied here).  It is added to the overall residual however, since it represents
 *  a mass flux.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void perturb_w(SMODEL *mod, int nodes_on_element, int ie, int dim, double *elem_mat, int DEBUG) {
    
#ifdef _DEBUG
    assert(dim == 2 || dim == 3);
#endif
    
    int i;
    double perturbation = 0.1;
    double elem_rhs_nobed[nodes_on_element];
    double elem_rhs_w_P[nodes_on_element]; // the residual resulting from a positive (P) perturbation of w
    double elem_rhs_w_M[nodes_on_element]; // the residual resulting from a negative (M) perturbation of w
    
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /* calculate perturbation  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    // get maximum velocity over the element
    // cjt :: a node perturbation is suppose to be constant, and this changes it for each element,
    // I don't think this is correct to do. Not that it matters here, since perturbation = 0.1 is big.
    int GlobalNodeID = UNSET_INT;
    double max_z_vel = 0.;
    for (i=0; i<nodes_on_element; i++) {
        if (dim == 3) {
            GlobalNodeID = mod->grid->elem3d[ie].nodes[i];
        } else {
            GlobalNodeID = mod->grid->elem2d[ie].nodes[i];
        }
        if (fabs(mod->sw->d3->vel[GlobalNodeID].z) > max_z_vel) max_z_vel = fabs(mod->sw->d3->vel[GlobalNodeID].z);
    }
    
    double epsilon = 0., epsilon_w = 0.;
    NUM_DIFF_EPSILON_GENERAL(epsilon, epsilon_w, max_z_vel, perturbation);
    
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /* calculate perturbed residuals ++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    for(i = 0; i < nodes_on_element; i++) {
        
        if (dim == 2) {
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) boundary perturbation of w-velocity +++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +w for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_w);
#endif
            fe_wvel_boundary_resid(mod,elem_rhs_w_P,elem_rhs_nobed,ie,epsilon,i,PERTURB_W,+1,DEBUG);
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) boundary perturbation of w-velocity +++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -w for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_w);
#endif
            fe_wvel_boundary_resid(mod,elem_rhs_w_M,elem_rhs_nobed,ie,epsilon,i,PERTURB_W,-1,DEBUG);
        } else if (dim == 3) {
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (+) body perturbation of w-velocity +++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing +w for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_w);
#endif
            fe_wvel_body_resid(mod,elem_rhs_w_P,ie,epsilon,i,PERTURB_W,+1,DEBUG);
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // (-) body perturbation of w-velocity +++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
            if (DEBUG) printf("\nload :: pertubing -w for %dd element %d :: perturbation: %20.10e\n",dim,ie,epsilon_w);
#endif
            fe_wvel_body_resid(mod,elem_rhs_w_M,ie,epsilon,i,PERTURB_W,-1,DEBUG);
        }
        
        /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        /* calculate residual gradient ++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        
        elem_matrix_deriv_1dof(i, nodes_on_element, elem_rhs_w_P, elem_rhs_w_M, elem_mat, epsilon_w);
        
    }
}

void print_wvel_mat_to_file(SSUPER_MODEL *sm, int imod) {
#ifdef _PETSC
    printf("In the wvel print function\n");
    SMODEL *mod = &(sm->submodel[imod]);

    int i = 0, j = 0, k = 0, l = 0, col = 0;
    FILE *fp;

    int pid, gid1, gid2, ifmap, jfmap, ierr;
    int nnodes, my_nnodes, macro_nnodes, resident_id, resident_pe, ghost_block_size, arr_index;
    
    nnodes = mod->grid->nnodes;
    my_nnodes = mod->grid->my_nnodes;
    macro_nnodes = mod->grid->macro_nnodes;
    
    int nelems = my_nnodes*macro_nnodes; // Could reduce this to just the 2 blocks in each submatrix
    int myid = sm->supersmpi->myid;
    int npes = sm->supersmpi->npes;
    PetscScalar v[nelems];
    int idxm[my_nnodes], idxn[macro_nnodes];

    // Initialize indices array
    for(pid = 0; pid < npes; pid++){
        if(myid == pid){
            //printf("Rank: %i\n",myid);
            //printf("Ownership range: %i\n",sm->ownership_range[pid+1]);
            for(i = 0; i < my_nnodes; i++){
                idxm[i] = 3*i + sm->Istart;
                //printf("idxm[%i]=%i\n",i,idxm[i]);
            }
            for(i = 0; i < macro_nnodes; i++){
                idxn[i] = 3*i;
                //printf("idxn[%i]=%i\n",i,idxn[i]);
            }
        }
        MPI_Barrier(PETSC_COMM_WORLD);
    }
    
    // Fill array with matrix values
    ierr = MatGetValues(sm->A,my_nnodes,idxm,macro_nnodes,idxn,v);

    //printf("myid is %i and npes is %i\n",myid,npes);

    // Each pe prints its array to file
    for(pid = 0; pid < npes; pid++){
        if(myid == pid){
            
            // Open file on single pe
            fp = fopen("wvel_petsc_matrix.txt","a");

            // Print header from pe 0
            if(myid == 0){
                fprintf(fp,"wvel matrix\n");
            }

            fprintf(fp,"======Rank %i======\n",myid);
           
            // Loop over (local) rows
            for(i = 0; i < my_nnodes; i++){
                ifmap = i;
                gid1 = mod->grid->node[i].gid;
                //printf("ifmap, gid = %i, %i\n",ifmap,gid1);
                
                // Loop over residential columns
                //printf("Resident columns\n");
                for(j = 0; j < my_nnodes; j++){
                    jfmap = j + sm->Istart/3;
                    gid2 = mod->grid->node[j].gid;
                    //printf("jfmap, gid = %i, %i\n",jfmap,gid2);
                    arr_index = jfmap + (i*macro_nnodes);
                    fprintf(fp,"%i %i %.8e - %i %i\n",gid1,gid2,v[arr_index],ifmap,jfmap);
                }
                
                // Loop over ghost-node columns
                //printf("Ghost columns\n");
                for(j = my_nnodes; j < nnodes; j++){
                    // Look up the global PETSc
                    resident_id = mod->grid->node[j].resident_id;
                    resident_pe = mod->grid->node[j].resident_pe;
                    gid2 = mod->grid->node[j].gid;
                    jfmap = resident_id + sm->ownership_range[resident_pe]/3;
                    //printf("jfmap, gid = %i, %i\n",jfmap,gid2);
                    //ghost_block_size = sm->ownership_range[resident_pe+1]-sm->ownership_range[resident_pe];
                    //ghost_block_size = ghost_block_size/3;
                    arr_index = jfmap + (i*macro_nnodes);
                    fprintf(fp,"%i %i %.8e - %i %i - ghost\n",gid1,gid2,v[arr_index],ifmap,jfmap);
                }
            }
            // Loop over all elements of array
            //for(i=0;i<nelems;i++){
            //    ifmap = i/ncols;
            //    jfmap = i%ncols;
            //    printf("maps: %i, %i\n",ifmap,jfmap);
            //    gid1 = mod->grid->node[ifmap].gid;
            //    gid2 = mod->grid->node[jfmap].gid;
            //    printf("gids: %i, %i\n",gid1,gid2);
            //    fprintf(fp,"%i %i %.8e\n",gid1,gid2,v[i]);
            //}
            fclose(fp);
        }
        MPI_Barrier(PETSC_COMM_WORLD);
    }
#else
    // TODO: Make work in parallel!!!
    SMODEL *mod = &(sm->submodel[imod]);

    int i = 0, j = 0, k = 0, l = 0, col = 0;
    FILE *fp = fopen("wvel_umf_matrix.txt","a");


    printf("mod nsys: %i\n",mod->nsys);
    printf("sm nsys: %i\n",sm->nsys);
    fprintf(fp,"wvel matrix:\n");
    for(i = 0; i < mod->grid->nnodes; i++){
        // Print matrix row diagonal
        for(j = 0; j < sm->nsys_sq; j++){
            fprintf(fp,"%i %i %.16e\n",i*sm->nsys+(j/sm->nsys)+1,i*sm->nsys+(j%sm->nsys)+1,sm->diagonal[i*sm->nsys_sq + j]);
        }
        for(j = 0; j < sm->matrix[i].size; j++){
            col = sm->matrix[i].index[j];
            for(k = j*sm->nsys_sq; k < (j+1)*sm->nsys_sq; k++){
                l = k - j*sm->nsys_sq;
                fprintf(fp,"%i %i %.16e %i\n",i*sm->nsys+(l/sm->nsys)+1,col*sm->nsys+(l%sm->nsys)+1,sm->matrix[i].value[k],j);
            }
        }

    }
    fclose(fp);

#endif
}
