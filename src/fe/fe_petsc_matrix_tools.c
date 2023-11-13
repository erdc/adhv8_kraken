#include "global_header.h"
// These routines should only be called for PETSC builds
#ifdef _PETSC
void fe_preallocate_elem_scalar(SGRID *grid, int nsys,
        int nodes_on_elem, int *GnodeID, int *fmap, Mat P);

void fe_preallocate_elem_sw2(SGRID *grid, int nsys,
        int nodes_on_elem, int *GnodeID, int *fmap, Mat P);

void fe_preallocate_elem_sw3(SGRID *grid, int nsys,
        int nodes_on_elem, int *GnodeID, int *fmap, Mat P);

void fe_fill_diagonal(Mat A, int total_nnodes_matrix,
        int total_my_nnodes, int max_nsys, int max_nsys_sq);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// This routine loops over all elements of a grid. It then loops
// over each node of the element and sets the corresponding entries
// of a MatPreallocator matrix called P. The purpose of this matrix
// is to store the location of the non-zero entries of the Jacobian
// for preallocation. The mechanism for specifying these non-zero
// entries is the MatSetValues() routine. However, the actual value
// is meaningless since we are just storing the location in a
// MatPreallocator object. Therefore, we just use arrays of zeros in
// setting the MatPreallocator P. At the end of the routine, we call
// the MatPreallocatorPreallocate() routine which uses the non-zero
// structure of the preallocator matrix P to preallocate memory for
// the Jacobian matrix A.
// 
// This code is built from a combination of the fe_transport_load
// routine to loop over the elements, and the routines in the
// fe_assemble.c file. In particular, the assemble_scalar and
// assemble_sw2 for the max_nsys=1 and max_nsys=3 cases respectively.
// As of 5/2/2022, the Navier-Stokes model will not be preallocated
// although in principle it could easily be added with a max_nsys=4
// case.
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_preallocate(SSUPER_MODEL *sm, int imod) {

    // aliases
    SMODEL *mod = &(sm->submodel[imod]);
    SGRID *grid = mod->grid;
    
    // local variable declarations
    int ie = UNSET_INT;
    int nnodes;
    PetscInt ierr;

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 3D ELEMENT CONTRIBUTIONS

    for (ie = 0; ie < grid->nelems3d; ie++) {

        nnodes = grid->elem3d[ie].nnodes;
        
        // Call assemble routine to allocate
        if(mod->max_nsys == 1){
            fe_preallocate_elem_scalar(grid, mod->max_nsys, nnodes,
                grid->elem3d[ie].nodes, mod->fmap, sm->P);
        }
        else if(mod->max_nsys == 3){
            if(mod->flag.SW3_FLOW == ON){
                fe_preallocate_elem_sw3(grid, mod->max_nsys, nnodes,
                    grid->elem3d[ie].nodes, mod->fmap, sm->P);
            }
            else{
                fe_preallocate_elem_sw2(grid, mod->max_nsys, nnodes,
                    grid->elem3d[ie].nodes, mod->fmap, sm->P);
            }
        }
        else{
            printf("Preallocation routine not set up for models with %i degrees of freedom.\n",mod->max_nsys);
            break;
        }
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 2D ELEMENT CONTRIBUTIONS

    for (ie = 0; ie < grid->nelems2d; ie++) {
        nnodes = grid->elem2d[ie].nnodes;
        if (grid->elem2d[ie].string != UNSET_INT) {
            if(mod->max_nsys == 1){
                fe_preallocate_elem_scalar(grid, mod->max_nsys, nnodes,
                    grid->elem2d[ie].nodes, mod->fmap, sm->P);
            }
            else if(mod->max_nsys == 3){
                if(mod->flag.SW3_FLOW == ON){
                    fe_preallocate_elem_sw3(grid, mod->max_nsys, nnodes,
                        grid->elem2d[ie].nodes, mod->fmap, sm->P);
                }
                else{
                    fe_preallocate_elem_sw2(grid, mod->max_nsys, nnodes,
                        grid->elem2d[ie].nodes, mod->fmap, sm->P);
                }
            }
            else{
                printf("Preallocation routine not set up for models with %i degrees of freedom.\n",mod->max_nsys);
                break;
            }
        }
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    // 1D ELEMENT CONTRIBUTIONS
   
    for (ie = 0; ie < grid->nelems1d; ie++) {

        nnodes = grid->elem1d[ie].nnodes;
        if (grid->elem1d[ie].string != UNSET_INT) {
            
            
            // assembles the element contributions
            if(mod->max_nsys == 1){
                fe_preallocate_elem_scalar(grid, mod->max_nsys, nnodes,
                    grid->elem1d[ie].nodes, mod->fmap, sm->P);
            }
            else if(mod->max_nsys == 3){
                if(mod->flag.SW3_FLOW == ON){
                    fe_preallocate_elem_sw3(grid, mod->max_nsys, nnodes,
                        grid->elem1d[ie].nodes, mod->fmap, sm->P);
                }
                else{
                    fe_preallocate_elem_sw2(grid, mod->max_nsys, nnodes,
                        grid->elem1d[ie].nodes, mod->fmap, sm->P);
                }
            }
            else{
                printf("Preallocation routine not set up for models with %i degrees of freedom.\n",mod->max_nsys);
                break;
            }
        }
    }
    
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    
    // Assemble the preallocator matrix
    // TODO: See if this Assembly step necessary for
    // MatPreallocator class. It may not be but doesn't
    // seem to hurt so leaving it for now.
    ierr = MatAssemblyBegin(sm->P, MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(sm->P, MAT_FINAL_ASSEMBLY);

    // Use the preallocator to preallocate space for A
    // Note that the second argument fills the matrix with zeros
    // so that the space will not be compressed out if it's not
    // used. For example, if we have a combination of transport
    // and sw2.
    // Note also that this automatically sets the matrix to flag
    // additional mallocs in MatSetValues with an error.
    // This option can be turned off if desired.
    ierr = MatPreallocatorPreallocate(sm->P,PETSC_TRUE,sm->A);
}

void fe_preallocate_elem_scalar(SGRID *grid, int nsys,
        int nodes_on_elem, int *GnodeID, int *fmap, Mat P)
{
    
    int nequations = 1;                  // CJT :: the number of equations per node (equal to nsys)
    int ndof = 1;                        // CJT :: ndof is the number of unknowns per node, currently always equal to nsys in AdH
    int block_size = nequations * ndof;  // CJT ::
    int ifmap, jfmap;                    // CJT :: The global matrix positions (different for monolithic runs)
    
    int i,j, index, ierr = 0;
    PetscScalar value[block_size];
    
    //printf("nodes_on_elem=%i\n",nodes_on_elem);
    // AVENGERS ASSEMBLE!
    for(i=0; i<nodes_on_elem; i++) {
        ifmap = fmap[GnodeID[i]];
#ifdef _MESSG
        if(grid->node[GnodeID[i]].resident_pe == grid->smpi->myid) {
#endif
            for(j=0; j<nodes_on_elem; j++) {
                jfmap = fmap[GnodeID[j]];
                
                index = nodes_on_elem * i + j;
                value[0] = 0; // Arbitrary value
                // Add values to matrix
                ierr = MatSetValuesLocal(P,1,&ifmap,1,&jfmap,value,ADD_VALUES);

            }
#ifdef _MESSG
        }
#endif
    }
}

void fe_preallocate_elem_sw2(SGRID *grid, int nsys,
        int nodes_on_elem, int *GnodeID, int *fmap, Mat P)
{
    int nequations = 3;                  // CJT :: the number of equations per node (equal to nsys)
    int ndof = 3;                        // CJT :: ndof is the number of unknowns per node, currently always equal to nsys in AdH
    int block_size = nequations * ndof;  // CJT :: should always be 9 here
    int ifmap, jfmap;                    // CJT :: The global matrix positions (different for monolithic runs)
    
    int i,j,index,ierr = 0;
    
    PetscScalar values[block_size];
    
    // AVENGERS ASSEMBLE!
    for (i=0; i<nodes_on_elem; i++) {
        ifmap = fmap[GnodeID[i]];
#ifdef _MESSG
        if(grid->node[GnodeID[i]].resident_pe == grid->smpi->myid){
#endif
            for (j=0; j<nodes_on_elem; j++) {
                jfmap = fmap[GnodeID[j]];
                index = nodes_on_elem * i + j;
                values[0] = 0.0;
                values[1] = 0.0;
                values[2] = 0.0;
                values[3] = 0.0;
                values[4] = 0.0;
                values[5] = 0.0;
                values[6] = 0.0;
                values[7] = 0.0;
                values[8] = 0.0;
                
                ierr = MatSetValuesBlockedLocal(P,1,&ifmap,1,&jfmap,values,ADD_VALUES);
            }
#ifdef _MESSG
        }
#endif
    }
}

// For preallocation of matrix for solving SW3 model
void fe_preallocate_elem_sw3(SGRID *grid, int nsys,
        int nodes_on_elem, int *GnodeID, int *fmap, Mat P)
{
    int nequations = 3;                  // CJT :: the number of equations per node (equal to nsys)
    int ndof = 3;                        // CJT :: ndof is the number of unknowns per node, currently always equal to nsys in AdH
    int block_size = nequations * ndof;  // CJT :: should always be 9 here
    int ifmap, jfmap;                    // CJT :: The global matrix positions (different for monolithic runs)
    
    int i,j,index,ierr = 0, iseg = UNSET_INT;
    
    PetscScalar values[block_size];

    // Get surface node IDs
    int GsurfnodeID[nodes_on_elem];
    for(i=0; i<nodes_on_elem; i++){
        iseg = find_vertical_segment(grid, GnodeID[i], grid->vertical_hash);
        GsurfnodeID[i] = grid->vertical_list[iseg]->id;
    }
    
    // Rows/columns for current element
    // AVENGERS ASSEMBLE!
    for (i=0; i<nodes_on_elem; i++) {
        ifmap = fmap[GnodeID[i]];
#ifdef _MESSG
        if(grid->node[GnodeID[i]].resident_pe == grid->smpi->myid){
#endif
            for (j=0; j<nodes_on_elem; j++) {
                jfmap = fmap[GnodeID[j]];
                index = nodes_on_elem * i + j;
                values[0] = 0.0;
                values[1] = 0.0;
                values[2] = 0.0;
                values[3] = 0.0;
                values[4] = 0.0;
                values[5] = 0.0;
                values[6] = 0.0;
                values[7] = 0.0;
                values[8] = 0.0;
                
                ierr = MatSetValuesBlockedLocal(P,1,&ifmap,1,&jfmap,values,ADD_VALUES);
            }
#ifdef _MESSG
        }
#endif
    }

    // Rows of surface element with columns of current element
    // AVENGERS ASSEMBLE!
    for (i=0; i<nodes_on_elem; i++) {
        ifmap = fmap[GsurfnodeID[i]];
#ifdef _MESSG
        if(grid->node[GsurfnodeID[i]].resident_pe == grid->smpi->myid){
#endif
            for (j=0; j<nodes_on_elem; j++) {
                jfmap = fmap[GnodeID[j]];
                index = nodes_on_elem * i + j;
                values[0] = 0.0;
                values[1] = 0.0;
                values[2] = 0.0;
                values[3] = 0.0;
                values[4] = 0.0;
                values[5] = 0.0;
                values[6] = 0.0;
                values[7] = 0.0;
                values[8] = 0.0;
                
                ierr = MatSetValuesBlockedLocal(P,1,&ifmap,1,&jfmap,values,ADD_VALUES);
            }
#ifdef _MESSG
        }
#endif
    }

    // Columns of surface element with rows of current element
    // AVENGERS ASSEMBLE!
    for (i=0; i<nodes_on_elem; i++) {
        ifmap = fmap[GnodeID[i]];
#ifdef _MESSG
        if(grid->node[GnodeID[i]].resident_pe == grid->smpi->myid){
#endif
            for (j=0; j<nodes_on_elem; j++) {
                jfmap = fmap[GsurfnodeID[j]];
                index = nodes_on_elem * i + j;
                values[0] = 0.0;
                values[1] = 0.0;
                values[2] = 0.0;
                values[3] = 0.0;
                values[4] = 0.0;
                values[5] = 0.0;
                values[6] = 0.0;
                values[7] = 0.0;
                values[8] = 0.0;
                
                ierr = MatSetValuesBlockedLocal(P,1,&ifmap,1,&jfmap,values,ADD_VALUES);
            }
#ifdef _MESSG
        }
#endif
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// This routine pads the diagonal of the PETSc matrix with 1's when the grid
// unrefines. That is, when total_nnodes_matrix > total_my_nnodes, the matrix
// will be larger than necessary for the grid and we need to make sure that there
// is at least one nonzero in every row, otherwise the solver will return an error.
//
// NOTE: Right now, this is not used. We simply destroy the matrix any time the grid
// changes at all and then create a new one. If at some point it is decided that we
// want to try to reuse "oversized" matrices, then this code could be used.
//
// ALSO NOTE: Using an explicit for-loop is probably less efficient than setting the
// whole (remaining) diagonal with a single call to MatSetValues(). I have commented
// out a sketch of how this should be done at the end of this routine.
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_fill_diagonal(Mat A, int total_nnodes_matrix,
        int total_my_nnodes, int max_nsys, int max_nsys_sq){

    int i=0, j=0, n, ierr;

    // Create block of size max_nsys x max_nsys with 1's on diagonal
    // and 0's elsewhere
    double values[max_nsys_sq];
    for(i = 0; i < max_nsys; i++){
        values[i*max_nsys] = 1.0;
        for(j = 1; j < max_nsys; j++){
            values[i*max_nsys+j] = 0.0;
        }
    }
    for(i=0;i<max_nsys_sq;i++){
        printf("values[%i]=%f\n",i,values[i]);
    }
   
    
    
    for(i = total_my_nnodes; i < total_nnodes_matrix; i++){
        printf("i=%i\n",i);
        ierr = MatSetValuesBlockedLocal(A, 1, &i, 1, &i, values, ADD_VALUES);
    }
    
   
    // TODO: Should convert to more efficiently setting the diagonal with
    // single call to MatSetValues without a for-loop; sketch of code below.
    // Need to change the values array above to have size
    // (total_nnodes_matrix-total_my_nnodes) * max_nsys_sq
    // so that it can be used one time.
    //
    // Create array of indices of block rows to fill with ID array above
    // Unset blocks are assumed to be a continguous chunk at the end of
    // each pe's range of rows with indices from total_my_nnodes to
    // (total_nnodes_matrix-1)
    //n = total_nnodes_matrix-total_my_nnodes;
    //int indices[n];
    //for(i = 0; i < n; i++){
    //    indices[i] = total_my_nnodes+i;
    //    printf("indices[%i]=%i\n",i,indices[i]);
    //}

    //ierr = MatSetValuesBlockedLocal(A, n, indices, n, indices, values, ADD_VALUES);
}
#endif // End the petsc if-def clause
