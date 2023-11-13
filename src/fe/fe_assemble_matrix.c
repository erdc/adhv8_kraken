/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  fe_assemble_matrix.c This file collects functions that assemble the globa Jacobi matrix. */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Assembles a global jacobi matrix from local matrices for a scalar equation
 *  \author    Jesus
 *  \author    Charlie Berger, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] matrix stores off-diagonal blocks of the Jacobi Matrix
 *  @param[in,out] diagonal stores the diagonal blocks of the Jacobi Matrix
 *  @param[in] grid the sw2 grid
 *  @param[in] nsys the number of equations per node to be simultaneously solved, should be 1 here
 *  @param[in] nnodes_on_elem the number of nodes on element
 *  @param[in] GnodeID the global IDs of the local elemental nodes
 *  @param[in] elem_mat stores the local elemental matrix Jacobian
 *
 * \details The scalar equation Jacobi Matrix is an nnode X nnode matrix with elements
 * \f$ J_{i,j} = \deriv{R^{i}}{v^{j}} \f$ \n
 *  where "v" a generic independent variable of the local residual.
 *
 *  A diagonal block here is one double, i.e.  \n
 *  diagonal[inode] = \f$ \deriv{R^{inode}}{v^{inode}} \f$ \n
 *
 *  The local elemental matrix assembled looks like (for triangular element): \n
 *  elem_mat[0] = {irow = phi_0 * R_elem, icol = inode = 0} \n
 *  elem_mat[1] = {irow = phi_0 * R_elem, icol = inode = 1} \n
 *  elem_mat[2] = {irow = phi_0 * R_elem, icol = inode = 2} \n
 *  elem_mat[3] = {irow = phi_1 * R_c_elem, icol = inode = 0} \n
 *  elem_mat[4] = {irow = phi_1 * R_c_elem, icol = inode = 1} \n
 *  elem_mat[5] = {irow = phi_1 * R_c_elem, icol = inode = 2} \n
 *  elem_mat[6] = {irow = phi_2 * R_c_elem, icol = inode = 0} \n
 *  elem_mat[7] = {irow = phi_2 * R_c_elem, icol = inode = 1} \n
 *  elem_mat[8] = {irow = phi_2 * R_c_elem, icol = inode = 2} \n
 *  \n \n
 *  The off-diagonal blocks are sparse, and are stored as a SPARSE_VECT.
 *  \note CJT \:: added PETSC interface :: MatSetValues(Mat A,PetscInt m,PetscInt *im,PetscInt n,PetscInt *in,PetscScalar *values,INSERT VALUES)
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_global_matrix_assemble_scalar(SGRID *grid, int nsys, int nodes_on_elem, int *GnodeID, int *fmap, double *elem_mat//, double *diagonal, SPARSE_VECT *matrix
#ifdef _PETSC
        , Mat A, int Istart, int Iend, const int *ownership_range
#else
        , double *diagonal, SPARSE_VECT *matrix
#endif
        ) {
    
#ifdef _DEBUG
    assert(nsys == 1);
#endif
    
    int nequations = 1;                  // CJT :: the number of equations per node (equal to nsys)
    int ndof = 1;                        // CJT :: ndof is the number of unknowns per node, currently always equal to nsys in AdH
    int block_size = nequations * ndof;  // CJT ::
    int ifmap, jfmap;                    // CJT :: The global matrix positions (different for monolithic runs)
    
    int i,j, index, ierr = 0;
#ifdef _PETSC
    // Hard code values array to be max block size for
    // models with multiple equations (e.g. transport+sw2).
    // Reason for this is that matrices are allocated based
    // on max_nsys which may be larger than 1 so we need to
    // make sure there are non-zeros on diagonal.
    // NOTE: needs to be updated to work with N-S which has
    // more degrees of freedom.
    // NOTE: May be more efficient to code this to be an input to
    // this routine so an array isn't allocated each time, would also
    // allow us to allocate an array of the correct size (e.g. array of
    // size 1 for pure transport with no other equations being solved)
    PetscScalar value[9];

    // Set "dummy" entries of values
    // Need 1's on diagonal, 0's elsewhere
    value[1] = 0;
    value[2] = 0;
    value[3] = 0;
    value[4] = 0;
    value[5] = 0;
    value[6] = 0;
    value[7] = 0;
    value[8] = 0;
#else
    double value[block_size];
#endif
    
    // AVENGERS ASSEMBLE!
    for(i=0; i<nodes_on_elem; i++) {
        ifmap = fmap[GnodeID[i]];
#ifdef _MESSG
        if(grid->node[GnodeID[i]].resident_pe == grid->smpi->myid) {
#endif
            for(j=0; j<nodes_on_elem; j++) {
                jfmap = fmap[GnodeID[j]];

                index = nodes_on_elem * i + j;
                value[0] = elem_mat[index];
                if(fabs(value[0]) > SMALL) { // cjt :: this is questionable to me
#ifdef _PETSC
                    // We need 1's on the main diagonal
                    // when max_nsys > nsys
                    if(ifmap == jfmap){
                        value[4] = 1;
                        value[8] = 1;
                    }
                    // Not on main diagonal
                    // TODO: works with fmap?
                    else{
                        value[4] = 0;
                        value[8] = 0;
                    }

                    // Add values to matrix
                    ierr = MatSetValuesBlockedLocal(A,1,&ifmap,1,&jfmap,value,ADD_VALUES);

#else
                    if(ifmap != jfmap)  {
                        spv_load(matrix+ifmap, jfmap, value, nsys, block_size, block_size);
                    } else {
                        diagonal[ifmap] += value[0];
                    }
#endif
                }
            }
#ifdef _MESSG
        }
#endif
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Assembles the global SW 2D jacobi matrix from local matrices
 *  \author    Charlie Berger, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] matrix stores off-diagonal blocks of the Jacobi Matrix
 *  @param[in,out] diagonal stores the diagonal blocks of the Jacobi Matrix
 *  @param[in] grid the sw2 grid
 *  @param[in] nsys the number of equations per node to be simultaneously solved, should be 3 here
 *  @param[in] nnodes_on_elem the number of nodes on element
 *  @param[in] GnodeID the global IDs of the local elemental nodes
 *  @param[in] elem_mat_u stores the Jacobian, local elemental matrix for the x-velocity perturbation
 *  @param[in] elem_mat_v stores the Jacobian, local elemental matrix for the y-velocity perturbation
 *  @param[in] elem_mat_h stores the Jacobian, local elemental matrix for the depth perturbation
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void fe_global_matrix_assemble_sw2(SGRID *grid, int nsys, int nodes_on_elem, int *GnodeID, int *fmap, DOF_3 *elem_mat_u, DOF_3 *elem_mat_v, DOF_3 *elem_mat_h,
#ifdef _PETSC
        Mat A, int Istart, int Iend, const int *ownership_range
#else
        double *diagonal, SPARSE_VECT *matrix
#endif
        ) {
    
#ifdef _DEBUG
    assert(nsys == 3);
#endif
    
    int nequations = 3;                  // CJT :: the number of equations per node (equal to nsys)
    int ndof = 3;                        // CJT :: ndof is the number of unknowns per node, currently always equal to nsys in AdH
    int block_size = nequations * ndof;  // CJT :: should always be 9 here
    int ifmap, jfmap;                    // CJT :: The global matrix positions (different for monolithic runs)
    
    int i,j,k,index,iloc,ierr = 0;
    
#ifdef _PETSC
    PetscScalar values[block_size];
#else
    double values[block_size];
#endif
    
    // AVENGERS ASSEMBLE!
    for (i=0; i<nodes_on_elem; i++) {
        ifmap = fmap[GnodeID[i]];
#ifdef _MESSG
        if(grid->node[GnodeID[i]].resident_pe == grid->smpi->myid){
#endif
            for (j=0; j<nodes_on_elem; j++) {
                jfmap = fmap[GnodeID[j]];
                index = nodes_on_elem * i + j;
                values[0] = elem_mat_u[index].x_eq;
                values[1] = elem_mat_v[index].x_eq;
                values[2] = elem_mat_h[index].x_eq;
                values[3] = elem_mat_u[index].y_eq;
                values[4] = elem_mat_v[index].y_eq;
                values[5] = elem_mat_h[index].y_eq;
                values[6] = elem_mat_u[index].c_eq;
                values[7] = elem_mat_v[index].c_eq;
                values[8] = elem_mat_h[index].c_eq;
#ifdef _PETSC
                ierr = MatSetValuesBlockedLocal(A,1,&ifmap,1,&jfmap,values,ADD_VALUES);
#else
                if (jfmap != ifmap) {
                    spv_load(matrix+ifmap, jfmap, values, nequations, block_size, block_size);
                } else {
                    iloc = ifmap * block_size;
                    for (k=0; k<block_size; k++) {diagonal[iloc + k] += values[k];}
                }
#endif
            }
#ifdef _MESSG
        }
#endif
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Assembles the global SW 3D jacobi matrix from local matrices
 *  \author    Charlie Berger, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] matrix stores off-diagonal blocks of the Jacobi Matrix
 *  @param[in,out] diagonal stores the diagonal blocks of the Jacobi Matrix
 *  @param[in] dim the dimensionality of the local element to be added
 *  @param[in] grid the sw2 grid
 *  @param[in] nsys the number of equations per node to be simultaneously solved, should be 3 here
 *  @param[in] nnodes_on_elem the number of nodes on element
 *  @param[in] ie the element id
 *  @param[in] elem_mat_u stores the Jacobian, local elemental matrix for the x-velocity perturbation
 *  @param[in] elem_mat_v stores the Jacobian, local elemental matrix for the y-velocity perturbation
 *  @param[in] elem_mat_dpl stores the Jacobian, local elemental matrix for the displacement perturbation
 *
 * \details The Jacobi Matrix is arranged as:
 * \f{eqnarray*}{ \AdHMatrix \f} \n
 *  where the diagonal blocks are enclosed by () and the off-diagonal blocks are enclosed with [].
 *  Here, the number of dofs (u,v,d) in the ith residual = the number of equations (mx,my,dc).
 *
 *  In AdH, the diagonal blocks are stored as an array going across rows first, i.e. \n
 *  diagonal[inode * nequations * ndof + 0] = \f$ \pderiv{R_{mx}^{inode}}{u^{inode}} \f$ \n
 *  diagonal[inode * nequations * ndof + 1] = \f$ \pderiv{R_{mx}^{inode}}{v^{inode}} \f$ \n
 *  diagonal[inode * nequations * ndof + 2] = \f$ \pderiv{R_{mx}^{inode}}{d^{inode}} \f$ \n
 *  diagonal[inode * nequations * ndof + 3] = \f$ \pderiv{R_{my}^{inode}}{u^{inode}} \f$ \n
 *  diagonal[inode * nequations * ndof + 4] = \f$ \pderiv{R_{my}^{inode}}{v^{inode}} \f$ \n
 *  diagonal[inode * nequations * ndof + 5] = \f$ \pderiv{R_{my}^{inode}}{d^{inode}} \f$ \n
 *  diagonal[inode * nequations * ndof + 6] = \f$ \pderiv{R_{c}^{inode}}{u^{inode}} \f$  \n
 *  diagonal[inode * nequations * ndof + 7] = \f$ \pderiv{R_{c}^{inode}}{v^{inode}} \f$  \n
 *  diagonal[inode * nequations * ndof + 8] = \f$ \pderiv{R_{c}^{inode}}{d^{inode}} \f$  \n
 *  \n
 *  The array "values" temporarily stores variables used to assemble the local matrices into the global.
 *  For a particular test function i and node j you get the following vector values
 *  value[0] = \f$ \pderiv{R_{mx}^{i}}{u^{j}} \f$ \n
 *  value[1] = \f$ \pderiv{R_{mx}^{i}}{v^{j}} \f$ \n
 *  value[2] = \f$ \pderiv{R_{mx}^{i}}{d^{j}} \f$ \n
 *  value[3] = \f$ \pderiv{R_{my}^{i}}{u^{j}} \f$ \n
 *  value[4] = \f$ \pderiv{R_{my}^{i}}{v^{j}} \f$ \n
 *  value[5] = \f$ \pderiv{R_{my}^{i}}{d^{j}} \f$ \n
 *  value[6] = \f$ \pderiv{R_{dc}^{i}}{u^{j}} \f$ \n
 *  value[7] = \f$ \pderiv{R_{dc}^{i}}{v^{j}} \f$ \n
 *  value[8] = \f$ \pderiv{R_{dc}^{i}}{d^{j}} \f$ \n
 *  \n
 *  These values are taken from the local elemental matrix, stored as: \n
 *  elem_mat_u[0].c_eq = {irow = phi_0 * R_c_elem, icol = inode = 0} \n
 *  elem_mat_u[1].c_eq = {irow = phi_0 * R_c_elem, icol = inode = 1} \n
 *  elem_mat_u[2].c_eq = {irow = phi_0 * R_c_elem, icol = inode = 2} \n
 *  elem_mat_u[3].c_eq = {irow = phi_0 * R_c_elem, icol = inode = 3} \n
 *  elem_mat_u[4].c_eq = {irow = phi_1 * R_c_elem, icol = inode = 0} \n
 *  elem_mat_u[5].c_eq = {irow = phi_1 * R_c_elem, icol = inode = 1} \n
 *  elem_mat_u[6].c_eq = {irow = phi_1 * R_c_elem, icol = inode = 2} \n
 *  elem_mat_u[7].c_eq = {irow = phi_1 * R_c_elem, icol = inode = 3} \n
 *  elem_mat_u[8].c_eq = {irow = phi_2 * R_c_elem, icol = inode = 0} \n
 *  elem_mat_u[9].c_eq = {irow = phi_2 * R_c_elem, icol = inode = 1} \n
 *  elem_mat_u[10].c_eq = {irow = phi_2 * R_c_elem, icol = inode = 2} \n
 *  elem_mat_u[11].c_eq = {irow = phi_2 * R_c_elem, icol = inode = 3} \n
 *  elem_mat_u[12].c_eq = {irow = phi_3 * R_c_elem, icol = inode = 0} \n
 *  elem_mat_u[13].c_eq = {irow = phi_3 * R_c_elem, icol = inode = 1} \n
 *  elem_mat_u[14].c_eq = {irow = phi_3 * R_c_elem, icol = inode = 2} \n
 *  elem_mat_u[15].c_eq = {irow = phi_3 * R_c_elem, icol = inode = 3} \n
 *  \n \n
 *  The off-diagonal blocks are sparse, and are stored as a SPARSE_VECT.
 *  \note CJT \:: the depth-averaged continuity residual Jacobian only contributes to surface nodes
 *  \note CJT \:: spv_load(matrix row - equation (phi_i * R) , dense sub-block column id (node), ...)
 *  \note CJT \:: added PETSC interface :: MatSetValues(Mat A,PetscInt m,PetscInt *im,PetscInt n,PetscInt *in,PetscScalar *values,INSERT VALUES)
 *
 *  \note CJT \::  CURRENTLY, THE WAY DISPLACEMENT PERTURBATIONS IN THE JACOBIAN ARE HANDLED IS QUESTIONABLE.
 *          IF THE DISPLACMENTS WERE CONSTANT DOWN A COLUMN, THE CURRENT ASSEMBLY WOULD BE FINE, BUT THEY ARE NOT.
 *          SUB-SURFACE DISPLACEMENTS ARE PROPORTIONAL TO SURFACE, BUT NOT IDENTICAL.  THEY ARE DEPENDENT THOUGH, SO IDK.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_global_matrix_assemble_sw3(int dim, SGRID *grid, int nsys, int nodes_on_elem, int ie, int *fmap, DOF_3 *elem_mat_u, DOF_3 *elem_mat_v, DOF_3 *elem_mat_dpl,
#ifdef _PETSC
        Mat A
#else
        double *diagonal, SPARSE_VECT *matrix
#endif
        ) {
    
#ifdef _DEBUG
    assert(nsys == 3);
#endif
  
    int nequations = 3;                  // CJT :: the number of equations per node (equal to nsys)
    int ndof = 3;                        // CJT :: ndof is the number of unknowns per node, currently always equal to nsys in AdH
    int block_size = nequations * ndof;  // CJT :: should always be 9 here
    int ifmap, jfmap;                    // CJT :: The global matrix positions (different for monolithic runs)
    
    int i, j, k, ierr = 0, index, iseg = UNSET_INT;
#ifdef _PETSC
    PetscScalar values[block_size];
    //PetscInt *block_row_start, *block_col_start;
#else
    double values[block_size];
    int iloc;
#endif
    
    // get global node and surface node IDs
    int GsurfnodeID[nodes_on_elem], GnodeID[nodes_on_elem];
    for(i=0; i<nodes_on_elem; i++) {
        if (dim == 3) {
            GnodeID[i] = grid->elem3d[ie].nodes[i];
        } else {
            GnodeID[i] = grid->elem2d[ie].nodes[i];
        }
        iseg = find_vertical_segment(grid, GnodeID[i], grid->vertical_hash);
        GsurfnodeID[i] = grid->vertical_list[iseg]->id;
    }
    
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                                DA-CONT DERIVATIVE W.R.T {U,V}
     *--------------------------------------------------------------------------------------------
     * This loops adds the depth-averaged continuity u,v derivatives to the Jacobian Matrix. \n
     * Residuals are added up a column, f$ R_I = R_i + R_{i+1} + R_{i+2} + ... \n
     * This sections only adds to surface residuals
     *********************************************************************************************/
    for(i=0; i<nodes_on_elem; i++) {
        ifmap = fmap[GsurfnodeID[i]]; // surface node!
#ifdef _MESSG
        if(grid->node[GnodeID[i]].resident_pe == grid->smpi->myid) {
#endif
            for(j=0; j<nodes_on_elem; j++) {
                jfmap = fmap[GnodeID[j]];
                index = nodes_on_elem * i + j;
                values[0] = 0.;                      // dRmx[inode]/du[inode]
                values[1] = 0.;                      // dRmx[inode]/dv[inode]
                values[2] = 0.;                      // dRmx[inode]/dh[inode]
                values[3] = 0.;                      // dRmy[inode]/du[inode]
                values[4] = 0.;                      // dRmy[inode]/dv[inode]
                values[5] = 0.;                      // dRmy[inode]/dh[inode]
                values[6] = elem_mat_u[index].c_eq;  // dRdc[inode]/du[inode]
                values[7] = elem_mat_v[index].c_eq;  // dRdc[inode]/dv[inode]
                values[8] = 0.;                      // dRdc[inode]/dh[inode]
#ifdef _PETSC
                ierr = MatSetValuesBlockedLocal(A,1,&(ifmap),1,&(jfmap),values,ADD_VALUES); 
#else
                if(jfmap != ifmap)  {
                    spv_load(matrix+ifmap, jfmap, values, nequations, block_size, block_size);
                } else {
                    iloc = ifmap * block_size;
                    for (k=0; k<block_size; k++) {diagonal[iloc + k] += values[k];}
                }
#endif
            }
#ifdef _MESSG
        }
#endif
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                            DA-CONT DERIVATIVE W.R.T DISPLACEMENT
     *--------------------------------------------------------------------------------------------
     * This loops adds the depth-averaged continuity displacement derivatives to the Jacobian Matrix. \n
     * Residuals are added up a column, f$ R_I = R_i + R_{i+1} + R_{i+2} + ... \n
     * This sections only adds to surface residuals
     *********************************************************************************************/
    for(i=0; i<nodes_on_elem; i++) {
        ifmap = fmap[GsurfnodeID[i]]; // Surface node!
#ifdef _MESSG
        if(grid->node[GnodeID[i]].resident_pe == grid->smpi->myid) {
#endif
            for(j=0; j<nodes_on_elem; j++) {
                jfmap = fmap[GsurfnodeID[j]]; // Surface node!
                index = nodes_on_elem * i + j;
                values[0] = 0.;
                values[1] = 0.;
                values[2] = 0.;
                values[3] = 0.;
                values[4] = 0.;
                values[5] = 0.;
                values[6] = 0.;
                values[7] = 0.;
                values[8] = elem_mat_dpl[index].c_eq; // dRdc[inode]/dh[inode]
#ifdef _PETSC
                ierr = MatSetValuesBlockedLocal(A,1,&(ifmap),1,&(jfmap),values,ADD_VALUES); 
#else
                if(jfmap != ifmap) {
                    spv_load(matrix+ifmap, jfmap, values, nequations, block_size, block_size);
                } else {
                    iloc = ifmap * block_size;
                    for (k=0; k<block_size; k++) {diagonal[iloc + k] += values[k];}
                }
#endif
            }
#ifdef _MESSG
        }
#endif
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                            MOMENTUM DERIVATIVES W.R.T {U,V}
     *--------------------------------------------------------------------------------------------
     * This loops adds the momentum velocity derivatives to the Jacobian Matrix. \n
     *********************************************************************************************/
    for(i=0; i<nodes_on_elem; i++) {
        ifmap = fmap[GnodeID[i]];
#ifdef _MESSG
        if(grid->node[GnodeID[i]].resident_pe == grid->smpi->myid) {
#endif
            for(j=0; j<nodes_on_elem; j++) {
                jfmap = fmap[GnodeID[j]];
                index = nodes_on_elem * i + j;
                values[0] = elem_mat_u[index].x_eq; // dRmx[inode]/du[inode]
                values[1] = elem_mat_v[index].x_eq; // dRmx[inode]/dv[inode]
                values[2] = 0.;
                values[3] = elem_mat_u[index].y_eq; // dRmy[inode]/du[inode]
                values[4] = elem_mat_v[index].y_eq; // dRmy[inode]/dv[inode]
                values[5] = 0.;
                values[6] = 0.;
                values[7] = 0.;
                values[8] = 0.;
#ifdef _PETSC
                ierr = MatSetValuesBlockedLocal(A,1,&(ifmap),1,&(jfmap),values,ADD_VALUES); 
#else
                if(jfmap != ifmap) { // off-diagonal block
                    spv_load(matrix+ifmap, jfmap, values, nequations, block_size, block_size);
                }
                else {
                    iloc = ifmap * block_size;
                    for (k=0; k<block_size; k++) {diagonal[iloc + k] += values[k];}
                }
#endif
            }
#ifdef _MESSG
        }
#endif
    }
    
    /*!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     *                           MOMENTUM DERIVATIVES W.R.T DISPLACEMENT
     *--------------------------------------------------------------------------------------------
     * This loops adds the momentum displacement derivatives to the Jacobian Matrix. \n
     *********************************************************************************************/
    for(i=0; i<nodes_on_elem; i++) {
        ifmap = fmap[GnodeID[i]];
#ifdef _MESSG
        if(grid->node[GnodeID[i]].resident_pe == grid->smpi->myid) {
#endif
            for(j=0; j<nodes_on_elem; j++) { // loop over node j dof
                jfmap = fmap[GsurfnodeID[j]]; // SURFACE NODE!
                index = nodes_on_elem * i + j;
                values[0] = 0.;
                values[1] = 0.;
                values[2] = elem_mat_dpl[index].x_eq; // Rmx[inode]/dh[inode]
                values[3] = 0.;
                values[4] = 0.;
                values[5] = elem_mat_dpl[index].y_eq; // Rmy[inode]/dh[inode]
                values[6] = 0.;
                values[7] = 0.;
                values[8] = 0.;
#ifdef _PETSC
                ierr = MatSetValuesBlockedLocal(A,1,&(ifmap),1,&(jfmap),values,ADD_VALUES); 
#else
                if(jfmap != ifmap) {
                    spv_load(matrix+ifmap, jfmap, values, nequations, block_size, block_size);
                }
                else {
                    iloc = ifmap * block_size;
                    for (k=0; k<block_size; k++) {diagonal[iloc + k] += values[k];}
                }
#endif
            }
        }
#ifdef _MESSG
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Assembles the global 3D Navier Stokes jacobi matrix from local matrices
 *  \author    Charlie Berger, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  @param[in,out] matrix stores off-diagonal blocks of the Jacobi Matrix
 *  @param[in,out] diagonal stores the diagonal blocks of the Jacobi Matrix
 *  @param[in] dim the dimensionality of the local element to be added
 *  @param[in] grid the sw2 grid
 *  @param[in] nsys the number of equations per node to be simultaneously solved, should be 3 here
 *  @param[in] nnodes_on_elem the number of nodes on element
 *  @param[in] ie the element id
 *  @param[in] elem_mat_u stores the Jacobian, local elemental matrix for the x-velocity perturbation
 *  @param[in] elem_mat_v stores the Jacobian, local elemental matrix for the y-velocity perturbation
 *  @param[in] elem_mat_dpl stores the Jacobian, local elemental matrix for the displacement perturbation
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void fe_global_matrix_assemble_ns3(int dim, SGRID *grid, int nsys, int nodes_on_elem, int ie, int *fmap, DOF_4 *elem_mat_u, DOF_4 *elem_mat_v, DOF_4 *elem_mat_w, DOF_4 *elem_mat_p, double *diagonal, SPARSE_VECT *matrix) {
    
#ifdef _DEBUG
    assert(nsys == 4);
#endif
    
    int nequations = 4;                  // CJT :: the number of equations per node (equal to nsys)
    int ndof = 4;                        // CJT :: ndof is the number of unknowns per node, currently always equal to nsys in AdH
    int block_size = nequations * ndof;  // CJT :: should always be 9 here
    int ifmap, jfmap;                    // CJT :: The global matrix positions (different for monolithic runs)
    
    int i, j, k, ierr = 0, index;
#ifdef PETSC
    PetscScalar values[block_size];
    PetscInt *block_row_start, *block_col_start;
#else
    double values[block_size];
    int iloc;
#endif
    
    
    // get global node and surface node IDs
    int GnodeID[nodes_on_elem];
    for(i=0; i<nodes_on_elem; i++) {
        if (dim == 3) {
            GnodeID[i] = grid->elem3d[ie].nodes[i];
        } else {
            GnodeID[i] = grid->elem2d[ie].nodes[i];
        }
    }

    for(i=0; i<nodes_on_elem; i++) {
        ifmap = fmap[GnodeID[i]];
#ifdef _MESSG
        if(grid->node[GnodeID[i]].resident_pe == grid->smpi->myid) {
#endif
            for(j=0; j<nodes_on_elem; j++) {
                jfmap = fmap[GnodeID[j]];
                index = nodes_on_elem * i + j;
                values[0]  = elem_mat_u[index].x_eq; // dRmx[inode]/du[inode]
                values[1]  = elem_mat_v[index].x_eq; // dRmx[inode]/dv[inode]
                values[2]  = elem_mat_w[index].x_eq; // dRmx[inode]/dw[inode]
                values[3]  = elem_mat_p[index].x_eq; // dRmx[inode]/dp[inode]
                
                values[4]  = elem_mat_u[index].y_eq; // dRmy[inode]/du[inode]
                values[5]  = elem_mat_v[index].y_eq; // dRmy[inode]/dv[inode]
                values[6]  = elem_mat_w[index].y_eq; // dRmy[inode]/dw[inode]
                values[7]  = elem_mat_p[index].y_eq; // dRmy[inode]/dp[inode]
                
                values[8]  = elem_mat_u[index].z_eq; // dRmz[inode]/du[inode]
                values[9]  = elem_mat_v[index].z_eq; // dRmz[inode]/dv[inode]
                values[10] = elem_mat_w[index].z_eq; // dRmz[inode]/dw[inode]
                values[11] = elem_mat_p[index].z_eq; // dRmz[inode]/dp[inode]
                
                values[12] = elem_mat_u[index].c_eq; // dC[inode]/du[inode]
                values[13] = elem_mat_v[index].c_eq; // dC[inode]/dv[inode]
                values[14] = elem_mat_w[index].c_eq; // dC[inode]/dw[inode]
                values[15] = elem_mat_p[index].c_eq; // dC[inode]/dp[inode]
#ifdef PETSC
                *block_row_start = ifmap * ndof;
                *block_col_start = jfmap * ndof;
                ierr = MatSetValuesBlocked(A,1,block_row_start,1,block_col_start,values,ADD_VALUES); 
#else
                if(jfmap != ifmap) { // off-diagonal block
                    spv_load(matrix+ifmap, jfmap, values, nequations, block_size, block_size);
                }
                else {
                    iloc = ifmap * block_size;
                    for (k=0; k<block_size; k++) {diagonal[iloc + k] += values[k];}
                }
#endif
            }
#ifdef _MESSG
        }
#endif
    }
    
}












