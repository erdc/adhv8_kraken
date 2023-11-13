/* block Jacobi / additive Schwarz preconditioning */

/* this routine makes a block-diagonal approximation to the global sparse matrix. */
/* Node-to-block assignments are supplied in the array "node_block". */
/* Coefficients of the global sparse matrix are included in    */
/* block matrix iff both nodes in same block .                 */

/* The global sparse matrix format consists of two arrays.     */
/* "diagonal" is a vector containing the diagonal coefficents  */
/* with length nnodes X ndof.                                  */
/* "matrix" is a compressed row structure - an array of pointers, each pointing to    */
/* a vector containing the off-diagonal coefficents corresponding to a node.                */
/* The length of each vector is size X ndof X ndof, where size is the number of other nodes */
/* connected to the current node. */

/* The block sparse matrix format is standard compressed row storage */

#include "global_header.h"

void solv_blk_load_sparse(int *node_block,
                          SPARSE_VECT * matrix, /* global sparse matrix format */
                          double *diagonal, /* global sparse matrix format */
                          int nodes,    /* number of nodes */
                          int blk_nsys, /* ndof per node */
                          int blk_nsys_sq,  /* ndof squared  */
                          int n, int nonzeroes, long *ia, long *ja, double *a)
{
    double val;
    int nnz, j, inode, jnode, ii, jj;
    SPARSE_VECT *my_sparse_row; /* structure for one row in global sparse matrix format  */
    
    /* build the block matrix.                                                                */
    /* load those coefficients from the full matrix for which node j and i are in same block. */
    /* order of the loops was chosen so rows of the block matrix are filled sequentially      */
    
    for (nnz = 0, inode = 0; inode < nodes; inode++) {
        my_sparse_row = matrix + inode;
        for (ii = 0; ii < blk_nsys; ii++) {
            
            for (jj = 0; jj < blk_nsys; jj++) {
                val = diagonal[inode * blk_nsys_sq + ii * blk_nsys + jj];
                if (fabs(val) > SMALL) {
                    a[nnz] = val;
                    ja[nnz] = inode * blk_nsys + jj;
                    ia[nnz] = inode * blk_nsys + ii;
                    nnz++;
                }
            }
            
            for (j = 0; j < my_sparse_row->size; j++) {
                jnode = my_sparse_row->index[j];
                if (node_block[jnode] == node_block[inode]) {
                    for (jj = 0; jj < blk_nsys; jj++) {
                        val = my_sparse_row->value[j * blk_nsys_sq + ii * blk_nsys + jj];
                        if (fabs(val) > SMALL) {
                            a[nnz] = val;
                            ja[nnz] = jnode * blk_nsys + jj;
                            ia[nnz] = inode * blk_nsys + ii;
                            nnz++;
                        }
                        
                    }   /* for jj */
                }   /* if     */
            }   /* for j  */
            
        }   /* for ii */
    }   /* for i  */
    
    if (nonzeroes != nnz){
#ifdef _MESSG
        int ierr_code, myid;
        ierr_code = MPI_Comm_rank(MPI_COMM_WORLD, &myid); // cjt :: only works for 1 grid CSTORM :: MPI_COMM_WORLD, &myid);
        if(ierr_code != MPI_SUCCESS) messg_err(ierr_code);
#else
        int myid = 0;
#endif
        printf("solv_blk_load: myid %d ERROR: nonzeroes<>nnz: %d<>%d \n", myid, nonzeroes, nnz);
    }
}
