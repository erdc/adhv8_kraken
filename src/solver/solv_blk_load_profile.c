/* load the profile matrix from the sparse row format */
/* the profile solver is a direct solver in which the 
   lower triangular matrix (absent the diagonal terms)
   is stored by rows in a long vector.  The upper triangular
   matrix is stored by columns, again in a long vector.
   This routine translates from the sparse row format
   to the profile matrix format.  The sparse format consists
   of diagonal: which is a long vector containing the block 
   diagonal terms ( number of nodes X number of degrees of 
   freedom per node), and matrix: which contains the off diagonal
   block for each row (number of degrees of freedom per node squared).
 */

#include "global_header.h"

void solv_blk_load_profile(int *node_block,
                           SPARSE_VECT * matrix,    /* the matrix */
                           double *diagonal,    /* the diagonal */
                           BAND_VECT * rows,    /* the rows of the matrix */
                           BAND_VECT * cols,    /* the columns of the matrix */
                           int blk_my_nnode,    /* the number of nodes I own */
                           int blk_nsys,    /* the number of equations being solved */
                           int blk_nsys_sq  /* the number of equations being solved squared */
    )
{
    int i, j;                   /* loop counters */
    int my_block;               /* the current nodes block number */
    int other_node;             /* the other node for the matrix entry */
    int other_block;            /* the block of the other node */
    int ishift, ilong_shift, jshift, oshift;    /* the shifts in the block structure for 
                                                   multiple equations per node problems */
    SPARSE_VECT *my_sparse_row; /* the current row being worked on */


    /* initialize the row and column boundaries */
    for (i = 0; i < blk_my_nnode * blk_nsys; i++) {
        rows[i].begin = i;
        rows[i].end = i;
        cols[i].begin = i;
        cols[i].end = i;
    }


    /* load the matrix */
    for (i = 0; i < blk_my_nnode; i++) {
        /* gets the row and block of the given node */
        my_block = node_block[i];
        my_sparse_row = matrix + i;

        /* loop over the entries in the row */
        for (j = 0; j < my_sparse_row->size; j++) {
            /* if the connected node is in the block, then check the 
               row and column beginnings */
            other_node = my_sparse_row->index[j];   /* the location node j in row i */
            other_block = node_block[other_node];   /* what block is node j in */
            if (other_block == my_block) {  /* if node j is in the same block as node i then load in profile */
                /* depending on the relationship of the block numbers, set the 
                   entries in the row or the column */
                if (other_node < i) {   /* if other_node < i then this is in the lower triangular matrix
                                           and is loaded by row */
                    if (blk_nsys == 1) {
                        bv_load(rows + i, my_sparse_row->value[j], other_node);
                    }
                    else if (blk_nsys == 4) {
                        ishift = i * blk_nsys;  /* in blk global matrix, how far to skip before reaching row i */
                        jshift = j * blk_nsys_sq;   /* in the sparse matrix, how far to skip before reaching the first j entry */
                        oshift = other_node * blk_nsys; /* in blk global matrix, how far to skip before reaching column */
                        /* bv_load actually puts the entry into the profile matrix
                           blk global row ,        entry value              , blk global column    */
                        bv_load(rows + ishift + 0, my_sparse_row->value[jshift + EQN4_ID11], oshift + 0);
                        bv_load(rows + ishift + 0, my_sparse_row->value[jshift + EQN4_ID12], oshift + 1);
                        bv_load(rows + ishift + 0, my_sparse_row->value[jshift + EQN4_ID13], oshift + 2);
                        bv_load(rows + ishift + 0, my_sparse_row->value[jshift + EQN4_ID14], oshift + 3);
                        bv_load(rows + ishift + 1, my_sparse_row->value[jshift + EQN4_ID21], oshift + 0);
                        bv_load(rows + ishift + 1, my_sparse_row->value[jshift + EQN4_ID22], oshift + 1);
                        bv_load(rows + ishift + 1, my_sparse_row->value[jshift + EQN4_ID23], oshift + 2);
                        bv_load(rows + ishift + 1, my_sparse_row->value[jshift + EQN4_ID24], oshift + 3);
                        bv_load(rows + ishift + 2, my_sparse_row->value[jshift + EQN4_ID31], oshift + 0);
                        bv_load(rows + ishift + 2, my_sparse_row->value[jshift + EQN4_ID32], oshift + 1);
                        bv_load(rows + ishift + 2, my_sparse_row->value[jshift + EQN4_ID33], oshift + 2);
                        bv_load(rows + ishift + 2, my_sparse_row->value[jshift + EQN4_ID34], oshift + 3);
                        bv_load(rows + ishift + 3, my_sparse_row->value[jshift + EQN4_ID41], oshift + 0);
                        bv_load(rows + ishift + 3, my_sparse_row->value[jshift + EQN4_ID42], oshift + 1);
                        bv_load(rows + ishift + 3, my_sparse_row->value[jshift + EQN4_ID43], oshift + 2);
                        bv_load(rows + ishift + 3, my_sparse_row->value[jshift + EQN4_ID44], oshift + 3);

                    }
                    else if (blk_nsys == 3) {
                        ishift = i * blk_nsys;  /* in blk global matrix, how far to skip before reaching row i */
                        jshift = j * blk_nsys_sq;   /* in the sparse matrix, how far to skip before reaching the first j entry */
                        oshift = other_node * blk_nsys; /* in blk global matrix, how far to skip before reaching column */
                        /* bv_load actually puts the entry into the profile matrix
                           blk global row ,        entry value              , blk global column    */
                        bv_load(rows + ishift + 0, my_sparse_row->value[jshift + EQN3_ID11], oshift + 0);
                        bv_load(rows + ishift + 0, my_sparse_row->value[jshift + EQN3_ID12], oshift + 1);
                        bv_load(rows + ishift + 0, my_sparse_row->value[jshift + EQN3_ID13], oshift + 2);
                        bv_load(rows + ishift + 1, my_sparse_row->value[jshift + EQN3_ID21], oshift + 0);
                        bv_load(rows + ishift + 1, my_sparse_row->value[jshift + EQN3_ID22], oshift + 1);
                        bv_load(rows + ishift + 1, my_sparse_row->value[jshift + EQN3_ID23], oshift + 2);
                        bv_load(rows + ishift + 2, my_sparse_row->value[jshift + EQN3_ID31], oshift + 0);
                        bv_load(rows + ishift + 2, my_sparse_row->value[jshift + EQN3_ID32], oshift + 1);
                        bv_load(rows + ishift + 2, my_sparse_row->value[jshift + EQN3_ID33], oshift + 2);

                    }
                    else if (blk_nsys == 2) {
                        ishift = i * blk_nsys;  /* in blk global matrix, how far to skip before reaching row i */
                        jshift = j * blk_nsys_sq;   /* in the sparse matrix, how far to skip before reaching the first j entry */
                        oshift = other_node * blk_nsys; /* in blk global matrix, how far to skip before reaching column */
                        /* bv_load actually puts the entry into the profile matrix
                           blk global row ,        entry value              , blk global column    */
                        bv_load(rows + ishift + 0, my_sparse_row->value[jshift + EQN2_ID11], oshift + 0);
                        bv_load(rows + ishift + 0, my_sparse_row->value[jshift + EQN2_ID12], oshift + 1);
                        bv_load(rows + ishift + 1, my_sparse_row->value[jshift + EQN2_ID21], oshift + 0);
                        bv_load(rows + ishift + 1, my_sparse_row->value[jshift + EQN2_ID22], oshift + 1);

                    }

                    else
                        tl_error("solv_blk_load: blk_nsys number not available.");
                }
                else if (other_node > i) {  /* this is part of the upper triangular matrix , except for diagonal which
                                               may contain entries in both upper and lower triangular matrix */
                    if (blk_nsys == 1) {
                        bv_load(cols + other_node, my_sparse_row->value[j], i);
                    }
                    else if (blk_nsys == 4) {
                        ishift = i * blk_nsys;
                        jshift = j * blk_nsys_sq;
                        oshift = other_node * blk_nsys;
                        bv_load(cols + oshift + 0, my_sparse_row->value[jshift + EQN4_ID11], ishift + 0);
                        bv_load(cols + oshift + 1, my_sparse_row->value[jshift + EQN4_ID12], ishift + 0);
                        bv_load(cols + oshift + 2, my_sparse_row->value[jshift + EQN4_ID13], ishift + 0);
                        bv_load(cols + oshift + 3, my_sparse_row->value[jshift + EQN4_ID14], ishift + 0);
                        bv_load(cols + oshift + 0, my_sparse_row->value[jshift + EQN4_ID21], ishift + 1);
                        bv_load(cols + oshift + 1, my_sparse_row->value[jshift + EQN4_ID22], ishift + 1);
                        bv_load(cols + oshift + 2, my_sparse_row->value[jshift + EQN4_ID23], ishift + 1);
                        bv_load(cols + oshift + 3, my_sparse_row->value[jshift + EQN4_ID24], ishift + 1);
                        bv_load(cols + oshift + 0, my_sparse_row->value[jshift + EQN4_ID31], ishift + 2);
                        bv_load(cols + oshift + 1, my_sparse_row->value[jshift + EQN4_ID32], ishift + 2);
                        bv_load(cols + oshift + 2, my_sparse_row->value[jshift + EQN4_ID33], ishift + 2);
                        bv_load(cols + oshift + 3, my_sparse_row->value[jshift + EQN4_ID34], ishift + 2);
                        bv_load(cols + oshift + 0, my_sparse_row->value[jshift + EQN4_ID41], ishift + 3);
                        bv_load(cols + oshift + 1, my_sparse_row->value[jshift + EQN4_ID42], ishift + 3);
                        bv_load(cols + oshift + 2, my_sparse_row->value[jshift + EQN4_ID43], ishift + 3);
                        bv_load(cols + oshift + 3, my_sparse_row->value[jshift + EQN4_ID44], ishift + 3);

                    }
                    else if (blk_nsys == 3) {
                        ishift = i * blk_nsys;
                        jshift = j * blk_nsys_sq;
                        oshift = other_node * blk_nsys;
                        bv_load(cols + oshift + 0, my_sparse_row->value[jshift + EQN3_ID11], ishift + 0);
                        bv_load(cols + oshift + 1, my_sparse_row->value[jshift + EQN3_ID12], ishift + 0);
                        bv_load(cols + oshift + 2, my_sparse_row->value[jshift + EQN3_ID13], ishift + 0);
                        bv_load(cols + oshift + 0, my_sparse_row->value[jshift + EQN3_ID21], ishift + 1);
                        bv_load(cols + oshift + 1, my_sparse_row->value[jshift + EQN3_ID22], ishift + 1);
                        bv_load(cols + oshift + 2, my_sparse_row->value[jshift + EQN3_ID23], ishift + 1);
                        bv_load(cols + oshift + 0, my_sparse_row->value[jshift + EQN3_ID31], ishift + 2);
                        bv_load(cols + oshift + 1, my_sparse_row->value[jshift + EQN3_ID32], ishift + 2);
                        bv_load(cols + oshift + 2, my_sparse_row->value[jshift + EQN3_ID33], ishift + 2);

                    }
                    else if (blk_nsys == 2) {
                        ishift = i * blk_nsys;
                        jshift = j * blk_nsys_sq;
                        oshift = other_node * blk_nsys;
                        bv_load(cols + oshift + 0, my_sparse_row->value[jshift + EQN2_ID11], ishift + 0);
                        bv_load(cols + oshift + 1, my_sparse_row->value[jshift + EQN2_ID12], ishift + 0);
                        bv_load(cols + oshift + 0, my_sparse_row->value[jshift + EQN2_ID21], ishift + 1);
                        bv_load(cols + oshift + 1, my_sparse_row->value[jshift + EQN2_ID22], ishift + 1);

                    }

                    else
                        tl_error("solv_blk_load: blk_nsys number not available.");

                }
                else
                    tl_error("Matching values found in solv_blk_load.\n");
            }
        }

        /* loads the diagonal */
        if (blk_nsys == 1) {
            bv_load(cols + i, diagonal[i], i);
        }
        else if (blk_nsys == 4) {
            ishift = i * blk_nsys;
            ilong_shift = i * blk_nsys_sq;
            bv_load(rows + ishift + 1, diagonal[ilong_shift + EQN4_ID21], ishift + 0);
            bv_load(rows + ishift + 2, diagonal[ilong_shift + EQN4_ID31], ishift + 0);
            bv_load(rows + ishift + 2, diagonal[ilong_shift + EQN4_ID32], ishift + 1);
            bv_load(rows + ishift + 3, diagonal[ilong_shift + EQN4_ID41], ishift + 0);
            bv_load(rows + ishift + 3, diagonal[ilong_shift + EQN4_ID42], ishift + 1);
            bv_load(rows + ishift + 3, diagonal[ilong_shift + EQN4_ID43], ishift + 2);
            bv_load(cols + ishift + 0, diagonal[ilong_shift + EQN4_ID11], ishift + 0);
            bv_load(cols + ishift + 1, diagonal[ilong_shift + EQN4_ID12], ishift + 0);
            bv_load(cols + ishift + 2, diagonal[ilong_shift + EQN4_ID13], ishift + 0);
            bv_load(cols + ishift + 3, diagonal[ilong_shift + EQN4_ID14], ishift + 0);
            bv_load(cols + ishift + 1, diagonal[ilong_shift + EQN4_ID22], ishift + 1);
            bv_load(cols + ishift + 2, diagonal[ilong_shift + EQN4_ID23], ishift + 1);
            bv_load(cols + ishift + 3, diagonal[ilong_shift + EQN4_ID24], ishift + 1);
            bv_load(cols + ishift + 2, diagonal[ilong_shift + EQN4_ID33], ishift + 2);
            bv_load(cols + ishift + 3, diagonal[ilong_shift + EQN4_ID34], ishift + 2);
            bv_load(cols + ishift + 3, diagonal[ilong_shift + EQN4_ID44], ishift + 3);
        }
        else if (blk_nsys == 3) {
            ishift = i * blk_nsys;
            ilong_shift = i * blk_nsys_sq;
            bv_load(rows + ishift + 1, diagonal[ilong_shift + EQN3_ID21], ishift + 0);
            bv_load(rows + ishift + 2, diagonal[ilong_shift + EQN3_ID31], ishift + 0);
            bv_load(rows + ishift + 2, diagonal[ilong_shift + EQN3_ID32], ishift + 1);
            bv_load(cols + ishift + 0, diagonal[ilong_shift + EQN3_ID11], ishift + 0);
            bv_load(cols + ishift + 1, diagonal[ilong_shift + EQN3_ID12], ishift + 0);
            bv_load(cols + ishift + 2, diagonal[ilong_shift + EQN3_ID13], ishift + 0);
            bv_load(cols + ishift + 1, diagonal[ilong_shift + EQN3_ID22], ishift + 1);
            bv_load(cols + ishift + 2, diagonal[ilong_shift + EQN3_ID23], ishift + 1);
            bv_load(cols + ishift + 2, diagonal[ilong_shift + EQN3_ID33], ishift + 2);
        }
        else if (blk_nsys == 2) {
            ishift = i * blk_nsys;
            ilong_shift = i * blk_nsys_sq;
            bv_load(rows + ishift + 1, diagonal[ilong_shift + EQN2_ID21], ishift + 0);
            bv_load(cols + ishift + 0, diagonal[ilong_shift + EQN2_ID11], ishift + 0);
            bv_load(cols + ishift + 1, diagonal[ilong_shift + EQN2_ID12], ishift + 0);
            bv_load(cols + ishift + 1, diagonal[ilong_shift + EQN2_ID22], ishift + 1);
        }

        else
            tl_error("solv_blk_load: blk_nsys number not available.");

    }
}
