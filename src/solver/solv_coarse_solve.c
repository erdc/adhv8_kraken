/* solves the coarse problem with an inverse matvec */

#include "global_header.h"

void solv_coarse_solve(int *bc_mask,    /* the bc mask for each dof */
                       double *result,  /* the solution vector */
                       int blk_my_nnode,    /* number of nodes I own */
                       int blk_nsys /* the number of equations being solved */
    )
{
    int i;                      /* loop counters */
    int coarse_ndof = 0;        /* the length of a coarse matrix row */
    int my_coarse_ndof = 0;     /* the length of the coarse matrix row on this processor */
    int my_block = UNSET_INT;   /* the current nodes block number */
    int my_local_block = UNSET_INT; /* the current blocks local number */
    int shift = 0;              /* keeps track of shift for N-S */

    /* determine the coarse matrix row lengths */
#ifdef _MESSG
    coarse_ndof = nblock * npes * blk_nsys;
#else
    coarse_ndof = nblock * blk_nsys;
#endif
    my_coarse_ndof = nblock * blk_nsys;

    /* initialize arrays */
    for (i = 0; i < my_coarse_ndof; i++)
        my_vr[i] = 0.0;
    for (i = 0; i < coarse_ndof; i++)
        vr[i] = 0.0;

    /* initialize the result vector */
    for (i = 0; i < coarse_ndof; i++)
        coarse_result[i] = 0.0;

    /* loop through the nodes I own and restrict the vector to the coarse grid */
    if (blk_nsys == 1) {
        for (i = 0; i < blk_my_nnode; i++) {
            /* if BC is not Dirichlet, sum entries in the same block */
            if (bc_mask[i] != YES) {
                my_block = node_block[i];
#ifdef _MESSG
                my_local_block = my_block - myid * nblock;  /* this is tied into the shift in node_order */
#else
                my_local_block = my_block;
#endif
                my_vr[my_local_block] += result[i];
            }
        }

    }
    else if (blk_nsys == 4) {
        /* "shift" moves the pointer down the vector by the size of the number
           of degrees of freedom to ensure that pressure terms are added
           to pressure terms and momentum terms are added to momentum
           terms */
        for (i = 0, shift = 0; i < blk_my_nnode; i++, shift += blk_nsys) {
            my_block = node_block[i];
#ifdef _MESSG
            /* the blocks are numbered globally, so a shift 
               must be done to get the local information */
            my_local_block = my_block - myid * nblock;  /* this is tied into the shift in node_order */
#else
            my_local_block = my_block;
#endif
            my_local_block *= blk_nsys;
            /* is BC is not Dirichlet, sum entries in the same block */
            if (bc_mask[shift + 0] != YES) {
                my_vr[my_local_block + 0] += result[shift + 0];
            }
            if (bc_mask[shift + 1] != YES) {
                my_vr[my_local_block + 1] += result[shift + 1];
            }
            if (bc_mask[shift + 2] != YES) {
                my_vr[my_local_block + 2] += result[shift + 2];
            }
            if (bc_mask[shift + 3] != YES) {
                my_vr[my_local_block + 3] += result[shift + 3];
            }

        }
    }
    else if (blk_nsys == 3) {
        /* "shift" moves the pointer down the vector by the size of the number
           of degrees of freedom to ensure that head terms are added
           to head terms and momentum terms are added to momentum
           terms */
        for (i = 0, shift = 0; i < blk_my_nnode; i++, shift += blk_nsys) {
            my_block = node_block[i];
#ifdef _MESSG
            /* the blocks are numbered globally, so a shift 
               must be done to get the local information */
            my_local_block = my_block - myid * nblock;  /* this is tied into the shift in node_order */
#else
            my_local_block = my_block;
#endif
            my_local_block *= blk_nsys;
            /* is BC is not Dirichlet, sum entries in the same block */
            if (bc_mask[shift + 0] != YES) {
                my_vr[my_local_block + 0] += result[shift + 0];
            }
            if (bc_mask[shift + 1] != YES) {
                my_vr[my_local_block + 1] += result[shift + 1];
            }
            if (bc_mask[shift + 2] != YES) {
                my_vr[my_local_block + 2] += result[shift + 2];
            }
        }
    }
    else if (blk_nsys == 2) {
        /* "shift" moves the pointer down the vector by the size of the number
           of degrees of freedom to ensure that head terms are added
           to head terms and momentum terms are added to momentum
           terms */
        for (i = 0, shift = 0; i < blk_my_nnode; i++, shift += blk_nsys) {
            my_block = node_block[i];
#ifdef _MESSG
            /* the blocks are numbered globally, so a shift 
               must be done to get the local information */
            my_local_block = my_block - myid * nblock;  /* this is tied into the shift in node_order */
#else
            my_local_block = my_block;
#endif
            my_local_block *= blk_nsys;
            /* is BC is not Dirichlet, sum entries in the same block */
            if (bc_mask[shift + 0] != YES) {
                my_vr[my_local_block + 0] += result[shift + 0];
            }
            if (bc_mask[shift + 1] != YES) {
                my_vr[my_local_block + 1] += result[shift + 1];
            }
        }
    }

    /* communicate the information on the restricted vector */
    messg_all_gather(my_vr, vr, my_coarse_ndof);

#ifdef _SuperLU_DIST
    solv_coarse_solve_sparse_parallel(coarse_ndof, vr);
#elif defined _UMFPACK
    solv_coarse_solve_sparse_serial(coarse_ndof, vr);
#else
    solv_coarse_solve_dense(coarse_ndof, vr);
#endif

    /* loop through the nodes I own; expand into the full vector */
    if (blk_nsys == 1) {
        for (i = 0; i < blk_my_nnode; i++) {
            if (bc_mask[i] != YES) {
                my_block = node_block[i];
                result[i] = coarse_result[my_block];
            }
        }
    }
    else if (blk_nsys == 4) {
        /* "shift" moves the pointer down the vector by the number of degrees
           of freedom so that pressure and momentum terms are expanded correctly */
        for (i = 0, shift = 0; i < blk_my_nnode; i++, shift += blk_nsys) {
            my_block = node_block[i] * blk_nsys;
            /* if BC is not Dirichlet, use the value */
            if (bc_mask[shift + 0] != YES) {
                result[shift + 0] = coarse_result[my_block + 0];
            }
            if (bc_mask[shift + 1] != YES) {
                result[shift + 1] = coarse_result[my_block + 1];
            }
            if (bc_mask[shift + 2] != YES) {
                result[shift + 2] = coarse_result[my_block + 2];
            }
            if (bc_mask[shift + 3] != YES) {
                result[shift + 3] = coarse_result[my_block + 3];
            }
        }
    }
    else if (blk_nsys == 3) {
        /* "shift" moves the pointer down the vector by the number of degrees
           of freedom so that head and momentum terms are expanded correctly */
        for (i = 0, shift = 0; i < blk_my_nnode; i++, shift += blk_nsys) {
            my_block = node_block[i] * blk_nsys;
            /* if BC is not Dirichlet, use the value */
            if (bc_mask[shift + 0] != YES) {
                result[shift + 0] = coarse_result[my_block + 0];
            }
            if (bc_mask[shift + 1] != YES) {
                result[shift + 1] = coarse_result[my_block + 1];
            }
            if (bc_mask[shift + 2] != YES) {
                result[shift + 2] = coarse_result[my_block + 2];
            }

        }
    }
    else if (blk_nsys == 2) {
        /* "shift" moves the pointer down the vector by the number of degrees
           of freedom so that head and momentum terms are expanded correctly */
        for (i = 0, shift = 0; i < blk_my_nnode; i++, shift += blk_nsys) {
            my_block = node_block[i] * blk_nsys;
            /* if BC is not Dirichlet, use the value */
            if (bc_mask[shift + 0] != YES) {
                result[shift + 0] = coarse_result[my_block + 0];
            }
            if (bc_mask[shift + 1] != YES) {
                result[shift + 1] = coarse_result[my_block + 1];
            }

        }
    }

}
