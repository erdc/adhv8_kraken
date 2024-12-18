#include "adh.h"
static int DEBUG = OFF;
static int DEBUG_WITH_PICKETS = OFF;

// File Prototypes
int count_vars(int score);
void load_vars(int score, int *elem_vars) ;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Counts the number of max DOFs for each element/node, stores them, and calculates fmap
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] supmod           (SUPER_MODEL *)  an AdH supermodels
 * @param[inout]       g           (SGRID *)  an ADH grid
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Column 1 -- surface water      ||  Colum 2 -- GW      ||  Column 3 -- TRANSPORT
//             0 = none                          0 = on                  0 - Ntrans
//             1 = SW 1D                         1 = off
//             2 = SW 2D
//             3 = SW 3D -- SPLIT {h,u,v}
//             4 = NS 3D
//             5 = NS 3D -- SPLIT {u,v,w}
//             6 = OF 2D
//             7 = SW 3D W
//             8 = NS 3D P
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//void smodel_super_elemVars(SMODEL_SUPER *sm, SGRID *g) {
//    
//    // ++++++++++++++++++++++++++++++++++
//    // allocate elemental variable arrays
//    // ++++++++++++++++++++++++++++++++++
//    sm->node_nvars = (int *) tl_alloc(sizeof(int), g->nnodes);
//    sm->node_vars = (int **) tl_alloc(sizeof(int *), g->nnodes);
//    
//    // Just do this in the FE-engine so no need to store
//    sm->elem_nvars = (int **) tl_alloc(sizeof(int *), 3); // 0 = 1D, 1 = 2D, 2 = 3D
//    sm->elem_nvars[0] = (int *) tl_alloc(sizeof(int), g->nelems1d);
//    sm->elem_nvars[1] = (int *) tl_alloc(sizeof(int), g->nelems2d);
//    sm->elem_nvars[2] = (int *) tl_alloc(sizeof(int), g->nelems3d);
//    
//    sm->elem_vars = (int ***) tl_alloc(sizeof(int **), 3); // 0 = 1D, 1 = 2D, 2 = 3D
//    sm->elem_vars[0] = (int **) tl_alloc(sizeof(int *), g->nelems1d);
//    sm->elem_vars[1] = (int **) tl_alloc(sizeof(int *), g->nelems2d);
//    sm->elem_vars[2] = (int **) tl_alloc(sizeof(int *), g->nelems3d);
//    
//    int i, j, k, ie, lnnodes, mat =  UNSET_INT, nelems = 0, ndim = 0, nnodes = 0, score = 0;
//    for (ndim=1; ndim<4; ndim++) {
//        nelems = UNSET_INT;
//        if (ndim == 1) {
//            nelems = g->nelems1d;
//        } else if (ndim == 2) {
//            nelems = g->nelems2d;
//        } else if (ndim == 3) {
//            nelems = g->nelems3d;
//        }
//        for (ie=0; ie<nelems; ie++) {
//            mat = UNSET_INT; lnnodes = UNSET_INT;
//            if (ndim == 1) {
//                mat = g->elem1d[ie].mat; // done in engine
//                lnnodes = 2;
//            } else if (ndim == 2) {
//                mat = g->elem2d[ie].mat;
//                lnnodes = g->elem2d[ie].nnodes;
//            } else if (ndim == 3) {
//                mat = g->elem3d[ie].mat;
//                lnnodes = g->elem3d[ie].nnodes;
//            }
//            
//            // first count independent variables on element
//            score = sm->physics_mat_code[mat]; // done in engine
//            sm->elem_nvars[ndim-1][ie] = count_vars(score);
//            
//            // find max nvars for a node over all elements connected
//            int flag[lnnodes]; sarray_init_int(flag,lnnodes);
//            for (i=0; i<lnnodes; i++) {
//                if (sm->node_nvars[g->node[i].gid] < sm->elem_nvars[ndim-1][ie]) {
//                    flag[i] = 1;
//                    sm->node_nvars[g->node[i].gid] = sm->elem_nvars[ndim-1][ie];
//                }
//            }
//            
//            // allocate element array of independent variables
//            sm->elem_vars[ndim-1][ie] = (int *) tl_alloc(sizeof(int),sm->elem_nvars[ndim-1][ie]);
//            
//            // store independent variables on element
//            load_vars(score,sm->elem_vars[ndim-1][ie]);
//            
//            // store independent variables on node
//            for (i=0; i<lnnodes; i++) {
//                if (flag[i] == 1) { // adding variables
//                    assert(sm->node_nvars[g->node[i].gid] == sm->elem_nvars[ndim-1][ie]);
//                    for (j=0; j<sm->node_nvars[g->node[i].gid]; j++) {
//                        sm->node_vars[g->node[i].gid][j] = sm->elem_vars[ndim-1][ie][j];
//                    }
//                }
//            }
//            
//            if (DEBUG) {
//                printf("ELEM_VARS\n");
//                for (i=0; i<sm->elem_nvars[ndim-1][ie]; i++) {
//                    printf("sm->elem_vars[ndim = %d][ie = %d][dof = %d]: %d\n",ndim,ie,i,sm->elem_vars[ndim-1][ie][i]);
//                }
//            }
//        
//        } // nelems
//    } //ndims
//    
//    // define fmap
//    // FOR RESIDUAL CALCULATIONS
//    // fmap_local[0:nnodes] --> for residual (ghost nodes are not a thing here)
//    // FOR MATRIX
//    // send to residential Pe, resident ID and get back the following:
//    // need sm->fmap[on resident pe of resident id] --> returns local row position of matrix, need to add
//    // need sm->node_nvars[my_nnodes:nnodes] from their residential PEs
//    // need sm->node_vars[my_nnodes:nnodes] from their residential PEs
//    // receive/send 3 integers ---> fmap, nvars, code
//    // sm->ghost[sum(all ghosts ndofs)]
//    
//    // sm->local_equation_range[min,max+1] --> starting/ending row (essentially sum of node_nvars with a start)
//    
//    // fmap --> dof_map
//    
//    
//    // need global fmap and variables for ghost nodes
//    
//    // build the local node to local equation number map
//    sm->dof_map_local = (int *) tl_alloc(sizeof(int), g->nnodes);
//    sarray_init_value_int(sm->dof_map_local,g->nnodes,UNSET_INT);
//    int count = 0;
//    for (i=0; i<g->nnodes; i++) {
//        sm->dof_map_local[i] = count;
//        count += sm->node_nvars[i];
//    }
//    
//    if (DEBUG) {
//        printf("NODE_VARS\n");
//        for (i=0; i<g->nnodes; i++) {
//            for (j=0; j<sm->node_nvars[i]; j++) {
//                printf("sm->node_vars[gid = %d][dof = %d]: %d\n",i,j,sm->node_vars[i][j]);
//            }
//        }
//        printf("\n\nFMAP\n");
//        for (i=0; i<g->nnodes; i++) {
//            printf("sm->fmap[gid = %d]: %d\n",i,sm->dof_map_local[i]);
//        }
//    }
//}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Given a physics score on an element, count the number of independent variables
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] score (int)  an ID that determines which physics are solved on this element
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// store at t0 in physics_mat[imat]
int count_vars(int score) {
    int count = 0, digit, div;
    for (div = 1; div <= score; div *= 10);

    // surface water
    div /= 10;
    digit = score/div;
    score %= div;
    printf("surface_water: %d\n",digit);
    if (digit == 1) {        count  += 2; // 1D SW {h,u_ad}
    } else if (digit == 2) { count  += 3; // 2D SW {h,u_da,v_da}
    } else if (digit == 3) { count  += 3; // 3D SW - SPLIT - HVEL -  {h,u,v}
    } else if (digit == 4) { count  += 4; // 3D NS {u,v,w,p}
    } else if (digit == 5) { count  += 3; // 3D NS - SPLIT - {u,v,w}
    } else if (digit == 6) { count  += 1; // 2D OF
    } else if (digit == 7) { count  += 1; // SW 3D - SPLIT - w
    } else if (digit == 8) { count  += 1; // 3D NS - SPLIT - p
    } else {
        tl_error("physics material error\n");
    }
    
    // ground water
    div /= 10;
    digit = score/div;
    score %= div;
    printf("ground_water: %d\n",digit);
    count += digit; // 1D/2D/3D GROUNDWATER
    
    // transport
    div /= 10;
    digit = score/div;
    score %= div;
    printf("transport: %d\n",digit);
    count += digit; // 1D/2D/3D TRANSPORT
    
    return count;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Given a physics score on an element, store theindependent variables
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] score (int)  an ID that determines which physics are solved on this element
 * @param[inout] elem_vars (int *)  stores which physics are solved on this element
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// store at t0 in physics_mat[imat]
void load_vars(int score, int *elem_vars) {
    int kk, count = 0, digit, div;
    
    for (div = 1; div <= score; div *= 10);
    
    // surface water
    div /= 10;
    digit = score/div;
    score %= div;
    if (digit == 1) {
        elem_vars[count] = PERTURB_H;   count++; // 1D SW
        elem_vars[count] = PERTURB_U;   count++; // 1D SW
    } else if (digit == 2) {
        elem_vars[count] = PERTURB_H;   count++; // 2D SW
        elem_vars[count] = PERTURB_U;   count++; // 2D SW
        elem_vars[count] = PERTURB_V;   count++; // 2D SW
    } else if (digit == 3) {
        elem_vars[count] = PERTURB_DPL; count++; // 3D SW - SPLIT
        elem_vars[count] = PERTURB_U;   count++; // 3D SW - SPLIT
        elem_vars[count] = PERTURB_V;   count++; // 3D SW - SPLIT
    } else if (digit == 4) {
        elem_vars[count] = PERTURB_U;   count++; // 3D NS
        elem_vars[count] = PERTURB_V;   count++; // 3D NS
        elem_vars[count] = PERTURB_W;   count++; // 3D NS
        elem_vars[count] = PERTURB_P;   count++; // 3D NS
    } else if (digit == 5) {
        elem_vars[count] = PERTURB_U;   count++; // 3D NS - SPLIT
        elem_vars[count] = PERTURB_V;   count++; // 3D NS - SPLIT
        elem_vars[count] = PERTURB_W;   count++; // 3D NS - SPLIT
    } else if (digit == 6) {
        elem_vars[count] = PERTURB_H;   count++; // 2D OF
    } else if (digit == 7) {
        elem_vars[count] = PERTURB_W;   count++; // 3D SW - SPLIT
    } else if (digit == 8) {
        elem_vars[count] = PERTURB_P;   count++; // 3D NS - SPLIT
    } else {
        tl_error("physics material error\n");
    }
    
    // ground water
    div /= 10;
    digit = score/div;
    score %= div;
    elem_vars[count] += PERTURB_H; count++; // 1D/2D/3D GROUNDWATER
    
    // transport
    div /= 10;
    digit = score/div;
    score %= div;
    for (kk=0; kk<digit; kk++) {
        elem_vars[count] += PERTURB_C; count++; // 1D/2D/3D TRANSPORT
    }
}
