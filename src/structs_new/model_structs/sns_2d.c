#include "global_header.h"
/***********************************************************/
/***********************************************************/
/***********************************************************/

// note: need nstations here to allocate node_contrib and wind/waves
void sns_2d_alloc_init(SNS_2D **ns2d, SGRID *grid, SIO *io, SFLAGS flags) {
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
void sns_2d_init(SNS_2D *ns, int nnodes_start, int nnodes_end) {
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
// reallocs when a new node increment is added
void sns_2d_realloc_init(SNS_2D *ns, int nnodes_old, int nnodes_new) {
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
// averages two nodes for a new node
void sns_2d_node_avg(SNS_2D *ns, int node_new, int node1, int node2) {
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
// averages two nodes for a  new node
void sns_2d_renumber(SNS_2D *ns, int max_nnode, int *new_numbers, int *order_tmp, double *dtmp, SVECT2D *vtmp) {
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
// check all arrays for INF or NAN
// NOTE :: darray and iarray are just utility arrays, no need to worry about them here
void sns_2d_checkall(SNS_2D ns, int nnodes) {
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void sns_2d_print_ts(SNS_2D *ns, SIO *info, SGRID *grid, double time, int outfact, SFILE sup, char *proj_name, int it1, int it2, int **ndata, int my_nnode_max, int *my_nnode_ext, int flag, SFILE_OUTPUT file_output) {
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void sns_2d_open_output(SMODEL * mod) {
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void sns_2d_free(SNS_2D *ns, SGRID *grid, SFLAGS flags) {
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

/* Approximates elemental error for the 2d shallow water equations */
void sns_2d_calculate_elem_error(SNS_2D *ns, SGRID *grid, SMAT *mat, double dt) {
}

