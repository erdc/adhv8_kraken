#include "global_header.h"

//**************************************************************************************************************//
//**************************************************************************************************************//

// noptions :: the number of integer options associated with the elem property
// lowerORhigher :: choice to choose either lower (-1) or higher (0) option at ambiguous corner node

void elem3d_to_node_int(SGRID *grid, int *nodal_var, int *elem3d_var, int noptions, int lowerORhigher) {

    int ie=0, lnode=0, ioption=0;
    int node[MAX_NNODES_ON_ELEM3D];

    int lower = -1;
    int higher = +1;

    assert(nodal_var);
    assert(elem3d_var);
    
    if (lowerORhigher != lower && lowerORhigher != higher) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> input lower or higher option is not valid\n");
    }

    for (ie=0; ie < grid->nelems3d; ie++) {
        for (lnode=0; lnode<grid->elem3d[ie].nnodes; lnode++) {
            node[lnode] = grid->elem3d[ie].nodes[lnode];
            nodal_var[node[lnode]] = UNSET_INT;
        }
    }

    for (ie=0; ie < grid->nelems3d; ie++) {                 // loop over elements
        for (ioption=0; ioption<noptions; ioption++) {        // loop over integer options
            if (elem3d_var[ie] == ioption) {
                for (lnode=0; lnode<grid->elem3d[ie].nnodes; lnode++) {     // loop over local nodes
                    node[lnode] = grid->elem3d[ie].nodes[lnode];
                    if (nodal_var[ie] == UNSET_INT) { 
                        nodal_var[node[lnode]] = ioption;
                    } else if ( (nodal_var[ie] < ioption) && (lowerORhigher == higher)) {
                        nodal_var[node[lnode]] = ioption;
                    } else if ( (nodal_var[ie] > ioption) && (lowerORhigher == lower)) {
                        nodal_var[node[lnode]] = ioption;
                    } else {
                        nodal_var[node[lnode]] = ioption;
                    }
                }
            }
        }
    }
}

//**************************************************************************************************************//
//**************************************************************************************************************//

void elem3d_to_node_double(SGRID *grid, STR_VALUE *string, double *node_var, double *elem3d_var) {
    
    int ie=0, inode=0, gnode=0;
    double weight = 0.;
    double elem_sfield[MAX_NNODES_ON_ELEM3D];
    
    assert(node_var);
    assert(elem3d_var);
    
    // allocate array to hold weight sum
    double *weight_sum = (double *) tl_alloc(sizeof(double), grid->nnodes);
    
    sarray_init_dbl(node_var, grid->nnodes);
    sarray_init_dbl(weight_sum, grid->nnodes);
    for(ie = 0; ie < grid->nelems3d; ie++) {
        if (grid->elem3d[ie].string != UNSET_INT) { // the 3d element is defined as some string
                if (string[grid->elem3d[ie].string].phys_flag == OFF) { // make sure physics is not turned off
                    continue;
                }
        }
            
        weight = 1. / 3. * grid->elem3d[ie].djac;
        for (inode=0; inode<grid->elem3d[ie].nnodes; inode++) {
            gnode = grid->elem3d[ie].nodes[inode]; // global node id
            node_var[gnode] += weight * elem3d_var[ie];
            weight_sum[gnode] += weight;
        }
    }
    
    // divide by the total weight
    for(gnode=0; gnode < grid->nnodes; gnode++)  {
        node_var[gnode] /= weight_sum[gnode];
    }
    
    // free arrays
    weight_sum = (double *) tl_free(sizeof(double), grid->nnodes, weight_sum);
}


