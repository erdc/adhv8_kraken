#include "global_header.h"

//**************************************************************************************************************//
//**************************************************************************************************************//

//**************************************************************************************************************//
//**************************************************************************************************************//

// noptions :: the number of integer options associated with the elem property
// lowerORhigher :: choice to choose either lower (-1) or higher (0) option at ambiguous corner node

void elem2d_to_node_int(SGRID *grid, int *nodal_var, int *elem2d_var, int noptions, int lowerORhigher) {

    int ie=0, lnode=0, ioption=0;
    int node[NDONTRI];

    int lower = -1;
    int higher = +1;

    assert(nodal_var);
    assert(elem2d_var);
    
    if (lowerORhigher != lower && lowerORhigher != higher) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> input lower or higher option is not valid\n");
    }

    for (ie=0; ie < grid->nelems2d; ie++) {
        for (lnode=0; lnode<NDONTRI; lnode++) {
            node[lnode] = grid->elem2d[ie].nodes[lnode];
            nodal_var[node[lnode]] = UNSET_INT;
        }
    }

    for (ie=0; ie < grid->nelems2d; ie++) {                 // loop over elements
        for (ioption=0; ioption<noptions; ioption++) {        // loop over integer options
            if (elem2d_var[ie] == ioption) {
                for (lnode=0; lnode<NDONTRI; lnode++) {     // loop over local nodes
                    node[lnode] = grid->elem2d[ie].nodes[lnode];
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

void elem2d_to_node_double(SGRID *grid, STR_VALUE *string, double *node_var, double *elem2d_var) {
    
    int ie=0, inode=0, gnode=0;
    double weight = 0.;
    double elem_sfield[NDONTRI];
    
    assert(node_var);
    assert(elem2d_var);
    
    // allocate array to hold weight sum
    double *weight_sum = (double *) tl_alloc(sizeof(double), grid->nnodes);
    
    sarray_init_dbl(node_var, grid->nnodes);
    sarray_init_dbl(weight_sum, grid->nnodes);
    for(ie = 0; ie < grid->nelems2d; ie++) {
        if (string[grid->elem2d[ie].string].phys_flag != OFF) {
            
            weight = 1. / 3. * grid->elem2d[ie].djac;
            
            for (inode=0; inode<NDONTRI; inode++) {
                gnode = grid->elem2d[ie].nodes[inode]; // global node id
                node_var[gnode] += weight * elem2d_var[ie];
                weight_sum[gnode] += weight;
            }
        }
    }
    
    // divide by the total weight
    for(gnode=0; gnode < grid->nnodes; gnode++)  {
        node_var[gnode] /= weight_sum[gnode];
    }
    
    // free arrays
    weight_sum = (double *) tl_free(sizeof(double), grid->nnodes, weight_sum);
}


//**************************************************************************************************************//
//**************************************************************************************************************//

// noptions :: the number of integer options associated with the elem property
// lowerORhigher :: choice to choose either lower (-1) or higher (0) option at ambiguous corner node

// Maps element material properties to nodal data on the bed
// can be used for sediment friction flag and coefficient

void elem2dbed_to_node_int(SGRID *grid, int *nodal_var, int *elem2d_var, int noptions, int lowerORhigher) {
    
    int ie=0, lnode=0, ioption=0, id2d_bed=UNSET_INT, id3d=UNSET_INT;
    int node[NDONTRI];
    
    int lower = -1;
    int higher = +1;
    
    assert(nodal_var);
    assert(elem2d_var);
    
    if (lowerORhigher != lower && lowerORhigher != higher) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> input lower or higher option is not valid\n");
    }
    
    // note:
    // nodal_var = nnodes_bed
    // elem2d_var = nelems2d_bed
    
    for (ie=0; ie < grid->nelems2d_bed; ie++) {
        for (lnode=0; lnode<NDONTRI; lnode++) {
            id3d = grid->elem2d[ grid->elem2d_bed[ie] ].nodes[lnode];
            id2d_bed = grid->nodeID_3d_to_2d_bed[id3d];
            node[lnode] = id2d_bed;
            nodal_var[node[lnode]] = UNSET_INT;
        }
    }
    
    for (ie=0; ie < grid->nelems2d_bed; ie++) {                 // loop over elements
        for (ioption=0; ioption<noptions; ioption++) {        // loop over integer options
            if (elem2d_var[ie] == ioption) {
                for (lnode=0; lnode<NDONTRI; lnode++) {     // loop over local nodes
                    id3d = grid->elem2d[ grid->elem2d_bed[ie] ].nodes[lnode];
                    id2d_bed = grid->nodeID_3d_to_2d_bed[id3d];
                    node[lnode] = id2d_bed;
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

void elem2dbed_to_node_double(SGRID *grid, STR_VALUE *string, double *node_var, double *elem2d_var) {
    
    int ie=0, inode=0, gnode=0, ie2d=UNSET_INT, id3d=UNSET_INT, id2d_bed=UNSET_INT;
    double weight = 0.;
    double elem_sfield[NDONTRI];
    
    assert(node_var);
    assert(elem2d_var);
    
    // note:
    // nodal_var = nnodes_bed
    // elem2d_var = nelems2d_bed
    
    // allocate array to hold weight sum
    double *weight_sum = (double *) tl_alloc(sizeof(double), grid->nnodes_bed);
    
    sarray_init_dbl(node_var, grid->nnodes_bed);
    sarray_init_dbl(weight_sum, grid->nnodes_bed);
    for(ie = 0; ie < grid->nelems2d_bed; ie++) {
        
        ie2d = grid->elem2d_bed[ie]; // map back to list of full 2d elements
        
        if (string[grid->elem2d[ie2d].string].phys_flag != OFF) {
            
            weight = 1. / 3. * fabs(grid->elem2d[ie2d].djac); // take absolute here because 2d bottom element has opposite node numbering
            //printf("ie: %d weight: %20.10f  elem2d_var: %20.10f\n",ie, weight,elem2d_var[ie]);
            
            for (inode=0; inode<NDONTRI; inode++) {
                id3d = grid->elem2d[ie2d].nodes[inode];
                id2d_bed = grid->nodeID_3d_to_2d_bed[id3d];
                node_var[id2d_bed] += weight * elem2d_var[ie];
                weight_sum[id2d_bed] += weight;
            }
        }
    }
    
    // divide by the total weight
    for(gnode=0; gnode < grid->nnodes_bed; gnode++)  {
        node_var[gnode] /= weight_sum[gnode];
    }
    
    // free arrays
    weight_sum = (double *) tl_free(sizeof(double), grid->nnodes_bed, weight_sum);
}
