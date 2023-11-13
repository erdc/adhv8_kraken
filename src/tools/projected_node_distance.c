
#include "global_header.h"

double projected_node_distance(SGRID *grid1, int n1, SGRID *grid2, int n2){
    return ( fabs(grid2->node[n2].x - grid1->node[n1].x) + fabs(grid2->node[n2].y - grid1->node[n1].y) );
}
