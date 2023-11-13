#include "global_header.h"

int is_my_elem(int node1, int node2, SGRID *g) {
    int i, found1=0,found2=0;

    for(i=0;i<g->nelems2d;i++){
        if((g->elem2d[i].nodes[0] == node1) || (g->elem2d[i].nodes[1] == node1) || (g->elem2d[i].nodes[2] == node1)) found1=1;

        if((g->elem2d[i].nodes[0] == node2) || (g->elem2d[i].nodes[1] == node2) || (g->elem2d[i].nodes[2] == node2)) found2=1;

        if((found1==1) && (found2==1)){
            break;
        }else{
            found1=0;
            found2=0;
        }
    }

    if((found1==1) && (found2==1)){
        return 1;
    }else{
        return -1;
    }
}
