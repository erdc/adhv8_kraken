#include "global_header.h"

/***********************************************************/
/***********************************************************/
/***********************************************************/

int is_my_node(int global_id, SGRID *g) {

    int local_index;

    for (local_index = 0; local_index < g->nnodes; local_index++)
    {
        if (g->node[local_index].gid == global_id){
            return local_index;
        }
    }
    return -1;
}
