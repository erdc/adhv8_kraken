/* this routine renumbers the int arrays */
#include "global_header.h"

void node_renumber_snode(
        int max_nnode,      // the newly incremented maximum number of nodes
        SNODE *nodes,       // the node array to be renumbered
        int *new_number,    // the new node numbers
        int *order_tmp) 
{
    int i;                 
    
        
#ifdef _MESSG_ADPT 
        /* reorder global node numbers */
        GLOBAL_NODE *gtmp = (GLOBAL_NODE *) tl_alloc(sizeof(GLOBAL_NODE), max_nnode);
        node_renumber_gn(node_pair, gtmp, new_node_number, order_tmp);
        node_renumber_gn(node_ladj, gtmp, new_node_number, order_tmp);
        node_renumber_gn(node_radj, gtmp, new_node_number, order_tmp);
        gtmp = (GLOBAL_NODE *) tl_free(sizeof(GLOBAL_NODE), max_nnode, gtmp);
#else

    // renumber parent ids
   /* 
    int *itmp_array = (int *) tl_alloc(sizeof(int), max_nnode);
    sarray_int_init(itmp_array, max_nnode);
    for (i = 0; i < max_nnode; i++) {
        if (new_number[i] != i) {
            itmp_array[new_number[i]] = nodes[new_number[i]].parent[0];
            if (order_tmp[i]==0) {
                nodes[new_number[i]].parent[0] = nodes[i].parent[0];
            } else {
                nodes[new_number[i]].parent[0] = itmp_array[i];
            }
        }
    }
    sarray_int_init(itmp_array, max_nnode);
    for (i = 0; i < max_nnode; i++) {
        if (new_number[i] != i) {
            itmp_array[new_number[i]] = nodes[new_number[i]].parent[1];
            if (order_tmp[i]==0) {
                nodes[new_number[i]].parent[1] = nodes[i].parent[1];
            } else {
                nodes[new_number[i]].parent[1] = itmp_array[i];
            }
        }
    }
   */ 



    // renumber node struct
    SNODE *tmp_array = (SNODE *) tl_alloc(sizeof(SNODE), max_nnode);
    snode_init_array(tmp_array, max_nnode);
    for (i = 0; i < max_nnode; i++) {
      if ((new_number[i] != i) && (new_number[i] >= 0) && (new_number[i] < max_nnode)) {
            snode_copy(&(tmp_array[new_number[i]]), nodes[new_number[i]]);
            if (order_tmp[i]==0) {
               snode_copy(&(nodes[new_number[i]]), nodes[i]);
            } else {
               snode_copy(&(nodes[new_number[i]]), tmp_array[i]);
            }
        }
    }

    // re-order node ids
    for (i = 0; i < max_nnode; i++) { 
      nodes[i].id = i;

    }

    tmp_array = (SNODE *) tl_free(sizeof(SNODE), max_nnode, tmp_array);
    
#endif

}
