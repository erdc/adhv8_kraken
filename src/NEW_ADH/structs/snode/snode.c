#include "adh.h"

/***********************************************************/
/***********************************************************/
/***********************************************************/
void snode_alloc_array(SNODE **node, int nnodes) {
    assert(nnodes > 0);
    (*node) = (SNODE *) tl_alloc(sizeof(SNODE), nnodes);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/* adding or removing variables here requires changes to node_new.c, node_pack_icnt.c, node_packi.c, and node_unpacki.c */
void snode_init(SNODE *node) {
    node->x = UNSET_FLT;
    node->y = UNSET_FLT;
    node->z = UNSET_FLT;
    node->string = UNSET_INT;
    node->edge_string =  NORMAL;
    node->node_string = NORMAL;
    node->els_flag = NORMAL;
    node->level = 0;
    node->block = UNSET_INT;
    node->area = UNSET_FLT;
    node->id = UNSET_INT;
    node->original_id = UNSET_INT;
    node->parent[0] = UNSET_INT;
    node->parent[1] = UNSET_INT;
    node->parent_res_id[0] = UNSET_INT;
    node->parent_res_id[1] = UNSET_INT;
    node->parent_res_pe[0] = UNSET_INT;
    node->parent_res_pe[1] = UNSET_INT;
    
    node->nelems_connected = 0;
    node->elemID = NULL;
    node->gid = UNSET_INT;
    node->global_surf_id = UNSET_INT;
    node->global_bed_id = UNSET_INT;
    node->bflag = UNSET_INT;

    node->myid = UNSET_INT;
#ifdef _MESSG
    node->resident_pe = UNSET_INT;
    node->resident_id = UNSET_INT;
#else
    node->resident_pe = 0;
    node->resident_id = UNSET_INT;
#endif
}

/***********************************************************/

void snode_init_array(SNODE *node, int nnodes) {
    int i = UNSET_INT;
    for (i = 0; i < nnodes; i++) {
        snode_init(&node[i]);
    }
}

/***********************************************************/
/***********************************************************/
void snode_init_alloc_array(SNODE **node, int nnodes) {
    snode_alloc_array(node, nnodes);
    snode_init_array(*node, nnodes);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
void snode_copy(SNODE *to, SNODE from) {
    to->x = from.x;
    to->y = from.y;
    to->z = from.z;
    to->string = from.string;
    to->edge_string = from.edge_string;
    to->node_string = from.node_string;
    to->els_flag = from.els_flag;
    to->level = from.level;
    to->block = from.block;
    to->area = from.area;
    to->id = from.id;
    to->original_id = from.original_id;
    to->parent[0] = from.parent[0];
    to->parent[1] = from.parent[1];
    to->parent_res_pe[0] = from.parent_res_pe[0];
    to->parent_res_pe[1] = from.parent_res_pe[1];
    to->parent_res_id[0] = from.parent_res_id[0];
    to->parent_res_id[1] = from.parent_res_id[1];
    to->nelems_connected = from.nelems_connected;
    to->elemID = from.elemID;
    to->global_surf_id = from.global_surf_id;
    to->global_bed_id = from.global_bed_id;
    to->bflag = from.bflag;

    to->myid = from.myid;
#ifdef _MESSG
    to->resident_pe = from.resident_pe;
    to->resident_id = from.resident_id;
    to->gid = from.gid;
#endif
}

/*---------------------------------------------------------*/
void snode_copy_array(SNODE *to, SNODE *from, int array_length) {
    int i;
    for (i=0; i<array_length; i++) {
        snode_copy(&to[i], from[i]);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
void snode_free(SNODE *node) {
    if (node->elemID != NULL) {
        node->elemID = tl_free(sizeof(int), node->nelems_connected, node->elemID);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
void snode_free_array(SNODE *node, int nnodes) {
    int i;
    for (i=0; i<nnodes; i++) {
        snode_free(&(node[i]));
    }
    node = (SNODE *) tl_free(sizeof(SNODE), nnodes, node);
}


/***********************************************************/
/***********************************************************/
/***********************************************************/

void snode_printScreen(SNODE nd) {
    printf("NODE ID: %d ----------------------------------------\n",nd.id);
    printf("global_id %d \n", nd.gid);
    printf("original ID: %d\n",nd.original_id);
    printf("global_surf_id %d \n", nd.global_surf_id);
    printf("global_bed_id %d \n", nd.global_bed_id);
    printf("position: {%10.5f,%10.5f,%10.5f}\n",nd.x,nd.y,nd.z);
    printf("string: %d \t edge string: %d\n",nd.string,nd.edge_string);
    printf("els_flag: %d\n",nd.els_flag);
    printf("level: %d \t block: %d\n",nd.level,nd.block);
    printf("parents %d %d \n", nd.parent[0],nd.parent[1]);
    printf("resident_pe: %d\n",nd.resident_pe);
    printf("myid: %d\n",nd.myid);
#ifdef _MESSG
    printf("parents res pe %d %d \n", nd.parent_res_pe[0],nd.parent_res_pe[1]);
    printf("parents res id %d %d \n", nd.parent_res_id[0],nd.parent_res_id[1]);
    printf("resident pe: %d \t resident id: %d\n",nd.resident_pe,nd.resident_id);
#endif
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void snode2svect(SNODE *nd, SVECT *v, int nnodes) {
    int i;
    for (i=0; i<nnodes; i++) {
        v[i].x = nd[i].x;
        v[i].y = nd[i].y;
        v[i].z = nd[i].z;
    }
}

void snode2svect2d(SNODE *nd, SVECT2D *v, int nnodes) {
    int i;
    for (i=0; i<nnodes; i++) {
        v[i].x = nd[i].x;
        v[i].y = nd[i].y;
    }
}
