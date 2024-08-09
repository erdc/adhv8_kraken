#include "adh.h"
/***********************************************************/
/***********************************************************/
/***********************************************************/

int add_node(NODE_LIST_ITEM ** head, int local_index, int global_index, SVECT coords, int sd, int rnode, int global_surf)
{
    NODE_LIST_ITEM *entry;
    
    entry = tl_alloc(sizeof(NODE_LIST_ITEM), 1);
    entry->rnode = rnode; /*local index on owning processor */
    entry->sd = sd; /* owning processor */
#ifdef _MESSG
    entry->global = global_index;
    entry->global_surf = global_surf;
    VT_3D_VCOPY(coords, entry->coords);
#endif
    entry->local = local_index;

    entry->next = *head;
    *head = entry;
    
    return true;
}

int add_node_tmp(NODE_LIST_ITEM ** head, int local_index, int global_index) {
    NODE_LIST_ITEM *entry;
    entry = tl_alloc(sizeof(NODE_LIST_ITEM), 1);
#ifdef _MESSG
    entry->global = global_index;
#endif
    entry->local = local_index;
    entry->next = *head;
    *head = entry;
    
    return true;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

int add_ghost_node(NODE_LIST_ITEM **head, int index, int my_nnode, int *num_ghosts)
{
    NODE_LIST_ITEM *entry;
    int local_index = UNSET_INT;
    int hashkey;
#ifdef _MESSG
    /*if ((local_index = search_ghost_node(head, index)) >= 0)
     {
     return local_index;
     }*/
    entry = tl_alloc(sizeof(NODE_LIST_ITEM), 1);
    local_index = my_nnode + *num_ghosts;
    entry->global = index;
    entry->local = local_index;
    entry->next = *head;
    *head = entry;
    (*num_ghosts)++;
#endif
    return local_index;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

int search_ghost_node(NODE_LIST_ITEM * head, int index)
{
#ifdef _MESSG
    NODE_LIST_ITEM *entry = head;
    while (entry->next != NULL)
    {
        if (entry->global == index)
        {
            return entry->local;
        }
        entry = entry->next;
    }
#endif
    return -1;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

int add_elem2d_list(ELEM2D_LIST_ITEM ** head, int index, int *node, int iel_mat)
{
    ELEM2D_LIST_ITEM *elem_list;
    int i;
    
    /* create an entry for this element */
    elem_list = tl_alloc(sizeof(ELEM2D_LIST_ITEM), 1);
    elem_list->ielem = index;
    elem_list->nd1 = node[0];
    elem_list->nd2 = node[1];
    elem_list->nd3 = node[2];
    elem_list->mat = iel_mat;
    
    elem_list->next = *head;
    *head = elem_list;
    
    return true;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

int add_elem3d_list(ELEM3D_LIST_ITEM ** head, int global_index, int *node, int iel_mat)
{
    ELEM3D_LIST_ITEM *elem_list;
    int i;
    /*  int hashkey;
     
     hashkey = global_index % NEWIO_HASHSIZE;
     */
    /* create an entry for this element */
    elem_list = tl_alloc(sizeof(ELEM3D_LIST_ITEM), 1);
    elem_list->ielem = global_index;
    elem_list->nd1 = node[0];
    elem_list->nd2 = node[1];
    elem_list->nd3 = node[2];
    elem_list->nd4 = node[3];
    elem_list->mat = iel_mat;
    
    /* elem_list->next = head[hashkey];
     head[hashkey] = elem_list;*/
    elem_list->next = *head;
    *head = elem_list;
    
    return true;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

int free_elem3d_list(ELEM3D_LIST_ITEM * elem3d_list)
{
    ELEM3D_LIST_ITEM *ptr, *next;
    int i;
    
    /*  for (i = 0; i < NEWIO_HASHSIZE; i++)
     {
     ptr = elem3d_list[i];*/
    ptr  = elem3d_list;
    while (ptr->next != NULL)
    {
        next = ptr->next;
        ptr = tl_free(sizeof(ELEM3D_LIST_ITEM), 1, ptr);
        ptr = next;
    }
    ptr = tl_free(sizeof(ELEM3D_LIST_ITEM), 1, ptr);
    /* }
     elem3d_list = tl_free(sizeof(ELEM3D_LIST_ITEM *), NEWIO_HASHSIZE, elem3d_list);
     */
    return true;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

int free_elem2d_list(ELEM2D_LIST_ITEM * elem2d_list)
{
    ELEM2D_LIST_ITEM *ptr, *next;
    int i;
    
    ptr = elem2d_list;
    while (ptr->next != NULL)
    {
        next = ptr->next;
        ptr = tl_free(sizeof(ELEM2D_LIST_ITEM), 1, ptr);
        ptr = next;
    }
    ptr = tl_free(sizeof(ELEM2D_LIST_ITEM), 1, ptr);
    
    return true;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

int free_node_list(NODE_LIST_ITEM * head)
{
    NODE_LIST_ITEM *next;
    
    while (head->next != NULL)
    {
        next = head->next;
        head = tl_free(sizeof(NODE_LIST_ITEM), 1, head);
        head = next;
    }
    head = tl_free(sizeof(NODE_LIST_ITEM), 1, head);
    
    return true;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

int free_ghost_nodes(NODE_LIST_ITEM *** ghost_nodes, int npes)
{
#ifdef _MESSG
    NODE_LIST_ITEM *ptr, *next;
    int i, j;
    
    for (i = 0; i < npes; i++)
    {
        for (j = 0; j < NEWIO_HASHSIZE; j++)
        {
            ptr = ghost_nodes[i][j];
            while (ptr->next != NULL)
            {
                next = ptr->next;
                ptr = tl_free(sizeof(NODE_LIST_ITEM), 1, ptr);
                ptr = next;
            }
            ptr = tl_free(sizeof(NODE_LIST_ITEM), 1, ptr);
        }
        ghost_nodes[i] = tl_free(sizeof(NODE_LIST_ITEM *), NEWIO_HASHSIZE, ghost_nodes[i]);
    }
    ghost_nodes = tl_free(sizeof(NODE_LIST_ITEM **), npes, ghost_nodes);
#endif
    return true;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

int add_surf_node(NODE_LIST_ITEM ** head, int local_index)
{
    NODE_LIST_ITEM *entry;

    entry = tl_alloc(sizeof(NODE_LIST_ITEM), 1);
    entry->local = local_index;
    entry->next = *head;
    *head = entry;

    return true;
}

