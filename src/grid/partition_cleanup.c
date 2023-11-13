/* deletes unnecessary nodes */

#include "global_header.h"

#ifdef _MESSG
void partition_cleanup(SGRID *g,int myid) {
    
    int i, j;        /* loop counter */
    int ie;			 /* loop counter over the elements */
    int *free_node;  /* the free nodes */
    
    /* allocate and initialize the free node array */
    free_node = (int *)tl_alloc(sizeof(int), g->nnodes);
    for(i = 0; i < g->nnodes; i++)
        free_node[i] = YES;
    /* find the unused elements and initialize them,
     find the used elements and set the nodes to be kept */
    j=0;
    for(ie = 0; ie < g->nelems3d; ie++)
    {
        if( g->elem3d[ie].mat != UNSET_INT) //elem3d_flags in old trunk
        {
            if((g->node[g->elem3d[ie].nodes[0]].resident_pe != myid) &&
               (g->node[g->elem3d[ie].nodes[1]].resident_pe != myid) &&
               (g->node[g->elem3d[ie].nodes[2]].resident_pe != myid) &&
               (g->node[g->elem3d[ie].nodes[3]].resident_pe != myid) )
            {j++;
                selem3d_init(&(g->elem3d[ie]));
                g->elem3d[ie].id = ie;
            }
            else
            {
                free_node[g->elem3d[ie].nodes[0]] = NO;
                free_node[g->elem3d[ie].nodes[1]] = NO;
                free_node[g->elem3d[ie].nodes[2]] = NO;
                free_node[g->elem3d[ie].nodes[3]] = NO;
            }
        }
    }
    for(ie = 0; ie < g->nelems2d; ie++)
    {
        
        if(g->elem2d[ie].string != UNSET_INT) //elem2d_flags in old trunk
        {
            if((g->node[g->elem2d[ie].nodes[0]].resident_pe != myid) &&
               (g->node[g->elem2d[ie].nodes[1]].resident_pe != myid) &&
               (g->node[g->elem2d[ie].nodes[2]].resident_pe != myid) )
            {
                selem2d_init(&(g->elem2d[ie]));
                g->elem2d[ie].id = ie;
            }
            else
            {
                free_node[g->elem2d[ie].nodes[0]] = NO;
                free_node[g->elem2d[ie].nodes[1]] = NO;
                free_node[g->elem2d[ie].nodes[2]] = NO;
            }
        }
    }
    for(ie = 0; ie < g->nelems1d; ie++)
    {
        if(g->elem1d[ie].string != UNSET_INT)
        {
            if((g->node[g->elem1d[ie].nodes[0]].resident_pe != myid) &&
               (g->node[g->elem1d[ie].nodes[1]].resident_pe != myid) )
            {
                selem1d_init(&(g->elem1d[ie]));
                g->elem1d[ie].id = ie;
            }
            else
            {
                free_node[g->elem1d[ie].nodes[0]] = NO;
                free_node[g->elem1d[ie].nodes[1]] = NO;
            }
        }
    }
    /* initialize the unused nodes */
    ie=0;
    for(i = 0; i < g->nnodes; i++) {
        if(free_node[i] == YES) {//printf("MYID %d FREE node %d gid %d x y z %10.5e %10.5e %10.5e \n", g->smpi->myid, i, g->node[i].gid,g->node[i].x,g->node[i].y,g->node[i].z);
            snode_init(&(g->node[i]));
            g->node[i].id = i;
        }
    }
    /* free the memory in the free_node array */
    free_node = (int *)tl_free(sizeof(int), g->nnodes, free_node);
}
#else
void partition_cleanup(
                       void
                       )
{
}
#endif
