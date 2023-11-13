/* this routine compares two nodes */

// From Lucas ::
/*
 In the coupled branch, the processor id in cstorm_comm is not necessarily the processor id for the model comm. This means we can't use the MPI_Comm_rank call as written in the old code.
 
 For qsort to work properly, the comparison function can't take in anything other than the data type being sorted. This means we can't send the correct model comm for MPI_Comm_rank or the myid for the processor directly. I had to include a myid in the node struct to get around this.
 
 You may be able to eliminate the myid variable and use node->myid directly. That part of the struct is set somewhere in the new renumber or new node routines.
 */


#include "global_header.h"

#define NODE1_LEAST -1
#define NODE2_LEAST  1

int node_cmp(
             SNODE *node1,  /* the first node */
             SNODE *node2   /* the second node */
)
{
    
    //printf("node_cmp: node1: %d  node2: %d\n",node1->id,node2->id);
    
    /* first check for unused nodes */
    if (node1->string != UNSET_INT && node2->string == UNSET_INT) {
        return (NODE1_LEAST);
    }
    else if (node1->string == UNSET_INT && node2->string != UNSET_INT) {
        return (NODE2_LEAST);
    }
    else if (node1->string == UNSET_INT && node2->string == UNSET_INT) {
        if (node1->id < node2->id)
            return (NODE1_LEAST);
        else
            return (NODE2_LEAST);
    }
    /* if both nodes are in use */
    else {
#ifdef _MESSG
        
        int myid, ierr_code;
        myid=node1->myid;
        
        if (node1->resident_pe == myid && node2->resident_pe != myid)
            return (NODE1_LEAST);
        else if (node1->resident_pe != myid && node2->resident_pe == myid)
            return (NODE2_LEAST);
        else if (node1->resident_pe == myid && node2->resident_pe == myid)
            if (node1->id < node2->id) // cjt :: (node1->gid < node2->gid)
                return (NODE1_LEAST);
            else
                return (NODE2_LEAST);
            else if (node1->resident_pe != node2->resident_pe)
                if (node1->resident_pe < node2->resident_pe)
                    return (NODE1_LEAST);
                else
                    return (NODE2_LEAST);
                else if (node1->resident_id < node2->resident_id) // cjt :: (node1->gid < node2->gid)
                    return (NODE1_LEAST);
                else
                    return (NODE2_LEAST);
        
#else
        if (node1->id < node2->id)
            return (NODE1_LEAST);
        else
            return (NODE2_LEAST);
#endif
    }
}
