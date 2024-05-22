/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  tl_assign_element2partition.c */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Give multiple PEs of partitioned nodes, created an arbitrary ruleset for assigned elements to PEs
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] pgrid          (SGRID *)  pointer to an AdH grid
 *
 * \note CJT: Used for calculating masses across PEs
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void sgrid_read(SGRID *g) {
    int ie, node_resid_count;
    
    for (ie=0; ie<g->nelems2d; ie++) {
        
        
        node_resid_count = 0;
        for (inode=0;inode<g->elem2d[ie].nnodes;inode++) {
            if (g->node[ g->elem2d[ie].node[inode] ].myid == g->smpi->myid) {
                node_resid_count++;
            }
        }
        
        // all nodes are resident
        if (node_resid_count == g->elem2d[ie].nnodes) {
            g->elem2d[ie].myid = g->smpi->myid;
        }
        
        // no nodes are resident
        if (node_resid_count == 0) {
            continue;
        }
        
        // at least one ghost node is on this element
        // assign to PE that owns the first node
        
        
    }
    
    
}
