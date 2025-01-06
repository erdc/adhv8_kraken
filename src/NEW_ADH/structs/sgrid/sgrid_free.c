#include "adh.h"
static int DEBUG = OFF;
static int DEBUG_WITH_PICKETS = OFF;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees and AdH grid
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] pgrid          (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void sgrid_free(SGRID *g) {
    int ie, inode;
    
    if (g->elem1d != NULL && g->nelems1d > 0) {selem1d_free_array(g->elem1d, g->max_nelems1d);}
    if (g->elem2d != NULL && g->nelems2d > 0) {selem2d_free_array(g->elem2d, g->max_nelems2d);}
    if (g->elem3d != NULL && g->nelems3d > 0) {selem3d_free_array(g->elem3d, g->max_nelems3d);}
    if (g->nnodes > 0) {
        for (inode=0; inode<g->max_nnodes; inode++) {snode_free(&(g->node[inode]));}
        g->node = (SNODE *) tl_free(sizeof(SNODE), g->max_nnodes, g->node);
    }
    if(g->inv_per_node!=NULL){g->inv_per_node = (int *) tl_free(sizeof(int), g->nnodes, g->inv_per_node);}
    g = (SGRID *) tl_free(sizeof(SGRID), 1, g);
}
