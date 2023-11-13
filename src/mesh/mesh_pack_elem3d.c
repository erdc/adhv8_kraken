/* pack the 3d elements in the mesh region indicated by the node list */
// CJT :: Hard Wired for TETS!!!!!
#include "global_header.h"

#ifdef _MESSG
void mesh_pack_elem3d(
                      int *pack_node,		/* the new partition of the nodes I own */
                      MESSG_BUFFER * ibuffer,		/* the integer buffer to be packed */
                      MESSG_BUFFER * dbuffer,   /* the double buffer to be packed */
                      SGRID *g
                      )
{
    int ie;			/* loop counter over the elements */
    int ibuff, dbuff;			/* position in buffer */
    int nelem3d_out;		/* the number of 3d elements going out */
    int nelem3d_idatum;		/* amount of integer data per 3d element */
    int nelem3d_ddatum;   /* amount of double data per 3d element */
    
    /* count the number of elements going out */
    for(ie = 0, nelem3d_out = 0; ie < g->nelems3d; ie++)
        if(pack_node[g->elem3d[ie].nodes[0]] == YES || pack_node[g->elem3d[ie].nodes[1]] == YES ||
           pack_node[g->elem3d[ie].nodes[2]] == YES || pack_node[g->elem3d[ie].nodes[3]] == YES)
            nelem3d_out++;
    
    /* allocate the element buffers */
    nelem3d_idatum = elem3d_pack_icnt();
    nelem3d_ddatum = 26; /* x y z coordinates for all 4 nodes */
    messg_buffer_alloc(nelem3d_out * nelem3d_idatum, sizeof(int), ibuffer);
    messg_buffer_alloc(nelem3d_out * nelem3d_ddatum, sizeof(double), dbuffer);
    
    /* set the buffer type */
    ibuffer->type = MESSG_INT;
    dbuffer->type = MESSG_DOUBLE;
    
    /* load the element buffers */
    for(ie = 0, ibuff = 0, dbuff = 0; ie < g->nelems3d; ie++){
        if(pack_node[g->elem3d[ie].nodes[0]] == YES || pack_node[g->elem3d[ie].nodes[1]] == YES ||
           pack_node[g->elem3d[ie].nodes[2]] == YES || pack_node[g->elem3d[ie].nodes[3]] == YES) {
            
            elem3d_packi(((int *)ibuffer->buffer) + ibuff, ie, g);
            ibuff += nelem3d_idatum;
            elem3d_packd(((double *)dbuffer->buffer) + dbuff, ie, g);
            dbuff += nelem3d_ddatum;
        }
    }
    
}
#else
void mesh_pack_elem3d(
                      void
                      )
{
}
#endif
