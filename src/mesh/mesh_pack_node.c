#include "global_header.h"

/*! 
 \brief Pack the given list of nodes
 */
void mesh_pack_node(int nnode_list, /* the number of nodes to be packed */
                    int *nodes_list,    /* the nodes to be packed */
                    MESSG_BUFFER * dbuffer, /* the double buffer to be packed */
                    MESSG_BUFFER * dbuffer_surface, /* the double buffer to be packed */
                    MESSG_BUFFER * dbuffer_bed, /* the double buffer to be packed */
                    MESSG_BUFFER * ibuffer,  /* the integer buffer to be packed */
                    SMODEL *mod
                    )
{
    int ii = 0;                   /* loop counters */
    int ibuff = 0;                /* position in buffer */
    int ibuff_bed = 0;
    int nd_ddatum = 0;            /* amount of double data per node */
    int nd_ddatum_surface = 0;            /* amount of surface double data per node */
    int nd_ddatum_bed = 0;            /* amount of bed double data per node */
    int nd_idatum = 0;            /* amount of integer data per node */
    int nint_node_buff = 0;       /* the number of integers sent for the nodes */
    int ndbl_node_buff = 0;       /* the number of doubles sent for the nodes */
    int ndbl_node_buff_surface = 0;       /* the number of surface doubles sent for the nodes */
    int ndbl_node_buff_bed = 0;       /* the number of bed doubles sent for the nodes */
    SGRID *g = mod->grid;
    
    /* set the packet size */
    nd_idatum = node_pack_icnt();
    nd_ddatum = node_pack_dcnt(mod);
    
    nint_node_buff = nnode_list * nd_idatum;
    ndbl_node_buff = nnode_list * nd_ddatum;
    if(g->type == COLUMNAR) {
        nd_ddatum_surface = node_pack_dcnt_surface(mod);
        nd_ddatum_bed = node_pack_dcnt_bed(mod);
        for(ii=0;ii<nnode_list;ii++){
            if(g->node[nodes_list[ii]].global_surf_id != UNSET_INT) ndbl_node_buff_surface += nd_ddatum_surface;
            if(g->node[nodes_list[ii]].global_bed_id != UNSET_INT) ndbl_node_buff_bed += nd_ddatum_bed;
        }
        ndbl_node_buff_surface += nd_ddatum_surface;
    }
    //tag(MPI_COMM_WORLD);
    
    /* allocate the element buffers */
    messg_buffer_alloc(nint_node_buff, sizeof(int), ibuffer);
    messg_buffer_alloc(ndbl_node_buff, sizeof(double), dbuffer);
    if((g->type == COLUMNAR) && (ndbl_node_buff_surface > 0)) messg_buffer_alloc(ndbl_node_buff_surface, sizeof(double), dbuffer_surface);
    if((g->type == COLUMNAR) && (ndbl_node_buff_bed > 0)) messg_buffer_alloc(ndbl_node_buff_bed, sizeof(double), dbuffer_bed);
    //tag(MPI_COMM_WORLD);
    
    /* set the buffer types */
    ibuffer->type = MESSG_INT;
    dbuffer->type = MESSG_DOUBLE;
    if((g->type == COLUMNAR) && (ndbl_node_buff_surface > 0)) dbuffer_surface->type = MESSG_DOUBLE;
    if((g->type == COLUMNAR) && (ndbl_node_buff_bed > 0)) dbuffer_bed->type = MESSG_DOUBLE;
    //tag(MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);
    
    /* load the nodal buffers */
    for (ii = 0, ibuff = 0; ii < nnode_list; ii++, ibuff += nd_idatum) {
        node_packi(((int *) ibuffer->buffer) + ibuff, nodes_list[ii], g);
    }
    
    //printf("myid: %d nnode_list: %d nd_ddatum: %d\n",g->smpi->myid,nnode_list,nd_ddatum);
    for (ii = 0, ibuff = 0; ii < nnode_list; ii++, ibuff += nd_ddatum){
        //printf("myid: %d node_list[%d]: %d\n",g->smpi->myid,ii,nodes_list[ii]);
        node_packd(((double *) dbuffer->buffer) + ibuff, nodes_list[ii], mod->sw,
#ifdef _ADH_GROUNDWATER
                mod->sgw,
#endif
                g, mod->ntransport, mod->con
#ifdef _SEDIMENT
                   , mod->sed,
                   mod->nconti,
                   mod->nsed,
                   mod->nlayers
#endif
                   );
    }
    
    if((g->type == COLUMNAR) && (ndbl_node_buff_surface > 0)) {
        ibuff=0;
        for(ii=0;ii<nnode_list;ii++){
            if(g->node[nodes_list[ii]].global_surf_id != UNSET_INT){
                node_packd_surface(((double *) dbuffer_surface->buffer) + ibuff, g->nodeID_3d_to_2d_sur[nodes_list[ii]], (mod->sw), g);
                ibuff += nd_ddatum_surface;
            }
            if(g->node[nodes_list[ii]].global_bed_id != UNSET_INT){
                node_packd_bed(((double *) dbuffer_bed->buffer) + ibuff_bed, g->nodeID_3d_to_2d_bed[nodes_list[ii]], (mod->sw), g
#ifdef _SEDIMENT
                               , mod->sed,
                               mod->nconti,
                               mod->nsed,
                               mod->nlayers
#endif
                               );
                ibuff_bed += nd_ddatum_bed;
            }
        }
    }
    
    return;
}
