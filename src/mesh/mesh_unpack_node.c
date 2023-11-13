#include "global_header.h"

/* unpacks the node buffers */
void mesh_unpack_node(MESSG_BUFFER * node_dbuffer,  /* the buffer for the double part of the node message */
                      MESSG_BUFFER * node_dbuffer_surface,  /* the buffer for the double part of the node message */
                      MESSG_BUFFER * node_dbuffer_bed,  /* the buffer for the double part of the node message */
                      MESSG_BUFFER * node_ibuffer,   /* the buffer for the integer part of the node message */
                      SMODEL *mod,
                      NODE_LIST_ITEM ** node_hashtab, /* tops of the linked lists for the node hash table */
                      int flag                        /* 0 for comm_node_data, 1 otherwise */
)
{
    int icnt_ibuff = 0;           /* counter in the integer buffer */
    int icnt_dbuff = 0;           /* counter in the double buffer */
    int icnt_dbuff_surface = 0;           /* counter in the double buffer */
    int icnt_dbuff_bed = 0;           /* counter in the double buffer */
    int local_id = 0;             /* the local_id node number */
    int nd_ddatum = 0;            /* amount of double data per node */
    int nd_ddatum_surface = 0;            /* amount of double surface data per node */
    int nd_ddatum_bed = 0;            /* amount of double bed data per node */
    int nd_idatum = 0;            /* amount of integer data per node */
    int *ibuffer = NULL;          /* Temporary Array */
    int i;                /* loop counter */
    int local_surface;
    int local_spot = 0;
    SGRID *g = mod->grid;
    int next_surf_spot=g->my_nnodes_sur;
    int next_bed_spot = g->my_nnodes_bed;
    int local_bed;
    SNODE gn;               /* the global node */
    /* set the packet size */
    nd_idatum = node_pack_icnt();
    nd_ddatum = node_pack_dcnt(mod);
    
    if(g->type == COLUMNAR){
        nd_ddatum_surface = node_pack_dcnt_surface(mod);
        nd_ddatum_bed = node_pack_dcnt_bed(mod);
    }
    snode_init(&gn);
    /* unpacks the buffers */
    ibuffer = (int *) node_ibuffer->buffer;
    for (icnt_ibuff = 0, icnt_dbuff = 0; icnt_ibuff < node_ibuffer->nitem; icnt_ibuff += nd_idatum, icnt_dbuff += nd_ddatum)
    {
        /* Each Time, look forward in the buffer to see what node is coming */
        /* We have to make sure that we keep up with the buffer location */
        gn.string = ibuffer[icnt_ibuff];
        gn.edge_string = ibuffer[icnt_ibuff+1];
        gn.node_string = ibuffer[icnt_ibuff+2];
        gn.original_id = ibuffer[icnt_ibuff+3];
        gn.parent[0] = ibuffer[icnt_ibuff+4];
        gn.parent[1] = ibuffer[icnt_ibuff+5];
        gn.level = ibuffer[icnt_ibuff+10];
        gn.block = ibuffer[icnt_ibuff+11];
        gn.els_flag = ibuffer[icnt_ibuff+12];
        gn.nelems_connected = ibuffer[icnt_ibuff+13];
        gn.gid = ibuffer[icnt_ibuff+14];
        gn.global_surf_id = ibuffer[icnt_ibuff+15];
        gn.global_bed_id = ibuffer[icnt_ibuff+16];
        gn.resident_pe = ibuffer[icnt_ibuff+17];
        gn.resident_id = ibuffer[icnt_ibuff+18];
        if(flag ==0){
            local_id = ibuffer[icnt_ibuff + 18];
        }
        else {
            local_id = node_get_local(gn, node_hashtab, mod);
        }
        
        node_unpacki(((int *) node_ibuffer->buffer) + icnt_ibuff, local_id, mod->grid);
        node_unpackd(((double *) node_dbuffer->buffer) + icnt_dbuff, local_id, mod->sw,
#ifdef _ADH_GROUNDWATER
                mod->sgw,
#endif
                mod->grid , mod->ntransport, mod->con
#ifdef _SEDIMENT
                     , mod->sed,
                     mod->nconti,
                     mod->nsed,
                     mod->nlayers
#endif
                     );
        if((g->type == COLUMNAR) && g->node[local_id].global_surf_id != UNSET_INT) {
            
            local_surface=g->nodeID_3d_to_2d_sur[local_id];
            node_unpackd_surface(((double *) node_dbuffer_surface->buffer) + icnt_dbuff_surface, local_surface, mod->sw);
            icnt_dbuff_surface += nd_ddatum_surface;
        }
        if((g->type == COLUMNAR) && g->node[local_id].global_bed_id != UNSET_INT) {
            
            local_bed=g->nodeID_3d_to_2d_bed[local_id];
            
            node_unpackd_bed(((double *) node_dbuffer_bed->buffer) + icnt_dbuff_bed, local_bed, mod->sw
#ifdef _SEDIMENT
                             , mod->sed,
                             mod->nconti,
                             mod->nsed,
                             mod->nlayers
#endif
                             );
            icnt_dbuff_bed += nd_ddatum_bed;
        }
        
    }
    
    return;
}
