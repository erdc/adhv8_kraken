/* this routine returns the node on an edge */

#include "global_header.h"

int adpt_get_node(
                  SMODEL *mod,
                  EDGE_LIST_ITEM * edge_pntr,	/* pointer to the edge */
                  EDGE_LIST_ITEM ** edge_hashtab	/* edge hash table */
)
{
    int new_node,ind_sur,ind_bed;			/* the new node */
    int nd1;			/* the 1st node on the edge */
    int nd2;			/* the 2nd node on the edge */
    EDGE_LIST_ITEM *flip_edge_pntr;	/* the pointer to the edge with the nodes flipped */
    
    SSW_3D *sw;
    SNS_3D *ns;
    if (mod->flag.SW_FLOW
#ifdef _ADH_GROUNDWATER
#ifdef _DWGW_COUPLING
                    && (mod->flag.GW_FLOW==OFF)
#endif
#endif
            ) {
        sw = mod->sw->d3;
    } else if (mod->flag.NS_FLOW) {
        ns = mod->ns->d3;
    }
    double depth;
    
    /* report an error if the edge is not found */
    if(edge_pntr == NULL) {
        tl_error("Edge not found in adpt_get_node.");
    }
    
    /* if the node exists, then return */
    if(edge_pntr->new_node != UNSET_INT) {
        return (edge_pntr->new_node);
    }
    
    /* the node does not exist, create the node and enter it in the hash table */
    new_node = node_new(mod);
    if(((mod->grid->node[edge_pntr->nd1].global_surf_id != UNSET_INT) && (mod->grid->node[edge_pntr->nd2].global_surf_id != UNSET_INT)) && (mod->grid->ndim==3)) {
        
        ind_sur = node_new_surface(mod);
        
        //mod->grid->old_global_surf[ind_sur] = gn.global_surf_id;
        mod->grid->nodeID_2d_to_3d_sur[ind_sur] = new_node;
        mod->grid->nodeID_3d_to_2d_sur[new_node] = ind_sur;
        
        if (mod->flag.SW_FLOW
#ifdef _ADH_GROUNDWATER
#ifdef _DWGW_COUPLING
                    && (mod->flag.GW_FLOW==OFF)
#endif
#endif
                    ) {
            ssw_3d_node_avg_sur(mod->sw->d3, ind_sur, mod->grid->nodeID_3d_to_2d_sur[edge_pntr->nd1], mod->grid->nodeID_3d_to_2d_sur[edge_pntr->nd2], mod->grid);
        } else if (mod->flag.NS_FLOW) {
            //sns_3d_node_avg_sur(mod->ns->d3, ind_sur, mod->grid->nodeID_3d_to_2d_sur[edge_pntr->nd1], mod->grid->nodeID_3d_to_2d_sur[edge_pntr->nd2], mod->grid);
        }
        
    }
    if((mod->grid->node[edge_pntr->nd1].global_bed_id != UNSET_INT) && (mod->grid->node[edge_pntr->nd2].global_bed_id != UNSET_INT) && (mod->grid->ndim==3)) {
        
        
        ind_bed = node_new_bed(mod);
        
        //mod->grid->old_global_bed[ind_bed] = gn.global_bed_id;
        mod->grid->nodeID_2d_to_3d_bed[ind_bed] = new_node;
        mod->grid->nodeID_3d_to_2d_bed[new_node] = ind_bed;
        
        if (mod->flag.SW_FLOW
#ifdef _ADH_GROUNDWATER
#ifdef _DWGW_COUPLING
                    && (mod->flag.GW_FLOW==OFF)
#endif
#endif
                    ) {
            ssw_3d_node_avg_bed(mod->sw->d3, ind_bed, mod->grid->nodeID_3d_to_2d_bed[edge_pntr->nd1], mod->grid->nodeID_3d_to_2d_bed[edge_pntr->nd2], mod->grid);
        } else if (mod->flag.NS_FLOW) {
            //sns_3d_node_avg_bed(mod->ns->d3, ind_bed, mod->grid->nodeID_3d_to_2d_bed[edge_pntr->nd1], mod->grid->nodeID_3d_to_2d_bed[edge_pntr->nd2], mod->grid);
        }
        
    }
    
    mod->grid->node[new_node].id = new_node;
    node_avg(mod, edge_pntr->nd1, edge_pntr->nd2, new_node);
    flip_edge_pntr = edge_hash_lookup(edge_pntr->nd2, edge_pntr->nd1, edge_hashtab);
    edge_pntr->new_node = new_node;
    flip_edge_pntr->new_node = new_node;
    
    /* set the adjacencies for the new node */
    nd1 = edge_pntr->nd1;
    nd2 = edge_pntr->nd2;
    
#ifdef _MESSG
    mod->grid->node[new_node].parent_res_id[0] = mod->grid->node[nd1].resident_id;
    mod->grid->node[new_node].parent_res_pe[0] = mod->grid->node[nd1].resident_pe;
    mod->grid->node[new_node].parent_res_id[1] = mod->grid->node[nd2].resident_id;
    mod->grid->node[new_node].parent_res_pe[1] = mod->grid->node[nd2].resident_pe;
#endif
    mod->grid->node[new_node].parent[0] = nd1;
    mod->grid->node[new_node].parent[1] = nd2;
    
    
    /* temporarily set the adjacencies for nd1 with "global ghost node numbers"  */
#ifdef _MESSG
    if((mod->grid->node[nd1].parent_res_pe[0] == mod->grid->node[nd2].resident_pe) &&
       (mod->grid->node[nd1].parent_res_id[0] == mod->grid->node[nd2].resident_id)){
        mod->grid->node[nd1].parent_res_pe[0] = mod->grid->smpi->myid;
        mod->grid->node[nd1].parent_res_id[0] = new_node;
    }else if((mod->grid->node[nd1].parent_res_pe[1] == mod->grid->node[nd2].resident_pe) &&
             (mod->grid->node[nd1].parent_res_id[1] == mod->grid->node[nd2].resident_id)){
        mod->grid->node[nd1].parent_res_pe[1] = mod->grid->smpi->myid;
        mod->grid->node[nd1].parent_res_id[1] = new_node;
    }
#endif
    
    if(mod->grid->node[nd1].parent[0] == nd2) {
        mod->grid->node[nd1].parent[0] = new_node;
    } else if(mod->grid->node[nd1].parent[1] == nd2) {
        mod->grid->node[nd1].parent[1] = new_node;
    }
    
    /* temporarily set the adjacencies for nd2 with "global ghost node numbers"  */
#ifdef _MESSG
    if((mod->grid->node[nd2].parent_res_pe[0] == mod->grid->node[nd1].resident_pe) &&
       (mod->grid->node[nd2].parent_res_id[0] == mod->grid->node[nd1].resident_id)){
        mod->grid->node[nd2].parent_res_pe[0] = mod->grid->smpi->myid;
        mod->grid->node[nd2].parent_res_id[0] = new_node;
    }else if((mod->grid->node[nd2].parent_res_pe[1] == mod->grid->node[nd1].resident_pe) &&
             (mod->grid->node[nd2].parent_res_id[1] == mod->grid->node[nd1].resident_id)){
        mod->grid->node[nd2].parent_res_pe[1] = mod->grid->smpi->myid;
        mod->grid->node[nd2].parent_res_id[1] = new_node;
    }
#endif
    if(mod->grid->node[nd2].parent[0] == nd1) {
        mod->grid->node[nd2].parent[0] = new_node;
    } else if(mod->grid->node[nd2].parent[1] == nd1) {
        mod->grid->node[nd2].parent[1] = new_node;
    }
    
    /* return the new node */
    return (new_node);
}
