/* this routine sets a node to be the average of two parent nodes */

#include "global_header.h"

void node_avg(SMODEL *mod, int nodeID_1, int nodeID_2, int the_node) {
    
    SNODE node1 = mod->grid->node[nodeID_1]; // alias
    SNODE node2 = mod->grid->node[nodeID_2]; // alias
    
    /* sets the node to be the average of the two parent nodes */
    mod->grid->node[the_node].x = (node1.x + node2.x) * 0.5;
    mod->grid->node[the_node].y = (node1.y + node2.y) * 0.5;
    mod->grid->node[the_node].z = (node1.z + node2.z) * 0.5;
    
    /* averages new nodes for physics arrays */
    if (mod->flag.SW2_FLOW) {ssw_2d_node_avg(mod->sw->d2, the_node, nodeID_1, nodeID_2);}
    if (mod->flag.SW3_FLOW) {ssw_3d_node_avg(mod->sw->d3, the_node, nodeID_1, nodeID_2, mod->grid);}
    if (mod->flag.NS3_FLOW) {sns_3d_node_avg(mod->ns->d3, the_node, nodeID_1, nodeID_2, mod->grid);}
#ifdef _ADH_GROUNDWATER
    if (mod->flag.GW_FLOW) {sgw_3d_node_avg(mod->sgw, the_node, nodeID_1, nodeID_2, mod->grid);}
#endif
    if (mod->flag.TRANSPORT) {scon_node_avg(mod->con, mod->ntransport, the_node, nodeID_1, nodeID_2);}
#ifdef _SEDIMENT
    if (mod->flag.SEDIMENT) {ssediment_node_avg(mod->sed, the_node, nodeID_1, nodeID_2, mod->grid);}
#endif
#ifdef _ADH_ICM
    //if (mod->flag.ICM) {}
#endif
    
    // handle node strings
    if (node1.string > NORMAL) {
        
        if (node2.string > NORMAL) {
            if (node1.string == node2.string ) {
                mod->grid->node[the_node].string = node1.string;
            } else {
                fprintf(stderr, "WARNING:  There are adjacent node strings!!\n");
                if (node1.string < node2.string)
                    mod->grid->node[the_node].string = node1.string;
                else
                    mod->grid->node[the_node].string = node2.string;
            }
        } else if (node2.string == NORMAL) {
            mod->grid->node[the_node].string  = NORMAL;
        } else {
            tl_error("Bad right parent node in node_new.");
        }
        
    } else if (node1.string == NORMAL) {
        
        if (node2.string > NORMAL) {
            mod->grid->node[the_node].string = NORMAL;
        }
        else if (node2.string == NORMAL) {
            mod->grid->node[the_node].string = NORMAL;
        }
        else {
            tl_error("Bad right parent node in node_new.");
        }
        
    } else {
        tl_error("Bad left parent node in node_new.");
    }
    
    
    // handle edge strings
    if (node1.edge_string > NORMAL) {
        
        if (node2.edge_string > NORMAL) {
            if (node1.edge_string == node2.edge_string) {
                mod->grid->node[the_node].edge_string = node1.edge_string;
            }
            else {
                fprintf(stderr, "WARNING:  There are adjacent node strings!!\n");
                if (node1.edge_string < node2.edge_string)
                    mod->grid->node[the_node].edge_string  = node1.edge_string;
                else
                    mod->grid->node[the_node].edge_string  = node2.edge_string;
            }
        }
        else if (node2.edge_string == NORMAL) {
            mod->grid->node[the_node].edge_string  = NORMAL;
        }
        else {
            tl_error("Bad right parent node in node_new.");
        }
        
    } else if (node1.edge_string== NORMAL) {
        
        if (node2.edge_string > NORMAL) {
            mod->grid->node[the_node].edge_string  = NORMAL;
        }
        else if (node2.edge_string == NORMAL) {
            mod->grid->node[the_node].edge_string  = NORMAL;
        }
        else {
            tl_error("Bad right parent node in node_new.");
        }
        
    } else {
        tl_error("Bad left parent node in node_new.");
    }
    
    mod->grid->node[the_node].id = the_node;
    
    /* the orig_nd_number is set to UNSET_INT */
    mod->grid->node[the_node].original_id = UNSET_INT;
    
    /* sets the node level */
    if (node1.level > node2.level)
        mod->grid->node[the_node].level = node1.level + 1;
    else
        mod->grid->node[the_node].level = node2.level + 1;
    
    /* sets the node block */
    mod->grid->node[the_node].block = UNSET_INT;
    
    /* determine what the new node is and which processor owns it */
#ifdef _MESSG
    /* set subdomain first */
    if (node1.resident_pe <= node2.resident_pe)
        mod->grid->node[the_node].resident_pe = node1.resident_pe;
    else
        mod->grid->node[the_node].resident_pe = node2.resident_pe;
    
    /* set rnode if I own the node, otherwise comm_update_edges will set rnode */
    if (mod->grid->node[the_node].resident_pe == mod->grid->smpi->myid){
        mod->grid->node[the_node].resident_id = the_node;
    }
    
#else
    mod->grid->my_nnodes++;
#endif
}
