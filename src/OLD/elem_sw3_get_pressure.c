#include "global_header.h"

// These list have names that are misleading
// Really, they list the nodes that make up a segment on the parents element
// The first four indices are just placeholders
// ie.  {top,bot} {node 0, node 1} --> one segment where a midpt is located

//int top_node_list[NDONTETQUAD]={0,0,0,0,0,0,0,1,1,2};  // on parent tet
//int bot_node_list[NDONTETQUAD]={0,0,0,0,1,2,3,2,3,3};  // on parent tet

//*****************************************************************************************
//*****************************************************************************************

void elem3d_get_local_pressures(SGRID *grid, SSW_3D_ELEM *elem, double *prs) {
    int i, nd1, nd2;
    double prs_value[5] = {0., 0., 0., 0., 0.};
    double ret_val = 0.;
    int **edge;

    // zero presssure values
    sarray_init_dbl(elem->pressure, MAX_NNODES_ON_ELEM3D_QUAD);

    // set vertix unperturbed element pressures (needs this either unpertubed called, for for non-column perturbed call nodes
    global_to_local_dbl(prs, elem->pressure, elem->elem3d->nodes, elem->elem3d->nnodes);
    //ELEM3D_GET_LOCAL(prs, elem->pressure, elem->elem3d->nodes);

// CJT :: may not be anything wrong with this
//    // get midpoint unperturbed pressures
//    int nnodes_quad = NDONTETQUAD;
//    edge = grid->nd_on_TetEdge;
//    if (elem->elem3d->nnodes == NDONPRISM) {
//        nnodes_quad = NDONPRISMQUAD;
//        edge = grid->nd_on_PrismEdge;
//    }
//    for(i=elem->elem3d->nnodes; i<nnodes_quad; i++) {
//            //nd1 = top_node_list[i];  nd2 = bot_node_list[i];
//            nd1 = edge[i - elem->elem3d->nnodes][0];
//            nd2 = edge[i - elem->elem3d->nnodes][1];
//            ret_val = get_value_midpt_list(grid, grid->midpt_list, elem->elem3d->nodes[nd1], elem->elem3d->nodes[nd2], prs_value);
//            elem->pressure[i] = prs_value[0];
//    }
    
    
    if (elem->elem3d->nnodes == NDONPRISM) {
        for(i=NDONPRISM; i<NDONPRISMQUAD; i++) {
            nd1 = grid->nd_on_PrismEdge[i - elem->elem3d->nnodes][0];
            nd2 = grid->nd_on_PrismEdge[i - elem->elem3d->nnodes][1];
            ret_val = get_value_midpt_list(grid, grid->midpt_list, elem->elem3d->nodes[nd1], elem->elem3d->nodes[nd2], prs_value);
            elem->pressure[i] = prs_value[0];
        }
    } else {
        for(i=NDONTET; i<NDONTETQUAD; i++) {
            nd1 = grid->nd_on_TetEdge[i - elem->elem3d->nnodes][0];
            nd2 = grid->nd_on_TetEdge[i - elem->elem3d->nnodes][1];
            ret_val = get_value_midpt_list(grid, grid->midpt_list, elem->elem3d->nodes[nd1], elem->elem3d->nodes[nd2], prs_value);
            elem->pressure[i] = prs_value[0];
        }
    }
    
}

//*****************************************************************************************
//*****************************************************************************************

void elem2d_get_local_pressures(SGRID *grid, SSW_3D_ELEM *elem, double *prs) {
    int i, k, target_node, nd1, nd2;
    double prs_value[5] = {0., 0., 0., 0., 0.};
    double ret_val = 0.;
    int **edge;

    // zero presssure values
    sarray_init_dbl(elem->pressure, MAX_NNODES_ON_ELEM3D_QUAD);

    // set vertix unperturbed element pressures (needs this either unpertubed called, for for non-column perturbed call nodes
    global_to_local_dbl(prs, elem->pressure, elem->elem2d->nodes, elem->elem2d->nnodes);
    //ELEM2D_GET_LOCAL(prs, elem->pressure, elem->elem2d->nodes);

    // get midpoint unperturbed pressures
    int nnodes_quad = 6; // number of quadratic vertices on triangle
    edge = grid->nd_on_TriEdge;
    if (elem->elem2d->nnodes == NDONQUAD) {
        nnodes_quad = 8;
        edge = grid->nd_on_QuadEdge;
    }
    
    for(i=elem->elem2d->nnodes; i<nnodes_quad; i++) {
            //nd1 = top_node_list[i];  nd2 = bot_node_list[i];
            nd1 = edge[i - elem->elem2d->nnodes][0];
            nd2 = edge[i - elem->elem2d->nnodes][1];
            ret_val = get_value_midpt_list(grid, grid->midpt_list, elem->elem2d->nodes[nd1], elem->elem2d->nodes[nd2], prs_value);
            elem->pressure[i] = prs_value[0];
    } 

/*
    for(i = 0; i < NDONTRI; i++) {
            k = i + 1;
            if (k == NDONTRI) k = 0;
            ret_val = get_value_midpt_list(grid, grid->midpt_list, elem->elem2d->nodes[i], elem->elem2d->nodes[k], prs_value);
            
            target_node = i + NDONTRI;
            elem->pressure[target_node] = prs_value[0];
    }
*/
}

//*****************************************************************************************
//*****************************************************************************************
//////////////////////////////////////////////////////////////////////////////////////////
// This files returns all perturbed pressure values over a tet column element 
// input ::
//              prs_flag :: 1 if this node is in a column that is pertubed, 0 if not
//                   prs :: the unperturbed vertex pressure values
//                 prs_p :: the perturbed vertex pressure values
//          surface_node :: the surface node at the top of a given nodes column
//                  flag :: denotes + (+1) or - (-1) perturbation
// output ::
//              pressure :: contains perturbed pressures
///////////////////////////////////////////////////////////////////////////////////////////

void elem3d_get_perturbed_pressures(int *node_in_column_flag, SGRID *grid, SSW_3D_ELEM *elem, int flag) {

    int i,j;
    int valueID1, valueID2, nd1, nd2;
    double prs_value[5] = {0., 0., 0., 0., 0.};
    double ret_val = 0.;
    int **edge;
    MIDPT_LIST_ITEM *ptr;
    
//    printf("elem3d_get_perturbed_pressures :: mdpt pressures :: total number of midpoints: %d\n",grid->num_midpts);
//    int k=0;
//    for (i=0; i<grid->num_midpts; i++) {
//        ptr = grid->midpt_list[i];
//        while (ptr->next != NULL) {
//            k++;
//            printf("surface index: %d || midpoint %d || node1: %d \t node2: %d  ",i,k,ptr->node1+1,ptr->node2+1);
//            for (j=0; j<5; j++) {
//                printf("%20.10e",ptr->value[j]);
//            }
//            printf("\n");
//            ptr = ptr->next;
//        }
//    }

    if (flag != -1 && flag != 1) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> perturbation flag must be -1 or 1.");
    }

    int nnodes_quad = NDONTETQUAD;
    edge = grid->nd_on_TetEdge;
    if (elem->elem3d->nnodes == NDONPRISM) {
        nnodes_quad = NDONPRISMQUAD;
        edge = grid->nd_on_PrismEdge;
    }
    
    /* Perturb the pressure at a vertex node */
    for(j=0; j<elem->elem3d->nnodes; j++) {
        if(node_in_column_flag[j] == 1) {
            if (flag == 1) {         // + perturbation
                elem->pressure[j] = elem->prs_plus[j];
            } else if (flag == -1) { // - perturbation
                elem->pressure[j] = elem->prs_minus[j];
            }
        }
    }

    /* fill in the pressure values at all midpoints that lie under the perturbed node */
    for(j=elem->elem3d->nnodes; j<nnodes_quad; j++) {
        //nd1 = top_node_list[i];  nd2 = bot_node_list[i];
        nd1 = edge[j - elem->elem3d->nnodes][0];
        nd2 = edge[j - elem->elem3d->nnodes][1];
        
        ret_val = get_value_midpt_list(grid, grid->midpt_list, elem->elem3d->nodes[nd1], elem->elem3d->nodes[nd2], prs_value);
        
        //printf("local midpoint node: %d [%d:%d] :: values: %20.10f %20.10f %20.10f %20.10f %20.10f \n",j+1,nd1+1,nd2+1,prs_value[0],prs_value[1],prs_value[2],prs_value[3],prs_value[4]);

        if (flag == 1) {
            valueID1 = 1; // + pertubation of left node ID
            valueID2 = 2; // + pertubation of right node ID
        } else {
            valueID1 = 3; // - perturbation of left node ID
            valueID2 = 4; // - perturbation of right node ID
        }

        // cjt :: I have no idea what this is doing ...
        if(elem->surface_nodeID[nd1] < elem->surface_nodeID[nd2]) {
            if(node_in_column_flag[nd1] == 1) {
                elem->pressure[j] = prs_value[valueID1];
            }
            else if (node_in_column_flag[nd2] == 1) {
                elem->pressure[j] = prs_value[valueID2];
            }
            else {
                elem->pressure[j] = prs_value[0];
            }
        } else { // also handles the case where surface nodes are equal (column) :: make sure this is good in that case
            if(node_in_column_flag[nd1] == 1) {
                elem->pressure[j] = prs_value[valueID2];
            }
            else if(node_in_column_flag[nd2] == 1) {
                elem->pressure[j] = prs_value[valueID1];
            }
            else {
                elem->pressure[j] = prs_value[0];
            }
        }
        assert(elem->pressure[j] > -1e-6);
    }
}

//*****************************************************************************************
//*****************************************************************************************

void elem2d_get_perturbed_pressures(int *node_in_column_flag, SGRID *grid, SSW_3D_ELEM *elem, int flag) {

    int i,j,k;
    int valueID1, valueID2, nd1, nd2, target_node = UNSET_INT;
    double prs_value[5] = {0., 0., 0., 0., 0.};
    double ret_val = 0.;
    int **edge;

    if (flag != -1 && flag != 1) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> perturbation flag must be -1 or 1.");
    }

    int nnodes_quad = 6; // number of quadratic vertices on triangle
    edge = grid->nd_on_TriEdge;
    if (elem->elem2d->nnodes == NDONQUAD) {
        nnodes_quad = 8;
        edge = grid->nd_on_QuadEdge;
    }
    
    /* Perturb the pressure at a vertex node */
    for(j=0; j<elem->elem2d->nnodes; j++) {
        if(node_in_column_flag[j] == 1) {
            if (flag == 1) {         // + perturbation
                elem->pressure[j] = elem->prs_plus[j];
            } else if (flag == -1) { // - perturbation
                elem->pressure[j] = elem->prs_minus[j];
            }
        }
    }

    /* fill in the pressure values at all midpoints that lie under the perturbed node */

    for(j=elem->elem2d->nnodes; j<nnodes_quad; j++) {
        nd1 = edge[j - elem->elem2d->nnodes][0];
        nd2 = edge[j - elem->elem2d->nnodes][1];
        
        ret_val = get_value_midpt_list(grid, grid->midpt_list, elem->elem2d->nodes[nd1], elem->elem2d->nodes[nd2], prs_value);
    
        if (flag == 1) {
            valueID1 = 1;
            valueID2 = 2;
        } else {
            valueID1 = 3;
            valueID2 = 4;
        }

        if(elem->surface_nodeID[nd1] < elem->surface_nodeID[nd2]) {
            if(node_in_column_flag[nd1] == 1) {
                elem->pressure[j] = prs_value[valueID1];
            } else if (node_in_column_flag[nd2] == 1) {
                elem->pressure[j] = prs_value[valueID2];
            } else {
                elem->pressure[j] = prs_value[0];
            }
        } else {
            if(node_in_column_flag[nd1] == 1) {
                elem->pressure[j] = prs_value[valueID2];
            } else if(node_in_column_flag[nd2] == 1) {
                elem->pressure[j] = prs_value[valueID1];
            } else {
                elem->pressure[j] = prs_value[0];
            }
        }
    }
}
