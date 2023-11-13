#include "global_header.h"
// shake a point off a node by some small amount, but keep it in the element (ride an edge)
SVECT shakeAndBake(SGRID *grid, int ielem, int lnode) {
    int i,iedge,nd=UNSET_INT,nd1=UNSET_INT, nd2=UNSET_INT;
    double PERTURBATION = 1e-6;

    // find an edge with the node
    for (i=0; i<3; i++) {
        if (grid->nd_on_TriEdge[i][0] == lnode || grid->nd_on_TriEdge[i][1] == lnode) {
            iedge = i;
            break;
        }
    }

    nd  = grid->elem2d[ielem].nodes[lnode]; // global ID of node with point overlay
    nd1 = grid->elem2d[ielem].nodes[grid->nd_on_TriEdge[iedge][0]]; // global ID of node1 on edge
    nd2 = grid->elem2d[ielem].nodes[grid->nd_on_TriEdge[iedge][1]]; // global ID of node2 on edge

    assert(nd == nd1 || nd == nd2);

    SVECT p;
    if (nd == nd1) {
        if (grid->node[nd1].x < grid->node[nd2].x) {
            p.x = grid->node[nd1].x + PERTURBATION;
            p.y = interpolate1D(p.x, grid->node[nd1].x, grid->node[nd2].x, grid->node[nd1].y, grid->node[nd2].y);
        } else if (grid->node[nd1].x > grid->node[nd2].x) {
            p.x = grid->node[nd1].x - PERTURBATION;
            p.y = interpolate1D(p.x, grid->node[nd1].x, grid->node[nd2].x, grid->node[nd1].y, grid->node[nd2].y);
        } else { // vertical line
            p.x = grid->node[nd1].x ;
            if (grid->node[nd1].y < grid->node[nd2].y) {
                p.y = grid->node[nd1].y + PERTURBATION;
            } else {
                p.y = grid->node[nd1].y - PERTURBATION;
            }
        }
    } else {
        if (grid->node[nd1].x < grid->node[nd2].x) {
            p.x = grid->node[nd2].x - PERTURBATION;
            p.y = interpolate1D(p.x, grid->node[nd1].x, grid->node[nd2].x, grid->node[nd1].y, grid->node[nd2].y);
        } else if (grid->node[nd1].x > grid->node[nd2].x) {
            p.x = grid->node[nd2].x + PERTURBATION;
            p.y = interpolate1D(p.x, grid->node[nd1].x, grid->node[nd2].x, grid->node[nd1].y, grid->node[nd2].y);
        } else { // vertical line
            p.x = grid->node[nd1].x ;
            if (grid->node[nd1].y < grid->node[nd2].y) {
                p.y = grid->node[nd2].y - PERTURBATION;
            } else {
                p.y = grid->node[nd2].y + PERTURBATION;
            }
        }
    }


    // make sure point is still in triangle
    //double weights[4];
    //assert(pointInTriangle(grid,ielem,p,weights) == 1);


    return p;
}
