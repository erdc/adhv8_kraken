#include "extrusion.h"

/* Create the 3D geometry
 *
 * Input:
 *   geo2d -- 2D data
 *   mat_layers -- the number of layers for each material
 *   head -- total head (from hot start file)
 *
 * Output:
 *   geo3d -- 3D data
 *   node_count -- the number of nodes in the vertical direction at a given node in
 *   the 2D mesh
 *   node_start -- the starting index of the node *below* the surface node.
 *   The nodes down the vertical line will be numbered sequentially,
 *   starting at this value
 *
 *   NEW: going to use this to map the 2D node numbers to renumber, so that node numbers
 *   go consecutively down columns... So that in the first vertical column there will be
 *   nodes 0, 1, ... n1  and the second column has nodes numbered n1 +1, n1 + 2, ...
 *
 *   NOTE:  We are keeping track of several different view of the geometry, so
 *   this can all get confusing.
 *   geo2d has the original 2D mesh
 *   geo3d has the extruded 3D mesh, with the nodes numbered to go consecutively
 *   from the surface node to the bottom.
 *
 *   To figure out the connection between geo2d and geo3d, it is useful to make
 *   use of node_start which maps the original node in geo2d to the corresponding
 *   (surface) node in the 3d mesh.
 *       node_start[my_node] takes the original node 'my_node' and will return
 *   the id of the node in the 3d mesh that corresponds to it.
 *
 *   node_count gives the number of nodes that are in the vertical line.
 *      node_count[my_node] will be the number of nodes that lie below the surface
 *   at that point.
 *
 */

//***********************************
// STANDARD PRISM LOCAL MAPPING *****
//
// (0) o---------o (2)
//     |\       /|
//     | \     / |
//     |  \   /  |
//     |   \ /   |
//     |    o(1) |
//     |    |    |
//     |    |    |
//     |    |    |
// (3) o----|----o (5)
//      \   |   /
//       \  |  /
//        \ | /
//         \|/
//          o(4)
//
//***********************************
//***********************************

static int DEBUG = 0;

//**************************************************************
//**************************************************************

int ccw(SVECT a, SVECT b, SVECT c) {
    double number = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y);
    if (number < 0) {
        return -1;
    } else if (number > 0) {
        return 1;
    } else {
        printf("ERROR :: ccw vectors are parallell\n");
        exit(-1);
    }
    return 0;
}

//**************************************************************
//**************************************************************

// NOTE: here I number nodes so that volume is always (+).
void add_triangular_prism(SGRID *geo3d, int *ndn, int ie, int mat, int *nelem3d) {
    
    int i,j;
    SVECT p[6];
    double x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6, z1, z2, z3, z4, z5, z6, d1, d2, d3, prism_term, t1, prism_volume = 0.;
    
    int eord[6] = {3,4,5,0,1,2};
    geo3d->elem3d[*nelem3d].nodes[0] = ndn[eord[0]];
    geo3d->elem3d[*nelem3d].nodes[1] = ndn[eord[1]];
    geo3d->elem3d[*nelem3d].nodes[2] = ndn[eord[2]];
    geo3d->elem3d[*nelem3d].nodes[3] = ndn[eord[3]];
    geo3d->elem3d[*nelem3d].nodes[4] = ndn[eord[4]];
    geo3d->elem3d[*nelem3d].nodes[5] = ndn[eord[5]];
    
    x1 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[0] ].x;
    x2 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[1] ].x;
    x3 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[2] ].x;
    x4 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[3] ].x;
    x5 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[4] ].x;
    x6 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[5] ].x;
    
    y1 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[0] ].y;
    y2 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[1] ].y;
    y3 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[2] ].y;
    y4 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[3] ].y;
    y5 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[4] ].y;
    y6 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[5] ].y;
    
    z1 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[0] ].z;
    z2 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[1] ].z;
    z3 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[2] ].z;
    z4 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[3] ].z;
    z5 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[4] ].z;
    z6 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[5] ].z;
    
    //printf("nodes :: %d %d %d %d %d %d \n",geo3d->elem3d[*nelem3d].nodes[0],geo3d->elem3d[*nelem3d].nodes[1],geo3d->elem3d[*nelem3d].nodes[2],geo3d->elem3d[*nelem3d].nodes[3],geo3d->elem3d[*nelem3d].nodes[4],geo3d->elem3d[*nelem3d].nodes[5]);
    //printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",x1,x2,x3,x4,x5,x6);
    //printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",y1,y2,y3,y4,y5,y6);
    //printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",z1,z2,z3,z4,z5,z6);
    
    // get point on surface/bottom triangle that is opposite the longest edge
    double dist12 = sqrt(pow(x1-x2,2) + pow(y1-y2,2));
    double dist23 = sqrt(pow(x2-x3,2) + pow(y2-y3,2));
    double dist31 = sqrt(pow(x3-x1,2) + pow(y3-y1,2));
    //if (dist12 > dist23 && dist12 > dist31 ) {
    if (fabs(dist12 - dist23)>1e-6 && fabs(dist12 - dist31)>1e-6) {
        p[0].x = x3; p[0].y = y3; p[0].z = z3; eord[0] = 2;
        // now move counter-clockwise
        p[1].x = x1; p[1].y = y1; p[1].z = z1; eord[1] = 0;
        p[2].x = x2; p[2].y = y2; p[2].z = z2; eord[2] = 1;
        if (ccw(p[0],p[1],p[2]) <= 0) {
            p[1].x = x2; p[1].y = y2; p[1].z = z2; eord[1] = 1;
            p[2].x = x1; p[2].y = y1; p[2].z = z1; eord[2] = 0;
            
            //printf("dist12: %20.10f \t dist23: %20.10f \t dist31: %20.10f\n",dist12,dist23,dist31);
            //printf("p[0].x: %10.5f \t p[0].y: %20.10f \t p[0].z: %20.10f\n",p[0].x,p[0].y,p[0].z);
            //printf("p[1].x: %10.5f \t p[1].y: %20.10f \t p[1].z: %20.10f\n",p[1].x,p[1].y,p[1].z);
            //printf("p[2].x: %10.5f \t p[2].y: %20.10f \t p[2].z: %20.10f\n",p[2].x,p[2].y,p[2].z);
            //printf("ccw1: %d\n",ccw(p[0],p[1],p[2]));
            
            assert(ccw(p[0],p[1],p[2]) > 0);
        }
    //} else if (dist23 > dist12 && dist23 > dist31 ) {
    } else if (fabs(dist23 - dist12)>1e-6 && fabs(dist23 - dist12)>1e-6) {
        p[0].x = x1; p[0].y = y1; p[0].z = z1; eord[0] = 0;
        // now move counter-clockwise
        p[1].x = x2; p[1].y = y2; p[1].z = z2; eord[1] = 1;
        p[2].x = x3; p[2].y = y3; p[2].z = z3; eord[2] = 2;
        if (ccw(p[0],p[1],p[2]) <= 0) {
            p[1].x = x3; p[1].y = y3; p[1].z = z3; eord[1] = 2;
            p[2].x = x2; p[2].y = y2; p[2].z = z2; eord[2] = 1;
            //printf("ccw2: %d\n",ccw(p[0],p[1],p[2]));
            assert(ccw(p[0],p[1],p[2]) > 0);
        }
    } else { // this could be equilateral, where it doesnt matter
        p[0].x = x2; p[0].y = y2; p[0].z = z2; eord[0] = 1;
        // now move counter-clockwise
        p[1].x = x1; p[1].y = y1; p[1].z = z1; eord[1] = 0;
        p[2].x = x3; p[2].y = y3; p[2].z = z3; eord[2] = 2;
        if (ccw(p[0],p[1],p[2]) <= 0) {
            p[1].x = x3; p[1].y = y3; p[1].z = z3; eord[1] = 2;
            p[2].x = x1; p[2].y = y1; p[2].z = z1; eord[2] = 0;
            //printf("p[0].x: %10.5f \t p[0].y: %20.10f \t p[0].z: %20.10f\n",p[0].x,p[0].y,p[0].z);
            //printf("p[1].x: %10.5f \t p[1].y: %20.10f \t p[1].z: %20.10f\n",p[1].x,p[1].y,p[1].z);
            //printf("p[2].x: %10.5f \t p[2].y: %20.10f \t p[2].z: %20.10f\n",p[2].x,p[2].y,p[2].z);
            //printf("ccw3: %d\n",ccw(p[0],p[1],p[2]));
            assert(ccw(p[0],p[1],p[2]) > 0);
        }
    }
    
    // find points above
    if (fabs(p[0].x - x4)<1e-6 && fabs(p[0].y - y4)<1e-6) {
        p[3].x = x4; p[3].y = y4; p[3].z = z4; eord[3] = 3;
    } else if (fabs(p[0].x - x5)<1e-6 && fabs(p[0].y - y5)<1e-6) {
        p[3].x = x5; p[3].y = y5; p[3].z = z5; eord[3] = 4;
    } else {
        p[3].x = x6; p[3].y = y6; p[3].z = z6; eord[3] = 5;
    }
    
    if (fabs(p[1].x - x4)<1e-6 && fabs(p[1].y - y4)<1e-6) {
        p[4].x = x4; p[4].y = y4; p[4].z = z4; eord[4] = 3;
    } else if (fabs(p[1].x - x5)<1e-6 && fabs(p[1].y - y5)<1e-6) {
        p[4].x = x5; p[4].y = y5; p[4].z = z5; eord[4] = 4;
    } else {
        p[4].x = x6; p[4].y = y6; p[4].z = z6; eord[4] = 5;
    }
    
    if (fabs(p[2].x - x4)<1e-6 && fabs(p[2].y - y4)<1e-6) {
        p[5].x = x4; p[5].y = y4; p[5].z = z4; eord[5] = 3;
    } else if (fabs(p[2].x - x5)<1e-6 && fabs(p[2].y - y5)<1e-6) {
        p[5].x = x5; p[5].y = y5; p[5].z = z5; eord[5] = 4;
    } else {
        p[5].x = x6; p[5].y = y6; p[5].z = z6; eord[5] = 5;
    }
    
    assert(p[3].z > p[0].z);
    assert(p[4].z > p[1].z);
    assert(p[5].z > p[2].z);
    
    // check if node numbering gives a prism volume that is (+)
    prism_volume = get_triprism_volume(p);
    
//    if(prism_volume <= 0) {
//        
//        eord[0] = 0;  eord[1] = 1;  eord[2] = 2;  eord[3] = 3;  eord[4] = 4;  eord[5] = 5;
//        geo3d->elem3d[*nelem3d].nodes[0] = ndn[eord[0]];
//        geo3d->elem3d[*nelem3d].nodes[1] = ndn[eord[1]];
//        geo3d->elem3d[*nelem3d].nodes[2] = ndn[eord[2]];
//        geo3d->elem3d[*nelem3d].nodes[3] = ndn[eord[3]];
//        geo3d->elem3d[*nelem3d].nodes[4] = ndn[eord[4]];
//        geo3d->elem3d[*nelem3d].nodes[5] = ndn[eord[5]];
//        
//        x1 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[0] ].x;
//        x2 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[1] ].x;
//        x3 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[2] ].x;
//        y1 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[0] ].y;
//        y2 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[1] ].y;
//        y3 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[2] ].y;
//        z1 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[0] ].z;
//        z2 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[1] ].z;
//        z3 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[2] ].z;
//        z4 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[3] ].z;
//        z5 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[4] ].z;
//        z6 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[5] ].z;
//        prism_volume = elem3d_get_prism_volume(x1,x2,x3,y1,y2,y3,z1,z2,z3,z4,z5,z6);
//    }
    
    
    if(prism_volume <= 0) {
        printf("ERROR: prism volume <= 0\n");
        printf("prism volume: %20.10f\n",prism_volume);
        printf("p values: \n");
        printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",p[0].x,p[1].x,p[2].x,p[3].x,p[4].x,p[5].x);
        printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",p[0].y,p[1].y,p[2].y,p[3].y,p[4].y,p[5].y);
        printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",p[0].z,p[1].z,p[2].z,p[3].z,p[4].z,p[5].z);
        exit(-1);
    }
    
    
    geo3d->elem3d[*nelem3d].nodes[0] = ndn[eord[3]];
    geo3d->elem3d[*nelem3d].nodes[1] = ndn[eord[4]];
    geo3d->elem3d[*nelem3d].nodes[2] = ndn[eord[5]];
    geo3d->elem3d[*nelem3d].nodes[3] = ndn[eord[0]];
    geo3d->elem3d[*nelem3d].nodes[4] = ndn[eord[1]];
    geo3d->elem3d[*nelem3d].nodes[5] = ndn[eord[2]];
    
    geo3d->elem3d[*nelem3d].mat = mat;
    geo3d->elem3d[*nelem3d].elem2d_sur = ie;
    geo3d->elem3d[*nelem3d].nnodes = NDONPRISM;
    
    x1 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[0] ].x;
    x2 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[1] ].x;
    x3 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[2] ].x;
    x4 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[3] ].x;
    x5 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[4] ].x;
    x6 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[5] ].x;
    
    y1 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[0] ].y;
    y2 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[1] ].y;
    y3 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[2] ].y;
    y4 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[3] ].y;
    y5 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[4] ].y;
    y6 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[5] ].y;
    
    z1 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[0] ].z;
    z2 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[1] ].z;
    z3 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[2] ].z;
    z4 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[3] ].z;
    z5 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[4] ].z;
    z6 = geo3d->node[ geo3d->elem3d[*nelem3d].nodes[5] ].z;
    
    //printf("DEBUG\n");
    //printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",p[0].x,p[1].x,p[2].x,p[3].x,p[4].x,p[5].x);
    //printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",p[0].y,p[1].y,p[2].y,p[3].y,p[4].y,p[5].y);
    //printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n\n",p[0].z,p[1].z,p[2].z,p[3].z,p[4].z,p[5].z);
    //printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",x1,x2,x3,x4,x5,x6);
    //printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",y1,y2,y3,y4,y5,y6);
    //printf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n\n",z1,z2,z3,z4,z5,z6);
    
    (*nelem3d)++;
}


//**************************************************************
//**************************************************************

void add_tetrahedron(SGRID *geo3d, int *ndn, int *eord, int reverse, int ie, int mat, int *nelem3d) {
    
    int i;
    int node_ids[NDONTET];
    for (i=0; i<NDONTET; i++) {node_ids[i] = ndn[eord[i]];}
    if (reverse) {
        node_ids[0] = ndn[eord[2]];
        node_ids[1] = ndn[eord[1]];
        node_ids[2] = ndn[eord[0]];
        node_ids[3] = ndn[eord[3]];
    }
    for (i=0; i<NDONTET; i++) {geo3d->elem3d[*nelem3d].nodes[i] = node_ids[i];}
    geo3d->elem3d[*nelem3d].mat = mat;
    geo3d->elem3d[*nelem3d].elem2d_sur = ie;
    geo3d->elem3d[*nelem3d].nnodes = NDONTET;
    
    // make sure nodes are properly numbered
    SNODE nodes[NDONTET];
    for (i=0; i<NDONTET; i++) {nodes[i] = geo3d->node[node_ids[i]];} // shallow copy
    if (get_tet_linear_djac(nodes,NULL) < 0.0) {
        fprintf(stderr, "Improperly numbered tetrahedron :: reverse :: ie2d: %d :: ie3d: %d  :: nodes :: %d, %d, %d, %d \n",
                ie,*nelem3d,geo3d->elem3d[*nelem3d].nodes[0],geo3d->elem3d[*nelem3d].nodes[1],geo3d->elem3d[*nelem3d].nodes[2],geo3d->elem3d[*nelem3d].nodes[3]);
        for (i=0; i<NDONTET; i++) {
            fprintf(stderr,"node %d {x,y,z} :: {%20.10f, %20.10f, %20.10f} \n",i,geo3d->node[geo3d->elem3d[(*nelem3d)].nodes[i]].x,geo3d->node[geo3d->elem3d[(*nelem3d)].nodes[i]].y,geo3d->node[geo3d->elem3d[(*nelem3d)].nodes[i]].z);
        }
        fflush(stdout);
        fflush(stderr);
        exit(0);
    }
    
    (*nelem3d)++;
}

//**************************************************************
//**************************************************************
void break_prism(SGRID *geo2d, SGRID *geo3d, int *ndn, int ie, int reverse, int *nelem3d) {
    
    int n = 0;
    int eord[3][4] = { {0, 1, 2, 5}, {0, 1, 5, 4}, {0, 4, 5, 3} };
    
    if (DEBUG) printf("breaking prism\n"); fflush(stdout);
    for (n=0; n<3; n++) {
        add_tetrahedron(geo3d, ndn, eord[n], reverse, ie, geo2d->elem2d[ie].mat, nelem3d);
    }
}

//**************************************************************
//**************************************************************
void break_pyramid(SGRID *geo2d, SGRID *geo3d, int *ndn, int ie, int reverse, int *nelem3d) {
    
    int bad_node = -3;
    int eord[4] = {0, 0, 0, 0};
    
    if (DEBUG) printf("breaking pyramid\n");
    fflush(stdout);
    fflush(stderr);
    
    // find which node in prism is missing
    int i;
    for (i=0; i<3; i++) {
        if (ndn[i+3] == UNSET_INT) bad_node = i+3;
    }
    
    if (bad_node == 3) {   // node 3 does not belong
        
        // first tet - two options are avail
        eord[0] = 0;
        eord[1] = 1;
        eord[2] = 2;
        eord[3] = 5;
        add_tetrahedron(geo3d, ndn, eord, reverse, ie, geo2d->elem2d[ie].mat, nelem3d);
        //printf("broke pyramid :: EC :: %d %d %d %d\n",geo3d->elem3d[*nelem3d-1].nodes[0],geo3d->elem3d[*nelem3d-1].nodes[1],geo3d->elem3d[*nelem3d-1].nodes[2],geo3d->elem3d[*nelem3d-1].nodes[3]);
        
        eord[0] = 0;
        eord[1] = 1;
        eord[2] = 5;
        eord[3] = 4;
        add_tetrahedron(geo3d, ndn, eord, reverse, ie, geo2d->elem2d[ie].mat, nelem3d);
        //printf("broke pyramid :: EC :: %d %d %d %d\n",geo3d->elem3d[*nelem3d-1].nodes[0],geo3d->elem3d[*nelem3d-1].nodes[1],geo3d->elem3d[*nelem3d-1].nodes[2],geo3d->elem3d[*nelem3d-1].nodes[3]);
    }
    else if (bad_node == 4) {  // node 4 does not belong
        
        // first tet - two options are avail
        eord[0] = 0;
        eord[1] = 1;
        eord[2] = 2;
        eord[3] = 5;
        add_tetrahedron(geo3d, ndn, eord, reverse, ie, geo2d->elem2d[ie].mat, nelem3d);
        //printf("broke pyramid :: EC :: %d %d %d %d\n",geo3d->elem3d[*nelem3d-1].nodes[0],geo3d->elem3d[*nelem3d-1].nodes[1],geo3d->elem3d[*nelem3d-1].nodes[2],geo3d->elem3d[*nelem3d-1].nodes[3]);
        
        eord[0] = 0;
        eord[1] = 1;
        eord[2] = 5;
        eord[3] = 3;
        add_tetrahedron(geo3d, ndn, eord, reverse, ie, geo2d->elem2d[ie].mat, nelem3d);
        //printf("broke pyramid :: EC :: %d %d %d %d\n",geo3d->elem3d[*nelem3d-1].nodes[0],geo3d->elem3d[*nelem3d-1].nodes[1],geo3d->elem3d[*nelem3d-1].nodes[2],geo3d->elem3d[*nelem3d-1].nodes[3]);
    }
    else if (bad_node == 5) {  // node 5 does not belong
        
        // first tet - two options are avail
        eord[0] = 0;
        eord[1] = 1;
        eord[2] = 2;
        eord[3] = 4;
        add_tetrahedron(geo3d, ndn, eord, reverse, ie, geo2d->elem2d[ie].mat, nelem3d);
        //printf("broke pyramid :: EC :: %d %d %d %d\n",geo3d->elem3d[*nelem3d-1].nodes[0],geo3d->elem3d[*nelem3d-1].nodes[1],geo3d->elem3d[*nelem3d-1].nodes[2],geo3d->elem3d[*nelem3d-1].nodes[3]);
        
        eord[0] = 0;
        eord[1] = 2;
        eord[2] = 3;
        eord[3] = 4;
        add_tetrahedron(geo3d, ndn, eord, reverse, ie, geo2d->elem2d[ie].mat, nelem3d);
        //printf("broke pyramid :: EC :: %d %d %d %d\n",geo3d->elem3d[*nelem3d-1].nodes[0],geo3d->elem3d[*nelem3d-1].nodes[1],geo3d->elem3d[*nelem3d-1].nodes[2],geo3d->elem3d[*nelem3d-1].nodes[3]);
    }
}

//**************************************************************
//**************************************************************
void break_tet(SGRID *geo2d, SGRID *geo3d, int *ndn, int ie, int reverse, int *nelem3d) {
    
    int good_node = -3;
    int eord[4] = {0, 0, 0, 0};
    
    if (DEBUG) printf("breaking tet\n");
    fflush(stdout);
    fflush(stderr);
    
    // find which node in prism is missing
    int i;
    for (i=0; i<3; i++) {
        if (ndn[i+3] != UNSET_INT) good_node = i+3;
    }
    
    eord[0] = 0;
    eord[1] = 1;
    eord[2] = 2;
    eord[3] = good_node;
    add_tetrahedron(geo3d, ndn, eord, reverse, ie, geo2d->elem2d[ie].mat, nelem3d);
    //printf("broke tet :: EC :: %d %d %d %d\n",geo3d->elem3d[*nelem3d-1].nodes[0],geo3d->elem3d[*nelem3d-1].nodes[1],geo3d->elem3d[*nelem3d-1].nodes[2],geo3d->elem3d[*nelem3d-1].nodes[3]);
}


//**************************************************************
//**************************************************************

int check_reverse(SGRID *geo3d, int *ndn) {
    SVECT vr, v1, v2, v3, vnode0, vnode1, vnode2, vnode3;
    
    vnode0.x = geo3d->node[ndn[0]].x; vnode0.y = geo3d->node[ndn[0]].y; vnode0.z = geo3d->node[ndn[0]].z;
    vnode1.x = geo3d->node[ndn[1]].x; vnode1.y = geo3d->node[ndn[1]].y; vnode1.z = geo3d->node[ndn[1]].z;
    vnode2.x = geo3d->node[ndn[2]].x; vnode2.y = geo3d->node[ndn[2]].y; vnode2.z = geo3d->node[ndn[2]].z;
    vnode3.x = geo3d->node[ndn[3]].x; vnode3.y = geo3d->node[ndn[3]].y; vnode3.z = geo3d->node[ndn[3]].z;
    
    v1 = svect_subtract(vnode1, vnode0);
    v2 = svect_subtract(vnode2, vnode0);
    v3 = svect_subtract(vnode3, vnode0);
    vr = svect_cross(v1,v2);
    if (svect_dotp(vr, v3) > 0.0) {
        return 0;
    } else {
        return 1;
    }
    printf("ERROR :: dot product in check reverse is 0\n");
    exit(-1);
    return -1;
}

//**************************************************************
//**************************************************************
//**************************************************************
//**************************************************************

void create_3d_geometry_adh(SGRID *geo2d, SGRID *geo3d, int *node_count, int *node_start, int *node_bin, double dz, double *head, zbins *bin, int bin_flag, int mesh_type) {
    
    int i, k, ll, n, nc, nn, ie, node_index, nt; /* counters and placeholders */
    int nelem3d = 0;
    int reverse;
    int ind[NDONTRI], iep[NDONTRI], node2d[2*NDONTRI], node3d[2*NDONTRI], node_layers[NDONTRI];
    int ilayer;
    double *surface;
    
    printf(">> creating 3d geometry\n");
    
    surface = (double *) tl_alloc(sizeof(double), geo2d->nnodes);
    
    /* Initialize 3D structure */
    geo3d->nnodes = 0;
    geo3d->nelems2d = 0;
    geo3d->nelems3d = 0;
    geo3d->nmat = geo2d->nmat;
    
    /* Set surface by adding elevation to head */
    for (i = 0; i < geo2d->nnodes; i++) {
        surface[i] = head[i] + geo2d->node[i].z;
    }
    
    /* Count up the number of nodes in the 3D mesh */
    geo3d->nnodes = 0;
    for (n = 0; n < geo2d->nnodes; n++) {
        geo3d->nnodes = geo3d->nnodes + node_count[n];
    }
    geo3d->node = (SNODE *) tl_alloc(sizeof(SNODE), geo3d->nnodes);
    
    /* Number the nodes */
    node_index = 0;
    for (nn = 0; nn < geo2d->nnodes; nn++) {
        node_start[nn] = node_index;
        geo3d->node[node_index].x = geo2d->node[nn].x;
        geo3d->node[node_index].y = geo2d->node[nn].y;
        geo3d->node[node_index].z = surface[nn];
        geo3d->node[node_index].global_surf_id = nn;
        
        node_index++;
        for (ll = 1; ll < node_count[nn]; ll++) {
            geo3d->node[node_index].x = geo3d->node[node_start[nn]].x;
            geo3d->node[node_index].y = geo3d->node[node_start[nn]].y;
            
            if (bin_flag) {
                //printf("node_bin[%d]: %d\n",nn,node_bin[nn]);
                dz = find_dz_bin(bin, geo3d->node[node_index-1].z, node_bin[nn]); // find bin
            }
            
            if (ll == node_count[nn]-1) { // bypass dz, just put node on bed
                geo3d->node[node_index].z = geo2d->node[nn].z;
                geo3d->node[node_index].global_surf_id = nn;
                node_index++;
                break;
            } else {
                geo3d->node[node_index].z = geo3d->node[node_index-1].z - dz;
                geo3d->node[node_index].global_surf_id = nn;
                node_index++;
            }
        }
    }
    
    
    if (DEBUG) {
        for (i=0; i<node_index; i++) {
            printf("node: %d {x,y,z} = {%10.5f, %10.5f, %10.5f}\n",i,geo3d->node[i].x,geo3d->node[i].y,geo3d->node[i].z);
        }
    }
    
    // sanity check :: make sure nodes are aligned in columns
    double xx,yy;
    node_index = 0;
    for (nn = 0; nn < geo2d->nnodes; nn++) {
        xx = geo3d->node[node_index].x;
        yy = geo3d->node[node_index].y;
        node_index++;
        for (ll = 1; ll < node_count[nn]; ll++) {
            if (fabs(xx - geo3d->node[node_index].x) > 1e-6 || fabs(yy - geo3d->node[node_index].y) > 1e-6) {
                printf("check 1 :: surface node: %d x,y = {%20.10f, %20.10f} :: subsurface node: %d x,y = {%20.10f, %20.10f}\n",node_index,xx,yy,node_index+1,geo3d->node[node_index].x,geo3d->node[node_index].y);
                fflush(stdout);
                fflush(stderr);
                exit(-1);
            }
            node_index++;
        }
    }
    
    /* Sanity check --- at this point, node_index should be the number of nodes in the 3D mesh */
    if (geo3d->nnodes != node_index) {
        printf("ERROR: geo3d->nnode != node_index\n");
        printf("geo3d->nnode = %d  node_index = %d\n", geo3d->nnodes, node_index);
        exit(0);
    }
    
    // count elements
    geo3d->nelems3d = 0;
    for (ie = 0; ie < geo2d->nelems2d; ie++) {
        //geo2d->elem2d[ie].elem2d_sur = ie;
        
        // Sort the indices by increasing node number
        for (n = 0; n < NDPRFC; n++) {
            ind[n] = 0;
            for (nc = 0; nc < NDPRFC; nc++) {
                if (node_start[geo2d->elem2d[ie].nodes[n]] > node_start[geo2d->elem2d[ie].nodes[nc]]) {
                    ind[n]++;
                }
            }
            iep[ind[n]] = n;
        }
        
        for (i=0; i<NDONTRI; i++) {
            node2d[i] = geo2d->elem2d[ie].nodes[iep[i]];
            node_layers[i] = node_count[node2d[i]] - 1;
        }
        
        // Traverse down the layers
        //geo3d->nelem3d = 0;
        for (ilayer = 1; ilayer<100000; ilayer++) {
            
            // count the number of horizonal nodes that can have this # of bed layers
            nt = 0;
            for (i=0; i<3; i++) {
                if (node_layers[i] >= ilayer) {
                    nt++;
                }
            }
            if (nt == 0) break;
            
            if (nt == 3) {
                // for a full layer (all three node columns have at least this # of layers), can be either a prism or a tet
                if (mesh_type != MIXED_ELEMENT_MESH) { // use triangular prism
                    geo3d->nelems3d += 3; // 3 tets in a full layer
                } else {
                    geo3d->nelems3d += 1; // 1 prism in a full layer
                }
            } else {
                // for a partial layer, can only use tets
                geo3d->nelems3d += nt;
            }
        }
    }
    
    selem3d_init_alloc_array(&(geo3d->elem3d), geo3d->nelems3d);
    for (i=0; i<geo3d->nelems3d; i++) {
        selem3d_alloc(&(geo3d->elem3d[i]), NDONPRISM);
    }
    
    //****************************************************************
    //----------------------------------------------------------------
    
    for (ie = 0; ie < geo2d->nelems2d; ie++) {
        //geo2d->elem2d[ie].column = ie;
        
        // Sort the indices by increasing node number
        for (n = 0; n < NDPRFC; n++) {
            ind[n] = 0;
            for (nc = 0; nc < NDPRFC; nc++) {
                if (node_start[geo2d->elem2d[ie].nodes[n]] > node_start[geo2d->elem2d[ie].nodes[nc]]) {
                    ind[n]++;
                }
            }
            iep[ind[n]] = n;
        }
        
        for (i=0; i<NDONTRI; i++) {
            node2d[i] = geo2d->elem2d[ie].nodes[iep[i]];
            node_layers[i] = node_count[node2d[i]] - 1;
        }
        
        
        // Traverse down the layers
        for (ilayer = 1; ilayer<100000; ilayer++) {
            
            //---------------------------------------------------
            // find nodes ---------------------------------------
            if (ilayer == 1) {
                nt = 3;
                for (i=0; i<3; i++) {
                    node3d[i] = node_start[node2d[i]];
                    node3d[i+3] = node3d[i] + 1;
                }
                reverse = check_reverse(geo3d, node3d); // just do at top
            } else {
                nt = 0;
                for (i=0; i<3; i++) {
                    // nodes3d already defined from last loop
                    if (node3d[i+3] != UNSET_INT) {
                        node3d[i] = node3d[i+3];
                    }
                    node3d[i+3] = UNSET_INT; // initially undefined
                    if (node_layers[i] >= ilayer) {
                        node3d[i+3] = node3d[i] + 1; // only increment if we can
                        nt++;
                    }
                }
            }
            if (nt == 0) break;
            
            //---------------------------------------------------
            //reverse = check_reverse(geo3d, node3d);
            //---------------------------------------------------
            
            if (nt == 3) { // break triangular prism into 3 tets (should always go here first!!! check for this later)
                if (mesh_type != MIXED_ELEMENT_MESH) {
                    break_prism(geo2d, geo3d, node3d, ie, reverse, &nelem3d);
                } else {
                    add_triangular_prism(geo3d, node3d, ie, geo2d->elem2d[ie].mat, &nelem3d);
                }
            } else if (nt == 2) { // break pyramid into 2 tets
                k=0;
                for (i=0; i<NDONPRISM; i++) {
                    if (node3d[i] != UNSET_INT) k++;
                }
                if (k != NDONPYRAMID) {
                    printf("k: %d \t NDONPYRAMID: %d\n",k,NDONPYRAMID);
                    exit(-1);
                }
                break_pyramid(geo2d, geo3d, node3d, ie, reverse, &nelem3d);
            } else if (nt == 1) { // tet is all that is left at bottom of grid
                k=0;
                for (i=0; i<NDONPRISM; i++) {
                    if (node3d[i] != UNSET_INT) k++;
                }
                if (k != NDONTET) {
                    printf("k: %d \t NDONPYRAMID: %d\n",k,NDONPYRAMID);
                    exit(-1);
                }
                break_tet(geo2d, geo3d, node3d, ie, reverse, &nelem3d);
            }
        }
    }
    
    if (geo3d->nelems3d != nelem3d) {
        printf("ERROR: CREATE 3D GEOMETRY :: geo3d->nelem3d: %d \t nelem3d: %d\n",geo3d->nelems3d,nelem3d);
        exit(-1);
    }

    int ntets=0,nprisms=0;
    for (i=0; i<nelem3d; i++) {
        if (geo3d->elem3d[i].nnodes == NDONTET) ntets++;
        if (geo3d->elem3d[i].nnodes == NDONPRISM) nprisms++;
    }
    
    //----------------------------------------------------------------
    //****************************************************************
    
    printf(">> ---- # of nodes: %d\n", geo3d->nnodes);
    printf(">> ---- # of 2D elems: %d\n", geo2d->nelems2d);
    printf(">> ---- # of 3D elems: %d\n", geo3d->nelems3d);
    printf(">> --------- # of tetrahedrons: %d\n",ntets);
    printf(">> --------- # of triangular prisms: %d\n",nprisms);
    printf(">> ---- # of materials: %d\n", geo3d->nmat);
    printf("*****************************************************\n\n");
    
    tl_free(sizeof(double), geo2d->nnodes, (void *) surface);
}
