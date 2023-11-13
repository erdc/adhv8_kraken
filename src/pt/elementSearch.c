#include "global_header.h"

int check_for_external_boundary(SGRID *grid, SVECT p1, SVECT p2, int *ielem, int lnd, int ledge, int lface, int DEBUG_SEARCH, SVECT *pi);

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

int findElementsAroundNode(SGRID *grid, SVECT p1, SVECT p2, int *ielem, int DEBUG_SEARCH, int nd, int elem_org, SVECT *pi,
                            int (*elementSearch)(SGRID *, SVECT, SVECT, int *i, bool, SVECT *)) {
    int i, flag = UNSET_INT, elem_store = *ielem;
    // Search all elements connected to this node for an intersect from p1 to p2
    //printf("grid->nc_nelems[%d]: %d\n",nd,grid->nc_nelems[nd]);
    if (grid->nc_nelems[nd] > 1) {
        for (i=0; i<grid->nc_nelems[nd]; i++) {
            if (grid->nc_elems[nd][i] == elem_org) continue;
            *ielem = grid->nc_elems[nd][i];
            if (DEBUG_SEARCH) printf("------ searching surrounding node %d connected element %d || elem_org: %d\n",nd,*ielem,elem_org);
            //exit(-1);
            flag = elementSearch(grid,p1,p2,ielem,false,pi); // non-iterative search
            //printf("flag: %d\n",flag);
            //exit(-1);
            if (flag !=0 ) break; // intersection found
        }
    } else {
        // if the intersect is through a node, and there are no more connected elements, it must be a boundary node
        *pi = p1;
        return -2;
    }
    //printf("findElementsAroundNode || after node loop || flag: %d\n",flag); exit(-1);
    if (flag == -2) {
       // *pi = p1;
        return -2; // hit external boundary
    }
    if (flag != 0) {
        // element connected to this node has an intersect, now carry on with iterative searching
        if (DEBUG_SEARCH) printf("findElementsAroundNode found node intersection\n");
        return elementSearch(grid,p1,p2,ielem,true,pi);
    }
    
    // if we're here, we did not find an element. If the node is an external node, this means the particle left the domain
    if (grid->node_flag[nd] == 1) {
        if (DEBUG_SEARCH) printf("particle wants to leave the domain \n");
        *ielem = elem_store;
        *pi = p1;
        return -2;
    }
    
    fflush(stdout);
    printf("ERROR :: FILE: %s LINE: %d :: flag: %d || ielem: %d || p1: %f %f %f || p2: %f %f %f || pi: %f %f %f || elem_org: %d \n",
           __FILE__,__LINE__,flag,*ielem,p1.x,p1.y,p1.z,p2.x,p2.y,p2.x,pi->x,pi->y,pi->z,elem_org);
    exit(-1);
    return 0;
}

int findElementsAroundEdge(SGRID *grid, SVECT p1, SVECT p2, int *ielem, int DEBUG_SEARCH, int iedge, int elem_org, SVECT *pi,
                            int (*elementSearch)(SGRID *, SVECT, SVECT, int *i, bool, SVECT *)) {
    assert(*ielem > -1 && *ielem < grid->nelems3d);
    int i, flag = UNSET_INT, elem_store = *ielem;
    // Search all elements connected to this node for an intersect from p1 to p2
    //printf("grid->nc_nelems[%d]: %d\n",nd,grid->nc_nelems[nd]);
    if (grid->elem3d[elem_store].nedge2elem[iedge] > 1) {
        for (i=0; i<grid->elem3d[elem_store].nedge2elem[iedge]; i++) {
            if (grid->elem3d[elem_store].edge2elem[iedge][i] == elem_org) continue;
            *ielem = grid->elem3d[elem_store].edge2elem[iedge][i];
            assert(*ielem > -1 && *ielem < grid->nelems3d);
            if (DEBUG_SEARCH) printf("------ searching surrounding edge %d connected element %d\n",iedge,*ielem);
            flag = elementSearch(grid,p1,p2,ielem,false,pi); // non-iterative search
            if (flag !=0 ) break; // intersection found
        }
    } else {
        // if the intersect is through a node, and there are no more connected elements, it must be a boundary node
        *pi = p1;
        return -2;
    }
    //printf("findElementsAroundNode || after node loop || flag: %d\n",flag);
    if (flag == -2) {
        //*pi = p1;
        return -2; // hit external boundary
    }
    if (flag != 0) {
        // element connected to this node has an intersect, now carry on with iterative searching
        return elementSearch(grid,p1,p2,ielem,true,pi);
    }
    
    // if we're here, we did not find an element. If the edge is external, this means the particle left the domain
    int nd1 = grid->elem3d[elem_store].nodes[ grid->nd_on_TetEdge[iedge][0] ];
    int nd2 = grid->elem3d[elem_store].nodes[ grid->nd_on_TetEdge[iedge][1] ];
    if (grid->node_flag[nd1] == 1 && grid->node_flag[nd2] == 1) {
        if (DEBUG_SEARCH) printf("particle wants to leave the domain through node \n");
        *ielem = elem_store;
        *pi = p1;
        return -2;
    }
    
    return 0;
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
// move all nodes to global IDS
int nodeLoopCheck2D(SGRID *grid, SVECT p1, SVECT p2, int elem_org, int DEBUG_SEARCH, int *ielem, SVECT *pi) {
    int i,ndt;
    bool isPointOnLine = false;
    
    for (i=0; i<3; i++) { // loop over nodes
        ndt = grid->elem2d[*ielem].nodes[i];
        isPointOnLine = pointOnLine2D(p1,p2,grid->node[ndt]);
        if (isPointOnLine) {
            //if (DEBUG_SEARCH) printf("------ point is on line \n");
            if (fabs(grid->node[ndt].x - p1.x) + fabs(grid->node[ndt].y- p1.y) + fabs(grid->node[ndt].z - p1.z) < 1e-12) continue;
            
            
//            if (fabs(grid->node[ndt].x/grid->xL - p1.x/grid->xL) + fabs(grid->node[ndt].y/grid->yL - p1.y/grid->yL)  < 1e-12) {
//                if (DEBUG_SEARCH) printf("------ skipping node %d\n",ndt);
//                continue;
//            }
            
            // is the intersected node on an external boundary?
//            if (grid->node_flag[ndt] == 1) {
//                if (DEBUG_SEARCH) printf("and intersects external boundary node %d \n",ndt);
//                *pi = grid->node[ndt];
//                return -2;
//            }
            // if not propagate the trajectory further
            return findElementsAroundNode(grid,grid->node[ndt],p2,ielem,DEBUG_SEARCH,ndt,elem_org,pi,elementSearch2D);
        }
    }
    return 0; // intersection point not found on element node
}

int nodeLoopCheck3D(SGRID *grid, SVECT p1, SVECT p2, int elem_org, int DEBUG_SEARCH, int *ielem, SVECT *pi) {
    int i,ndt;
    bool isPointOnLine = false;
    
    for (i=0; i<4; i++) { // loop over nodes
        ndt = grid->elem3d[*ielem].nodes[i];
        isPointOnLine = pointOnLine(p1,p2,grid->node[ndt]);
        if (isPointOnLine) {
            if (fabs(grid->node[ndt].x - p1.x) + fabs(grid->node[ndt].y - p1.y) + fabs(grid->node[ndt].z - p1.z) < 1e-12) continue;
            
            // is the intersected node on an external boundary?
//            if (grid->node_flag[ndt] == 1) {
//                if (DEBUG_SEARCH) printf("and intersects external boundary node %d \n",ndt);
//                *pi = grid->node[ndt];
//                return -2;
//            }
            // if not propagate the trajectory further
            return findElementsAroundNode(grid,grid->node[ndt],p2,ielem,DEBUG_SEARCH,ndt,elem_org,pi,elementSearch3D);
        }
    }
    return 0; // intersection point not found on element node
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

int edgeLoopCheck2D(SGRID *grid, SVECT p1, SVECT p2, int elem_org, int DEBUG_SEARCH, int *ielem, SVECT *pi) {
    int i,nd1_ID,nd2_ID, intersect_flag = UNSET_INT, elem_store = *ielem;
    SVECT nd1, nd2;
    
    for (i=0; i<3; i++) { // loop over edges
        nd1_ID = grid->elem2d[*ielem].nodes[grid->nd_on_TriEdge[i][0]];
        nd2_ID = grid->elem2d[*ielem].nodes[grid->nd_on_TriEdge[i][1]];
        nd1 = grid->node[nd1_ID];
        nd2 = grid->node[nd2_ID];
        intersect_flag = intersect2D_SegSeg(p1,p2,nd1,nd2,pi);
        //printf("\nedge nodes: %d %d || intersect_flag: %d || p1: %f %f %f || p2: %f %f %f || nd1: %f %f %f || nd2: %f %f %f || pi: %f %f %f\n",nd1_ID,nd2_ID,intersect_flag,p1.x,p1.y,p1.z,p2.x,p2.y,p2.z,nd1.x,nd1.y,nd1.z,nd2.x,nd2.y,nd2.z,pi->x,pi->y,pi->z);
        if (fabs(pi->x - p1.x) + fabs(pi->y - p1.y) + fabs(pi->z - p1.z) < 1e-12) continue;
        // segments are colinear along external boundary edge
        if (intersect_flag == 1) {
            // is the intersected edge on an external boundary?
            if (grid->elem2d[*ielem].edge_flag[i] == UNSET_INT) {
                //if (intersect_flag == -1) {
                //    *pi = p2;
                //}
                return -2;
            }
            // if not propagate the trajectory further
            *ielem = grid->elem2d[*ielem].edge_flag[i];
            if (DEBUG_SEARCH) printf(" and ends/travels through adjacent element %d (edgeloop) \n",*ielem);
            return elementSearch2D(grid,*pi,p2,ielem,true,pi);
        }
    }
    *ielem = elem_store;
    
    return 0; // intersection point not found on element edge
}

int edgeLoopCheck3D(SGRID *grid, SVECT p1, SVECT p2, int elem_org, int DEBUG_SEARCH, int *ielem, SVECT *pi) {
    int i,nd1_ID,nd2_ID, intersect_flag = UNSET_INT;
    SVECT nd1, nd2;
    
    for (i=0; i<6; i++) { // loop over edges
        nd1_ID = grid->elem3d[*ielem].nodes[grid->nd_on_TetEdge[i][0]];
        nd2_ID = grid->elem3d[*ielem].nodes[grid->nd_on_TetEdge[i][1]];
        nd1 = grid->node[nd1_ID];
        nd2 = grid->node[nd2_ID];
        intersect_flag = intersect3D_SegSeg(p1,p2,nd1,nd2,pi);
        if (fabs(pi->x - p1.x) + fabs(pi->y - p1.y) + fabs(pi->z - p1.z) < 1e-12) continue;
        if (intersect_flag == 1) {
            // is the intersected edge on an external boundary?
            if (grid->elem3d[*ielem].edge_flag[i] == 1) return -2;
            // if not propagate the trajectory further
            if (DEBUG_SEARCH) printf(" and ends/travels through adjacent element %d (edgeloop) \n",*ielem);
            return findElementsAroundEdge(grid,*pi,p2,ielem,DEBUG_SEARCH,i,elem_org,pi,elementSearch3D);
            //return findElementsAroundNode(grid,*pi,p2,ielem,DEBUG_SEARCH,nd1_ID,elem_org,pi,elementSearch3D);
        }
    }
    return 0; // intersection point not found on element edge
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

int faceLoopCheck(SGRID *grid, SVECT p1, SVECT p2, int elem_org, int DEBUG_SEARCH, int *ielem, SVECT *pi, int check) {
    int i,j,k,nd1_ID,nd2_ID,nd3_ID,intersect_flag = UNSET_INT;
    SVECT nd1, nd2, nd3, pi_org = *pi;
    bool isPointOnLine = false;
    
    for (i=0; i<4; i++) { // loop over faces
        //printf("FACE %d\n",i);
        nd1_ID = grid->elem3d[*ielem].nodes[nd_on_fc[i][0]];
        nd2_ID = grid->elem3d[*ielem].nodes[nd_on_fc[i][1]];
        nd3_ID = grid->elem3d[*ielem].nodes[nd_on_fc[i][2]];
        nd1 = grid->node[nd1_ID];
        nd2 = grid->node[nd2_ID];
        nd3 = grid->node[nd3_ID];
        intersect_flag = intersect3D_SegTriangle(p1,p2,nd1,nd2,nd3,pi);
        if (fabs(pi->x - p1.x) + fabs(pi->y - p1.y) + fabs(pi->z - p1.z) < 1e-12) continue;
        if (intersect_flag == true) {
            //printf("in faceLoopCheck :: found intersection: intersect_flag: %d :: face_flag: %d :: pi: {%20.15e %20.15e %20.15e}\n",intersect_flag,grid->elem3d[*ielem].face_flag[i],pi->x,pi->y,pi->z);
            // did the particle trajectory intersect an external face?
            if (grid->elem3d[*ielem].face_flag[i] == UNSET_INT) {
                return -2;
            }
            // if not propagate the trajectory further
            *ielem = grid->elem3d[*ielem].face_flag[i];
            return elementSearch3D(grid,*pi,p2,ielem,true,pi);
        }
    }
    
    return 0; // intersection point not found on element face
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
// returns:
//              0 - particle trajectory never crosses element
//             -2 - particle trajectory crosses an external boundary node or edge
//             -1 = particle trajectory ends on element
int elementSearch2D(SGRID *grid, SVECT p1, SVECT p2, int *ielem, bool cycle, SVECT *pi) {
    
    //RUN_VERBOSE = true;
    
    int DEBUG_SEARCH = 0;
    if (DEBUG_SEARCH) {
        printf("\n---- ELEMENT SEARCH || p1 elem: %d || p1: {%20.15e,%20.15e,%20.15e) || p2: {%20.15e,%20.15e,%20.15e}\n",
               *ielem,p1.x,p1.y,p1.z,p2.x,p2.y,p2.z);
        int i;
        for (i=0; i<3; i++) {
            printf("---- ELEMENT SEARCH || element node %d || gid: %d || node_flag: %d || {x,y,z} || {%20.15e %20.15e %20.15e}\n",
                   i,grid->elem2d[*ielem].nodes[i],grid->node_flag[grid->elem2d[*ielem].nodes[i]],
                   grid->node[grid->elem2d[*ielem].nodes[i]].x,
                   grid->node[grid->elem2d[*ielem].nodes[i]].y,
                   grid->node[grid->elem2d[*ielem].nodes[i]].z);
        }
    }
    
    int elem_org = *ielem;
    double weights[3] = {0., 0., 0.};
    int check = UNSET_INT, flag = UNSET_INT, lnd = UNSET_INT, ledge = UNSET_INT, lface = UNSET_INT;
    
//    // sanity check, p1 better start in ielem
//    check = UNSET_INT; flag = UNSET_INT; lnd = UNSET_INT; ledge = UNSET_INT; lface = UNSET_INT;
//    flag = compute_interpolate_weights_2D_triangle(grid, *ielem, p1.x, p1.y, weights, &lnd, &ledge);
//    assert(flag==1);
    
    //--------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------
    // check if particle trajectory ends (p2) on element
    check = UNSET_INT; flag = UNSET_INT; lnd = UNSET_INT; ledge = UNSET_INT; lface = UNSET_INT;
    flag = compute_interpolate_weights_2D_triangle(grid, *ielem, p2.x, p2.y, weights, &lnd, &ledge);
    if(flag == 1) { // particle final position found in element
        //printf("------ checking if p2 falls on external boundary of original element || returns: %d\n",
        //       check_for_external_boundary(grid,p1,p2,ielem,lnd,ledge,lface,DEBUG_SEARCH,pi));
        return check_for_external_boundary(grid,p1,p2,ielem,lnd,ledge,lface,DEBUG_SEARCH,pi);
    }
    
    //--------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------
    // make sure particle trajectory starts (p1) on original element, and flag if it starts on node, edge, etc.
    check = UNSET_INT; flag = UNSET_INT; lnd = UNSET_INT; ledge = UNSET_INT; lface = UNSET_INT;
    flag = compute_interpolate_weights_2D_triangle(grid, *ielem, p1.x, p1.y, weights, &lnd, &ledge);
    if(flag != 1) {
    if (DEBUG_SEARCH) printf("-- trajectory does not start (p1) on original element \n");
        return 0; // can happen when trajectory was on edge and node connectivity is used
    }
    
    //--------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------
    int intersect_flag = UNSET_INT, nd, nd_ID, nd1_ID, nd2_ID;
    
    //--------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------
    // calculate intersection point
    //--------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------
    if (lnd != UNSET_INT) {
        //--------------------------------------------------------------------------------------
        //-- case 1 :: particle starts on node
        //--------------------------------------------------------------------------------------
        nd_ID = grid->elem2d[*ielem].nodes[lnd]; // starting node global ID
        if (DEBUG_SEARCH) printf("-- starts on node %d \n",nd_ID);
        
        // case 1-a :: intersection falls on another node
        if (DEBUG_SEARCH) printf("---- going to nodeLoopCheck2D || ielem: %d || p2: %f %f %f || pi: %f %f %f\n",*ielem,p2.x,p2.y,p2.z,pi->x,pi->y,pi->z);
                    //if (cycle == false) exit(-1);
        intersect_flag = nodeLoopCheck2D(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi);
                    //if (cycle == false) exit(-1);
        if (DEBUG_SEARCH) printf("---- finished nodeLoopCheck2D || ielem: %d || p2: %f %f %f || pi: %f %f %f\n",*ielem,p2.x,p2.y,p2.z,pi->x,pi->y,pi->z);
        if (intersect_flag != 0) return intersect_flag;


        // case 1-b :: intersection falls on an edge (only opposite edge is possible)
        if (DEBUG_SEARCH) printf("---- going to edgeLoopCheck2D || ielem: %d || p2: %f %f %f || pi: %f %f %f\n",*ielem,p2.x,p2.y,p2.z,pi->x,pi->y,pi->z);
        intersect_flag = edgeLoopCheck2D(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi);
        if (DEBUG_SEARCH) printf("---- finished edgeLoopCheck2D || ielem: %d || p2: %f %f %f || pi: %f %f %f\n",*ielem,p2.x,p2.y,p2.z,pi->x,pi->y,pi->z);
        if (intersect_flag != 0) return intersect_flag;
        
        // case 1-c :: intersection *NOT* found on this element
        if (cycle == false) return 0; // do not keep iterating
        // look at elements connected to this node
        if (DEBUG_SEARCH) printf("-- ends/travels through adjacent element %d \n",*ielem);
        return findElementsAroundNode(grid,p1,p2,ielem,DEBUG_SEARCH,nd_ID,elem_org,pi,elementSearch2D);
        
    } else if (ledge != UNSET_INT) {
        //--------------------------------------------------------------------------------------
        //-- case 1 :: particle starts on edge
        //--------------------------------------------------------------------------------------
        nd1_ID = grid->elem2d[*ielem].nodes[grid->nd_on_TriEdge[ledge][0]]; // global edge node ID 1
        nd2_ID = grid->elem2d[*ielem].nodes[grid->nd_on_TriEdge[ledge][1]]; // global eege node ID 2
        if (DEBUG_SEARCH) printf("---- starts on edge with nodes %d %d \n",nd1_ID,nd2_ID);
        
        // case 2-a :: intersection falls on opposite node
        if (DEBUG_SEARCH) printf("---- going to nodeLoopCheck2D || ielem: %d || p2: %f %f %f || pi: %f %f %f\n",*ielem,p2.x,p2.y,p2.z,pi->x,pi->y,pi->z);
        intersect_flag = nodeLoopCheck2D(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi);
        if (DEBUG_SEARCH) printf("---- finished nodeLoopCheck2D intersect_flag: %d || ielem: %d || p2: %f %f %f || pi: %f %f %f\n",intersect_flag,*ielem,p2.x,p2.y,p2.z,pi->x,pi->y,pi->z);
        if (intersect_flag != 0) return intersect_flag;
        
        // case 2-b :: intersection falls on an edge (2 other edges possible)
        if (DEBUG_SEARCH) printf("---- going to edgeLoopCheck2D || ielem: %d || p2: %f %f %f || pi: %f %f %f\n",*ielem,p2.x,p2.y,p2.z,pi->x,pi->y,pi->z);
        intersect_flag = edgeLoopCheck2D(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi);
        if (DEBUG_SEARCH) printf("---- finished edgeLoopCheck2D || ielem: %d || p2: %f %f %f || pi: %f %f %f\n",*ielem,p2.x,p2.y,p2.z,pi->x,pi->y,pi->z);
        if (intersect_flag != 0) return intersect_flag;
        
        // case 2-c :: intersection falls in adjacent element; take a small step along the trajectory from p1 and find element
        if (cycle == false) return 0; // do not keep looking at adjacent elements
        if (grid->elem2d[*ielem].edge_flag[ledge] == UNSET_INT)  {
            // particle does not cross node or other edge, so it must have exited grid
            return -2;
        }
        *ielem = grid->elem2d[*ielem].edge_flag[ledge];
        if (DEBUG_SEARCH) printf("-- ends/travels through adjacent element %d \n",*ielem);
        return elementSearch2D(grid,p1,p2,ielem,true,pi);
        
    } else {
        //--------------------------------------------------------------------------------------
        //-- case 1 :: particle starts in element
        //--------------------------------------------------------------------------------------
        
        // -- case 3 :: particle starts within element
        if (DEBUG_SEARCH) printf(" starts within element %d \n",*ielem);
        
        // case 3-a :: intersection falls on another node
        if (DEBUG_SEARCH) printf(" -- going to nodeLoopCheck2D || ielem: %d || p2: %f %f %f || pi: %f %f %f\n",*ielem,p2.x,p2.y,p2.z,pi->x,pi->y,pi->z);
        intersect_flag = nodeLoopCheck2D(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi);
        if (DEBUG_SEARCH) printf(" -- finished nodeLoopCheck2D || ielem: %d || p2: %f %f %f || pi: %f %f %f\n",*ielem,p2.x,p2.y,p2.z,pi->x,pi->y,pi->z);
        if (intersect_flag != 0) return intersect_flag;
        
        // case 3-b :: intersection falls on an edge
        if (DEBUG_SEARCH) printf(" -- going to edgeLoopCheck2D || ielem: %d || p2: %f %f %f || pi: %f %f %f\n",*ielem,p2.x,p2.y,p2.z,pi->x,pi->y,pi->z);
        intersect_flag = edgeLoopCheck2D(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi);
        if (DEBUG_SEARCH) printf(" -- finished edgeLoopCheck2D || ielem: %d || p2: %f %f %f || pi: %f %f %f\n",*ielem,p2.x,p2.y,p2.z,pi->x,pi->y,pi->z);
        if (intersect_flag != 0) return intersect_flag;
    }
    
    if (DEBUG_SEARCH) printf(" never crosses element %d \n",*ielem);
    return 0; // particle trajectory never crosses this element
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

int elementSearch3D(SGRID *grid, SVECT p1, SVECT p2, int *ielem, bool cycle, SVECT *pi) {
    
    assert(*ielem > -1 && *ielem < grid->nelems3d);
    
    int DEBUG_SEARCH = 0;
    if (DEBUG_SEARCH) {
        printf("\n---- ELEMENT SEARCH || p1 elem: %d || p1: {%20.15e,%20.15e,%20.15e) || p2: {%20.15e,%20.15e,%20.15e}\n",
               *ielem,p1.x,p1.y,p1.z,p2.x,p2.y,p2.z);
        int i;
        for (i=0; i<4; i++) {
            printf("---- ELEMENT SEARCH || element node %d || gid: %d || node_flag: %d || {x,y,z} || {%20.15e %20.15e %20.15e}\n",
                   i,grid->elem3d[*ielem].nodes[i],grid->node_flag[grid->elem3d[*ielem].nodes[i]],
                   grid->node[grid->elem3d[*ielem].nodes[i]].x,
                   grid->node[grid->elem3d[*ielem].nodes[i]].y,
                   grid->node[grid->elem3d[*ielem].nodes[i]].z);
        }
    }
    
    int elem_org = *ielem;
    double weights[4] = {0., 0., 0., 0.};
    
    // sanity check, p1 better start in ielem
    int check = UNSET_INT, flag = UNSET_INT, lnd = UNSET_INT, ledge = UNSET_INT, lface = UNSET_INT;
    //flag = compute_interpolate_weights_3D_tetrahedron(grid, *ielem, p1.x, p1.y, p1.z, weights, &lnd, &ledge, &lface);
    //assert(flag==1);
    
    //--------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------
    // check if particle trajectory ends (p2) on element
    check = UNSET_INT; flag = UNSET_INT; lnd = UNSET_INT; ledge = UNSET_INT; lface = UNSET_INT;
    flag = compute_interpolate_weights_3D_tetrahedron(grid, *ielem, p2.x, p2.y, p2.z, weights, &lnd, &ledge, &lface);
    if(flag == 1) { // particle final position found in element
        //printf("calling check_for_external_boundary\n");
        return check_for_external_boundary(grid,p1,p2,ielem,lnd,ledge,lface,DEBUG_SEARCH,pi);
    }
    
    //--------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------
    // make sure particle trajectory starts (p1) on original element
    check = UNSET_INT; flag = UNSET_INT; lnd = UNSET_INT; ledge = UNSET_INT; lface = UNSET_INT;
    flag = compute_interpolate_weights_3D_tetrahedron(grid, *ielem, p1.x, p1.y, p1.z, weights, &lnd, &ledge, &lface);
    //printf("\n ---- flag: %d\n",flag);
    if(flag != 1) {
        if (DEBUG_SEARCH) printf(" Could not find particle start on given element!");
        return 0; // can happen when trajectory was on edge and node connectivity is used
    }
    
    //--------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------
    int i, j, k, intersect_flag = UNSET_INT, nd, nd_ID, nd1_ID, nd2_ID, nd3_ID, nd1_1_ID, nd2_1_ID, nd3_1_ID, ndt;
    int node_search = UNSET_INT;
    SVECT nd1, nd2, nd3;
    bool isPointOnLine = false;
    
    //--------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------
    // calculate intersection point
    //--------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------
    if (lnd != UNSET_INT) {
        //--------------------------------------------------------------------------------------
        //-- case 1 :: particle starts on node
        //--------------------------------------------------------------------------------------
        
        nd_ID = grid->elem3d[*ielem].nodes[lnd]; // starting node global ID
        if (DEBUG_SEARCH) printf(" starts on node %d ",nd_ID);
        
        // case 1-a :: intersection falls on another node
        intersect_flag = nodeLoopCheck3D(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi);
        if (intersect_flag != 0) return intersect_flag;
        
        // case 1-b :: intersection falls on one of three opposite edges
        intersect_flag = edgeLoopCheck3D(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi);
        if (intersect_flag != 0) return intersect_flag;
        
        // case 1-c :: intersection falls on an face (can only intersect face opposite node)
        intersect_flag = faceLoopCheck(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi,3);
        if (intersect_flag != 0) return intersect_flag;
        
        // case 1-d :: no intersection in this element
        if (cycle == false) return 0; // intersection *NOT* found, do not keep looking at adjacent elements
        return findElementsAroundNode(grid,p1,p2,ielem,DEBUG_SEARCH,nd_ID,elem_org,pi,elementSearch3D);
        
        
    } else if (ledge != UNSET_INT) {
        
        //--------------------------------------------------------------------------------------
        //-- case 1 :: particle starts on edge
        //--------------------------------------------------------------------------------------
        
        nd1_ID = grid->elem3d[*ielem].nodes[grid->nd_on_TetEdge[ledge][0]]; // global edge node ID 1
        nd2_ID = grid->elem3d[*ielem].nodes[grid->nd_on_TetEdge[ledge][1]]; // global eege node ID 2
        if (DEBUG_SEARCH) printf(" starts on local edge of element %d with nodes %d %d ",
                                 *ielem,nd1_ID,nd2_ID);
        
        // case 2-a :: intersection falls on node
        if (DEBUG_SEARCH) printf("ledge :: intersection of node check\n");
        intersect_flag = nodeLoopCheck3D(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi);
        if (intersect_flag != 0) return intersect_flag;
        
        // case 2-b :: intersection falls on another edge
        if (DEBUG_SEARCH) printf("ledge :: intersection of edge check\n");
        intersect_flag = edgeLoopCheck3D(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi);
        if (intersect_flag != 0) return intersect_flag;

        // case 2-c :: intersection falls on an face
        if (DEBUG_SEARCH) printf("ledge :: intersection of face check\n");
        intersect_flag = faceLoopCheck(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi,2) ;
        if (intersect_flag != 0) return intersect_flag;
        
        // case 2-d :: intersection falls on another edge connected element
        if (cycle == false) return 0;  // intersection *NOT* found, do not keep looking at adjacent elements
        if (DEBUG_SEARCH) printf("ledge :: intersection on another edge connected element check\n");
        return findElementsAroundNode(grid,p1,p2,ielem,DEBUG_SEARCH,nd1_ID,elem_org,pi,elementSearch3D);
  
        
    } else if (lface != UNSET_INT) {
        
        //--------------------------------------------------------------------------------------
        //-- case 3 :: particle starts on face
        //--------------------------------------------------------------------------------------
        
        nd1_ID = grid->elem3d[*ielem].nodes[nd_on_fc[lface][0]];
        nd2_ID = grid->elem3d[*ielem].nodes[nd_on_fc[lface][1]];
        nd3_ID = grid->elem3d[*ielem].nodes[nd_on_fc[lface][2]];
        if (DEBUG_SEARCH) printf(" starts on local face %d of element %d with nodes %d %d %d ",lface,*ielem,nd1_ID,nd2_ID,nd3_ID);
        
        // case 3-a :: intersection falls node
        intersect_flag = nodeLoopCheck3D(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi);
        if (intersect_flag != 0) return intersect_flag;

        // case 3-b :: intersection falls on an edge not on this face
        intersect_flag = edgeLoopCheck3D(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi);
        if (intersect_flag != 0) return intersect_flag;
        
        // case 3-c :: intersection falls on another face
        intersect_flag = faceLoopCheck(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi,1) ;
        if (intersect_flag != 0) return intersect_flag;
        
        // case 3-d :: intersection falls in adjacent element
        if (cycle == false) return 0;  // intersection *NOT* found, do not keep looking at adjacent elements
        
        if (grid->elem3d[*ielem].face_flag[lface] == UNSET_INT)  {
            // particle does not cross node or other edge, so it must have exited grid
            *pi = p1;
            return -2;
        }
        // if not propagate the trajectory further
        *ielem = grid->elem3d[*ielem].face_flag[lface];
        if (DEBUG_SEARCH) printf("and ends/travels through adjacent element %d \n",*ielem);
        return elementSearch3D(grid,p1,p2,ielem,true,pi);
        
    } else {
        
        //--------------------------------------------------------------------------------------
        //-- case 4 :: particle starts within element
        //--------------------------------------------------------------------------------------
        if (DEBUG_SEARCH) printf(" starts within element ");
        
        // case 4-a :: intersection falls on another node
        intersect_flag = nodeLoopCheck3D(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi);
        if (intersect_flag != 0) return intersect_flag;
        
        // case 4-b :: intersection falls on an edge
        intersect_flag = edgeLoopCheck3D(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi);
        if (intersect_flag != 0) return intersect_flag;
        
        // case 4-c :: intersection falls on an face
        intersect_flag = faceLoopCheck(grid,p1,p2,elem_org,DEBUG_SEARCH,ielem,pi,0) ;
        if (intersect_flag != 0) return intersect_flag;

    }
    
    if (DEBUG_SEARCH) printf(" particle trajectory never crosses this element ");
    return 0;
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

// when p2 ends on element face, edge or node, check if it's an external edge for future velocity projections
int check_for_external_boundary(SGRID *grid, SVECT p1, SVECT p2, int *ielem, int lnd, int ledge, int lface, int DEBUG_SEARCH, SVECT *pi) {
    
    int check = -999, nd = -999;
    
    if (lface != UNSET_INT) {
        if (grid->elem3d[*ielem].face_flag[lface] == UNSET_INT) {
            if (DEBUG_SEARCH) printf("------ particle trajectory ends on external local face %d of element %d\n",lface,*ielem);
            // obtain intersection point
            int nd1_ID = grid->elem3d[*ielem].nodes[nd_on_fc[lface][0]];
            int nd2_ID = grid->elem3d[*ielem].nodes[nd_on_fc[lface][1]];
            int nd3_ID = grid->elem3d[*ielem].nodes[nd_on_fc[lface][2]];
            SVECT nd1 = grid->node[nd1_ID];
            SVECT nd2 = grid->node[nd2_ID];
            SVECT nd3 = grid->node[nd3_ID];
            int intersect_flag = intersect3D_SegTriangle(p1,p2,nd1,nd2,nd3,pi);
            if (intersect_flag != 1 && intersect_flag != 2) {
                printf("ERROR :: ray segement is off 3D face!\n");
                printf("p1: %20.15f %20.15f %20.15f\n",p1.x,p1.y,p1.z);
                printf("p2: %20.15f %20.15f %20.15f\n",p2.x,p2.y,p2.z);
                printf("nd1: %20.15f %20.15f %20.15f\n",nd1.x,nd1.y,nd1.z);
                printf("nd2: %20.15f %20.15f %20.15f\n",nd2.x,nd2.y,nd2.z);
                printf("nd3: %20.15f %20.15f %20.15f\n",nd3.x,nd3.y,nd3.z);
                printf("intersect_flag: %d\n",intersect_flag);
                exit(-1);
            }
            assert(intersect_flag == 1 || intersect_flag == 2);
            if (intersect_flag == 2) {
                *pi = p2; // segments are colinear along external boundary edge
            }
//            printf("right before returning -2\n");
            return -2; // external boundary crossed
        } else {
            if (DEBUG_SEARCH)  printf("------ particle trajectory ends inside element %d\n",*ielem);
            return -1;
        }
    }
    
    if (ledge != UNSET_INT) {
        if (grid->ndim == 2) {
            check = grid->elem2d[*ielem].edge_flag[ledge];
            if (check == UNSET_INT) {
                // obtain intersection point
                int nd1_ID = grid->elem2d[*ielem].nodes[grid->nd_on_TriEdge[ledge][0]];
                int nd2_ID = grid->elem2d[*ielem].nodes[grid->nd_on_TriEdge[ledge][1]];
                if (DEBUG_SEARCH)  printf("------ particle trajectory ends on external 2D local edge %d of element %d with global nodes %d %d \n",ledge,*ielem,nd1_ID,nd2_ID);
                SVECT nd1 = grid->node[nd1_ID];
                SVECT nd2 = grid->node[nd2_ID];
                int intersect_flag = intersect2D_SegSeg(p1,p2,nd1,nd2,pi);
                assert(intersect_flag == 1 || intersect_flag == -1); // make sure lines intersect or are colinear
                if (intersect_flag == -1) {
                    *pi = p2; // segments are colinear along external boundary edge
                }
                return -2; // external boundary crossed
            }
        } else  if (grid->ndim == 3) {
            check = grid->elem3d[*ielem].edge_flag[ledge];
            if (check == 1) {
                // obtain intersection point
                int nd1_ID = grid->elem3d[*ielem].nodes[grid->nd_on_TetEdge[ledge][0]];
                int nd2_ID = grid->elem3d[*ielem].nodes[grid->nd_on_TetEdge[ledge][1]];
                if (DEBUG_SEARCH)  printf("------ particle trajectory ends on external 3D local edge %d of element %d with global nodes %d %d \n",ledge,*ielem,nd1_ID,nd2_ID);
                SVECT nd1 = grid->node[nd1_ID];
                SVECT nd2 = grid->node[nd2_ID];
                int intersect_flag = intersect3D_SegSeg(p1,p2,nd1,nd2,pi);
                assert(intersect_flag == 1 || intersect_flag == -1); // make sure lines intersect or are colinear
                if (intersect_flag == -1) {
                    *pi = p2; // segments are colinear along external boundary edge
                }
                return -2; // external boundary crossed
            }
        }
        else {
            if (DEBUG_SEARCH)  printf("------ particle trajectory ends inside element %d\n",*ielem);
            return -1;
        }
    }
    
    if (lnd != UNSET_INT) {
        if (grid->ndim == 3) {
            nd = grid->elem3d[*ielem].nodes[lnd];
        } else {
            nd = grid->elem2d[*ielem].nodes[lnd];
        }
        if (grid->node_flag[nd] == 1) {
            if (DEBUG_SEARCH) printf("------ particle trajectory ends on external local node %d of element %d\n",nd,*ielem);
            *pi = grid->node[nd];
            return -2; // external boundary crossed
        } else {
            if (DEBUG_SEARCH)  printf("------ particle trajectory ends inside element %d\n",*ielem);
            return -1;
        }
    }
    if (DEBUG_SEARCH) printf("------ particle found in element %d\n",*ielem);
    return -1;
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
// dealing with particles moving elements after wgrid displacement
// made for tetrahedrons!
// assumes no bottom node displacement!
// returns 1 if the partticle is displaced off grid, 0 otherwise
int checkElementAfterDisplacement(SGRID *grid, SVECT *p, int *ielem) {
    
    //printf("**** entering checkElementAfterDisplacement\n"); fflush(stdout);
    
    int elem_store = *ielem;
    
    SVECT pi; // node/edge/face intersection point
    int i,check = UNSET_INT, lnd = UNSET_INT, ledge = UNSET_INT, lface = UNSET_INT, flag = 0;
    double weights[4];

    // is the particle on the element given?
    //printf("ielem: %d r: %10.5f %10.5f %10.5f\n",*ielem,p->x,p->y,p->z); fflush(stdout);
    check = UNSET_INT; flag = UNSET_INT; lnd = UNSET_INT; ledge = UNSET_INT; lface = UNSET_INT;
    flag = compute_interpolate_weights_3D_tetrahedron(grid, *ielem, p->x, p->y, p->z, weights, &lnd, &ledge, &lface);
    //if (flag == 1) printf("Still in same element\n");
    if(flag == 1) {
        //printf("**** exiting checkElementAfterDisplacement\n");
        return 0; // particle never left element after grid displacement
    }
    
    SVECT p1, p2; // ray endpoints

    SVECT nd1, nd2, nd3;
    double denom = 0., z, zmin = 1e12;
    for (i=0; i<4; i++) { // loop over faces
        nd1 =  grid->node[grid->elem3d[*ielem].nodes[nd_on_fc[i][0]]];
        nd2 =  grid->node[grid->elem3d[*ielem].nodes[nd_on_fc[i][1]]];
        nd3 =  grid->node[grid->elem3d[*ielem].nodes[nd_on_fc[i][2]]];
        
        denom = (nd1.x*nd2.y - nd1.x*nd3.y - nd2.x*nd1.y + nd2.x*nd3.y + nd3.x*nd1.y - nd3.x*nd2.y);
        if (denom > 1e-10) {
            z = (nd1.x*nd2.y*nd3.z - nd1.x*nd3.y*nd2.z - nd2.x*nd1.y*nd3.z + nd2.x*nd3.y*nd1.z + nd3.x*nd1.y*nd2.z - nd3.x*nd2.y*nd1.z -
                 p->x*(nd1.y*nd2.z - nd1.y*nd3.z - nd2.y*nd1.z + nd2.y*nd3.z + nd3.y*nd1.z - nd3.y*nd2.z) +
                 p->y*(nd1.x*nd2.z - nd1.x*nd3.z - nd2.x*nd1.z + nd2.x*nd3.z + nd3.x*nd1.z - nd3.x*nd2.z))/denom;
            if (z < zmin) zmin  = z;
        }
    }
    
    // particle is still on grid and has changed elements. Find new element
    // by starting on center of original element and tracking to element position on displaced grid
    *ielem = elem_store;
    p1.x = p->x;
    p1.y = p->y;
    p1.z = zmin;
    pi = p1;
    flag = elementSearch3D(grid,p1,*p,ielem,true,&pi);
    if (flag == -2 && pi.z < p->z) p->z = pi.z;
    if (flag == 0) {
        printf("ERROR: flag: %d p1: {%20.15e %20.15e %20.15f} || p: {%20.15e %20.15e %20.15f} \n",flag,p1.x,p1.y,p1.z,p->x,p->y,p->z);
        exit(-1);
    }
    assert(flag != 0); // better have found an element
    //printf("**** exiting checkElementAfterDisplacement\n");
    if (elem_store != *ielem) return 1;
    else return 0;
}
