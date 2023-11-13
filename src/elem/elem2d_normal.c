#include "global_header.h"

SVECT elem2d_normal(SNODE *node) {

  SVECT nrml;            /* the normal */
  SVECT side1, side2;        /* two sides - the cross product is the normal */

    /* computes two side vectors */
    VT_3D_INIT(side1); /* initialize to zero */
    VT_3D_INIT(side2); /* initialize to zero */
    side1.x = node[1].x - node[0].x;
    side1.y = node[1].y - node[0].y;
    side1.z = node[1].z - node[0].z;
    side2.x = node[2].x - node[0].x;
    side2.y = node[2].y - node[0].y;
    side2.z = node[2].z - node[0].z;

    /* computes their cross product */
    //nrml[iface].x = side1.y * side2.z - side1.z * side2.y;
    //nrml[iface].y = side1.z * side2.x - side1.x * side2.z;
    //nrml[iface].z = side1.x * side2.y - side1.y * side2.x;

    return nrml;
}
