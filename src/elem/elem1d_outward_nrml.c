/* in the overland flow calculations we sometimes need to apply boundary conditions via the weak statement */
/* to do this we need the outward normal vector */
/* this routine finds the vector outward normal projected into the z plane */

#include "global_header.h"

void elem1d_outward_nrml(SGRID *grid)
{

    double distance;            /* used as a distance meas and also the dot product */
    double x, x1, x2;           /* x measures */
    double y, y1, y2;           /* y measures */
    double delx;                /* x difference */
    double dely;                /* y difference */
    int ie;                     /* loop counter over the elements */

    /* the outward normal in the z plane is needed to do the overland flow calculations */
    /* right now this is only run in the groundwater flow type */

    /* loop through the 1d elements */
    for (ie = 0; ie < grid->nelems1d; ie++) {

        if (grid->elem1d[ie].string != UNSET_INT) {

            /* here we calculate the z plane outward normal to the 1d element */
            x1 = grid->node[grid->elem1d[ie].nodes[0]].x;
            x2 = grid->node[grid->elem1d[ie].nodes[1]].x;
            y1 = grid->node[grid->elem1d[ie].nodes[0]].y;
            y2 = grid->node[grid->elem1d[ie].nodes[1]].y;
            delx = x2 - x1;
            dely = y2 - y1;
            distance = sqrt(delx * delx + dely * dely);
            grid->elem1d[ie].nrml.x = dely / distance;
            grid->elem1d[ie].nrml.y = -delx / distance;

            /* we only know that we are normal to the 1d element */
            /* we must now make sure that the normal points outward */
            /* to do this we find the 2d centroid of the associated 2d element */
            /* then we contruct a vector from the first node in the element to this centroid */
            /* the dot product of this with our normal vector should be negative */
            /* if not reverse the sign on the normal */

            x = one_3 * (
                    grid->node[grid->elem2d[grid->elem1d[ie].elem2d].nodes[0]].x + 
                    grid->node[grid->elem2d[grid->elem1d[ie].elem2d].nodes[1]].x + 
                    grid->node[grid->elem2d[grid->elem1d[ie].elem2d].nodes[2]].x);

            y = one_3 * (
                    grid->node[grid->elem2d[grid->elem1d[ie].elem2d].nodes[0]].y + 
                    grid->node[grid->elem2d[grid->elem1d[ie].elem2d].nodes[1]].y + 
                    grid->node[grid->elem2d[grid->elem1d[ie].elem2d].nodes[2]].y);

            delx = x - x1;
            dely = y - y1;

            distance = delx * grid->elem1d[ie].nrml.x + dely * grid->elem1d[ie].nrml.y;

            if (distance > 0.) {
                grid->elem1d[ie].nrml.x *= -1.;
                grid->elem1d[ie].nrml.y *= -1.;
            }

        }

    }

}
