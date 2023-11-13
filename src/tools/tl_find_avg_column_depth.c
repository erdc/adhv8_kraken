/* JLH -- Attempting to fix this */
/* Looks like at first it was assumed that the surface nodes would be passed in -- but
 * in fe_sw3_hvel_load, this is being called for *every* 2D element
 */

/* Given a list of 4 integers, one of which is identical to one of the others, return a list of  */
/* the 3 unique integers */

#include "global_header.h"

double tl_find_avg_column_depth(SGRID *grid, int *sur_node, double *displacement) {

    ID_LIST_ITEM *ptr;
    double surface, bottom;     /* the average surface and bottom elevations for a column */
    double ele_top, ele_bot;    /* the top and bottom elevations of a node */
    double depth;               /* the average  */
    int i;                      /* counter */
    int iseg;                   /* segment number */
    int nd;                     /* node number */

    surface = 0.;
    bottom = 0.;
    for (i = 0; i < NDONTRI; i++) {
        iseg = find_vertical_segment(grid, sur_node[i], grid->vertical_hash);
        ptr = grid->vertical_list[iseg];
        nd = ptr->id;
        ele_top = grid->node[nd].z + displacement[nd];
        surface += one_3 * ele_top;

        while (ptr->next != NULL) {
            nd = ptr->id;
            ele_bot = grid->node[nd].z + displacement[nd];
            ptr = ptr->next;
        }
        bottom += one_3 * ele_bot;
    }
    depth = surface - bottom;
    return (depth);

}
