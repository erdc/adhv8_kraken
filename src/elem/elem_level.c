/* this routine calculates the level of a 3d element */
/* the element level is the level of the lowest node */

#include "global_header.h"

/* Gajanan gkc - shortening these functions along with addition of prism support. */

//********************************************************************//
//********************************************************************//
//********************************************************************//

int elem3d_level(SELEM_3D elem3d) {
    int elem_lev = 0;		/* temp used for sorting */
    int i;

    for (i=0; i<elem3d.nnodes; i++){
        if (elem3d.levels[i] > elem_lev) elem_lev = elem3d.levels[i];
    }
    
    /* returns the element level */
    return (elem_lev);
}

//********************************************************************//
//********************************************************************//
//********************************************************************//

int elem2d_level(SELEM_2D elem2d) {
    int elem_lev = 0;		/* temp used for sorting */
    int i;

    for (i=0; i<elem2d.nnodes; i++){
        if (elem2d.levels[i] > elem_lev) elem_lev = elem2d.levels[i];
    }

    /* returns the element level */
    return (elem_lev);
}

//********************************************************************//
//********************************************************************//
//********************************************************************//

int elem1d_level(SELEM_1D elem1d) {
    int elem_lev = 0;		/* temp used for sorting */
    int i;

    for (i=0; i<elem1d.nnodes; i++){
        if (elem1d.levels[i] > elem_lev) elem_lev = elem1d.levels[i];
    }

    /* returns the element level */
    return (elem_lev);
}
