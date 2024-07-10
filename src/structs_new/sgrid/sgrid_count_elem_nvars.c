#include "global_header.h"
static int DEBUG = OFF;
static int DEBUG_WITH_PICKETS = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Counts the number of max DOFs for each element and stores the count in SELEMs
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[in] supmod           (SUPER_MODEL *)  an AdH supermodels
 * @param[inout]       g           (SGRID *)  an ADH grid
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Column 1 -- surface water      ||  Colum 2 -- GW      ||  Column 3 -- TRANSPORT
//             0 = none                          0 = on                  0 - Ntrans
//             1 = SW 1D                         1 = off
//             2 = SW 2D
//             3 = SW 3D
//             4 = NS 3D
//             5 = OF 2D

void sgrid_countElemNvars(SUPER_MODEL *smod, SGRID *g) {
    int score, div, digit;

    for (ie=0; ie<g->nelems3d; ie++) {

        score = smod->physics[g->elem3d[ie].mat];
        for (div = 1; div <= score; div *= 10);

        g->elem3d[ie].nvars = 0;

        // surface water
        div /= 10;
        digit = score/div;
        score %= div;
        printf("surface_water: %d\n",digit);
        if (digit == 1) g->elem3d[ie].nvars  += 2; // 1D SW
        if (digit == 2) g->elem3d[ie].nvars  += 3; // 2D SW
        if (digit == 3) g->elem3d[ie].nvars  += 3; // 3D SW
        if (digit == 4) g->elem3d[ie].nvars  += 3; // 3D NS
        if (digit == 5) g->elem3d[ie].nvars  += 1; // 2D OF

        // ground water
        div /= 10;
        digit = score/div;
        score %= div;
        printf("ground_water: %d\n",digit);
        g->elem3d[ie].nvars += digit; // 1D/2D/3D GROUNDWATER

        // transport
        div /= 10;
        digit = score/div;
        score %= div;
        printf("transport: %d\n",digit);
        g->elem3d[ie].nvars += digit; // 1D/2D/3D TRANSPORT

        printf("g->elem3d[%d].nvars: %d\n",ie,g->elem3d[ie].nvars);
    }

    for (ie=0; ie<g->nelems2d; ie++) {

        score = smod->physics[g->elem2d[ie].mat];
        for (div = 1; div <= score; div *= 10);

        g->elem2d[ie].nvars = 0;

        // surface water
        div /= 10;
        digit = score/div;
        score %= div;
        printf("surface_water: %d\n",digit);
        if (digit == 1) g->elem2d[ie].nvars  += 2; // 1D SW
        if (digit == 2) g->elem2d[ie].nvars  += 3; // 2D SW
        if (digit == 3) g->elem2d[ie].nvars  += 3; // 3D SW
        if (digit == 4) g->elem2d[ie].nvars  += 3; // 3D NS
        if (digit == 5) g->elem2d[ie].nvars  += 1; // 2D OF

        // ground water
        div /= 10;
        digit = score/div;
        score %= div;
        printf("ground_water: %d\n",digit);
        g->elem2d[ie].nvars += digit; // 1D/2D/3D GROUNDWATER

        // transport
        div /= 10;
        digit = score/div;
        score %= div;
        printf("transport: %d\n",digit);
        g->elem2d[ie].nvars += digit; // 1D/2D/3D TRANSPORT

        printf("g->elem2d[%d].nvars: %d\n",ie,g->elem2d[ie].nvars);
    }

    for (ie=0; ie<g->nelems1d; ie++) {

        score = smod->physics[g->elem1d[ie].mat];
        for (div = 1; div <= score; div *= 10);

        g->elem1d[ie].nvars = 0;

        // surface water
        div /= 10;
        digit = score/div;
        score %= div;
        printf("surface_water: %d\n",digit);
        if (digit == 1) g->elem1d[ie].nvars  += 2; // 1D SW
        if (digit == 2) g->elem1d[ie].nvars  += 3; // 2D SW
        if (digit == 3) g->elem1d[ie].nvars  += 3; // 3D SW // CJT - technically not needed
        if (digit == 4) g->elem1d[ie].nvars  += 3; // 3D NS // CJT - technically not needed
        if (digit == 5) g->elem1d[ie].nvars  += 1; // 2D OF

        // ground water
        div /= 10;
        digit = score/div;
        score %= div;
        printf("ground_water: %d\n",digit);
        g->elem1d[ie].nvars += digit; // 1D/2D/3D GROUNDWATER

        // transport
        div /= 10;
        digit = score/div;
        score %= div;
        printf("transport: %d\n",digit);
        g->elem1d[ie].nvars += digit; // 1D/2D/3D TRANSPORT

        printf("g->elem1d[%d].nvars: %d\n",ie,g->elem1d[ie].nvars);
    }

}

