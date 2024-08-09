#include "adh.h"
static int DEBUG = OFF;
static int DEBUG_WITH_PICKETS = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Print an AdH grid to screen
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] pgrid          (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void sgrid_printScreen(SGRID *g) {
    int i;

    int npes = g->smpi->npes; // alias
    int myid = g->smpi->myid; // alias
    int buffer1[npes],buffer2[npes],buffer3[npes];
    double dbuffer1[npes],dbuffer2[npes],dbuffer3[npes];

    if (myid == 0) {
        printf("\n");
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        printf("----- Grid: \"%s\" statistics || type: %s || highest dimension: %d\n", g->filename, g->type, g->ndim);
        printf("---------- total # global nodes: %d || original: %d\n", g->macro_nnodes,g->orig_macro_nnodes);
        printf("---------- total # global surface nodes: %d || original: %d\n", g->macro_nnodes_sur,g->orig_macro_nnodes_sur);
        printf("---------- total # global bed nodes: %d || original: %d\n", g->macro_nnodes_bed,g->orig_macro_nnodes_bed);
        printf("---------- total # global 1D elements: %d || original: %d\n", g->macro_nelems1d,g->orig_macro_nelems1d);
        printf("---------- total # global 2D elements: %d || original: %d\n", g->macro_nelems2d,g->orig_macro_nelems2d);
        printf("---------- total # global 3D elements: %d || original: %d\n", g->macro_nelems3d,g->orig_macro_nelems3d);

        printf("---------- total # global Triangle elements: %d\n", g->macro_nTris);
        printf("---------- total # global Quadrilateral elements: %d\n", g->macro_nQuads);
        printf("---------- total # global Prism elements: %d\n", g->macro_nPrisms);
        printf("---------- total # global Tetrahedral elements: %d\n", g->macro_nTets);

        //printf("---------- total 1D grid edge length: %20.10e\n", g->mesh_length);
        printf("---------- total 2D grid body area: %10.2f\n",   g->mesh_area);
        printf("---------- total 3D grid surface area: %10.2f\n", g->mesh_area_surface);
        printf("---------- total 3D grid bed area: %10.2f\n", g->mesh_area_bed);
        printf("---------- total 3D grid sidewall area (TRIANGLES ONLY): %10.2f\n", g->mesh_area_sidewalls);
        printf("---------- total 3D grid volume: %10.2f\n", g->mesh_volume);
        printf("---------- total grid x-bounds: [%8.2f,%8.2f]\n", g->x_min,g->x_max);
        printf("---------- total grid y-bounds: [%8.2f,%8.2f]\n", g->y_min,g->y_max);
        printf("---------- total grid z-bounds: [%8.2f,%8.2f]\n", g->z_min,g->z_max);
    }

    if (npes > 1) {
#ifdef _MESSG

        // print local nnode values in PE order
        sarray_init_int(buffer1, UNSET_INT); messg_gather_int(0, &g->my_nnodes,  buffer1, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer2, UNSET_INT); messg_gather_int(0, &g->max_nnodes, buffer2, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer3, UNSET_INT); messg_gather_int(0, &g->nnodes_old, buffer3, 1, g->smpi->ADH_COMM);
        if (myid == 0) {
            for (i=0; i<npes; i++) {
                printf("------------ PE: %d || my_nnodes: %d || max_nnodes: %d || nnodes_old: %d\n",i,buffer1[i],buffer2[i],buffer3[i]);
            }
        }

        // print local nnode values in PE order
        sarray_init_int(buffer1, UNSET_INT); messg_gather_int(0, &g->my_nnodes_sur,  buffer1, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer2, UNSET_INT); messg_gather_int(0, &g->max_nnodes_sur, buffer2, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer3, UNSET_INT); messg_gather_int(0, &g->nnodes_sur_old, buffer3, 1, g->smpi->ADH_COMM);
        if (myid == 0) {
            for (i=0; i<npes; i++) {
                printf("------------ PE: %d || surface || my_nnodes: %d || max_nnodes: %d || nnodes_old: %d\n",i,buffer1[i],buffer2[i],buffer3[i]);
            }
        }

        // print local nnode values in PE order
        sarray_init_int(buffer1, UNSET_INT); messg_gather_int(0, &g->my_nnodes_bed,  buffer1, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer2, UNSET_INT); messg_gather_int(0, &g->max_nnodes_bed, buffer2, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer3, UNSET_INT); messg_gather_int(0, &g->nnodes_bed_old, buffer3, 1, g->smpi->ADH_COMM);
        if (myid == 0) {
            for (i=0; i<npes; i++) {
                printf("------------ PE: %d || bed || my_nnodes: %d || max_nnodes: %d || nnodes_old: %d\n",i,buffer1[i],buffer2[i],buffer3[i]);
            }
        }

        // print local 1D element values in PE order
        sarray_init_int(buffer1, UNSET_INT); messg_gather_int(0, &g->nelems1d,  buffer1, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer2, UNSET_INT); messg_gather_int(0, &g->max_nelems1d, buffer2, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer3, UNSET_INT); messg_gather_int(0, &g->nelems1d_old, buffer3, 1, g->smpi->ADH_COMM);
        if (myid == 0) {
            for (i=0; i<npes; i++) {
                printf("------------ PE: %d || nelems1d: %d || max_nelems1d: %d || nelems1d_old: %d\n",i,buffer1[i],buffer2[i],buffer3[i]);
            }
        }

        // print local 2D element values in PE order
        sarray_init_int(buffer1, UNSET_INT); messg_gather_int(0, &g->nelems2d,  buffer1, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer2, UNSET_INT); messg_gather_int(0, &g->max_nelems2d, buffer2, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer3, UNSET_INT); messg_gather_int(0, &g->nelems2d_old, buffer3, 1, g->smpi->ADH_COMM);
        if (myid == 0) {
            for (i=0; i<npes; i++) {
                printf("------------ PE: %d || nelems2d: %d || max_nelems2d: %d || nelems2d_old: %d\n",i,buffer1[i],buffer2[i],buffer3[i]);
            }
        }

        // print local 3D element values in PE order
        sarray_init_int(buffer1, UNSET_INT); messg_gather_int(0, &g->nelems3d,  buffer1, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer2, UNSET_INT); messg_gather_int(0, &g->max_nelems3d, buffer2, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer3, UNSET_INT); messg_gather_int(0, &g->nelems3d_old, buffer3, 1, g->smpi->ADH_COMM);
        if (myid == 0) {
            for (i=0; i<npes; i++) {
                printf("------------ PE: %d || nelems3d: %d || max_nelems3d: %d || nelems3d_old: %d\n",i,buffer1[i],buffer2[i],buffer3[i]);
            }
        }


        // print local element types for 1D
        sarray_init_int(buffer1, UNSET_INT); messg_gather_int(0, &g->my_nelems1d,  buffer1, 1, g->smpi->ADH_COMM);
        if (myid == 0) {
            for (i=0; i<npes; i++) {
                printf("------------ PE: %d || my_nelems1d: %d \n",i,buffer1[i]);
            }
        }
        // print local element types for 2D
        sarray_init_int(buffer1, UNSET_INT); messg_gather_int(0, &g->my_nelems2d,  buffer1, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer2, UNSET_INT); messg_gather_int(0, &g->my_nQuads, buffer2, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer3, UNSET_INT); messg_gather_int(0, &g->my_nTris, buffer3, 1, g->smpi->ADH_COMM);
        if (myid == 0) {
            for (i=0; i<npes; i++) {
                printf("------------ PE: %d || my_nelems2d: %d || my_nQuads: %d || my_nTris: %d\n",i,buffer1[i],buffer2[i],buffer3[i]);
            }
        }
        // print local element types for 3D
        sarray_init_int(buffer1, UNSET_INT); messg_gather_int(0, &g->my_nelems3d,  buffer1, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer2, UNSET_INT); messg_gather_int(0, &g->my_nTets, buffer2, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer3, UNSET_INT); messg_gather_int(0, &g->my_nPrisms, buffer3, 1, g->smpi->ADH_COMM);
        if (myid == 0) {
            for (i=0; i<npes; i++) {
                printf("------------ PE: %d || my_nelems3d: %d || my_nTets: %d || my_nPrisms: %d\n",i,buffer1[i],buffer2[i],buffer3[i]);
            }
        }


        // print local mesh areas
        sarray_init_dbl(dbuffer1, 0); messg_gather_dbl(0, &g->my_mesh_area_surface,   dbuffer1, 1, g->smpi->ADH_COMM);
        sarray_init_dbl(dbuffer2, 0); messg_gather_dbl(0, &g->my_mesh_area_bed,       dbuffer2, 1, g->smpi->ADH_COMM);
        sarray_init_dbl(dbuffer3, 0); messg_gather_dbl(0, &g->my_mesh_area_sidewalls, dbuffer3, 1, g->smpi->ADH_COMM);
        if (myid == 0) {
            for (i=0; i<npes; i++) {
                printf("------------ PE: %d || surface area: %10.2f || bed area: %10.2f || sidewall area: %10.2f\n",i,dbuffer1[i],dbuffer2[i],dbuffer3[i]);
            }
        }

        // print local mesh volumes
        sarray_init_dbl(dbuffer2, 0); messg_gather_dbl(0, &g->my_mesh_area,   dbuffer2, 1, g->smpi->ADH_COMM);
        sarray_init_dbl(dbuffer3, 0); messg_gather_dbl(0, &g->my_mesh_volume, dbuffer3, 1, g->smpi->ADH_COMM);
        if (myid == 0) {
            for (i=0; i<npes; i++) {
                printf("------------ PE: %d || 2D grid body area: %10.2f || 3D grid body volume: %10.2f\n",i,dbuffer2[i],dbuffer3[i]);
            }
        }
#endif

    } else {

        printf("---------- nnodes: %d || max_nnodes: %d || nnodes_old: %d\n",g->nnodes,g->max_nnodes,g->nnodes_old);
        printf("---------- nelems1d: %d || max_nelems1d: %d || nelems1d_old: %d\n",g->nelems1d,g->max_nelems1d,g->nelems1d_old);
        printf("---------- nelems2d: %d || max_nelems2d: %d || nelems2d_old: %d\n",g->nelems2d,g->max_nelems2d,g->nelems2d_old);
        printf("---------- nelems3d: %d || max_nelems3d: %d || nelems3d_old: %d\n",g->nelems3d,g->max_nelems3d,g->nelems3d_old);


    }
    if (myid == 0) {
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    }

//    if (myid == 1) {
//        int i;
//        for (i=0; i<g->nnodes; i++) {
//            snode_printScreen(g->node[i]);
//        }
//        for (i=0; i<g->nelems1d; i++) {
//            selem1d_printScreen(&g->elem1d[i]);
//        }
//        for (i=0; i<g->nelems2d; i++) {
//            selem2d_printScreen(&g->elem2d[i]);
//        }
//        for (i=0; i<g->nelems3d; i++) {
//            selem3d_printScreen(&g->elem3d[i]);
//        }
//    }


}

