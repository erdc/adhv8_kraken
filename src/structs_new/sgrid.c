/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  sgrid.c This file collects methods of the SGRID structure */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//#include "global_header.h"
#include "local_header.h"

// file scope variables
static int DEBUG = OFF;
static int DEBUG_WITH_PICKETS = OFF;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Allocates an AdH grid
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] pgrid          (SGRID **)  a double pointer to an AdH grid
 * @param[in]  nnodes        (int) the number of residential + ghost nodes in the grid
 * @param[in]  nelems1d    (int) the number of residential + ghost 1D elements on the grid
 * @param[in]  nelems2d    (int) the number of residential + ghost 2D elements on the grid
 * @param[in]  nelems3d    (int) the number of residential + ghost 3D elements on the grid
 * @param[in]  type             (char *)  type of grid (unstructured, columnar, mixed, etc.)
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void sgrid_alloc_init(SGRID **pgrid, char *filename, char *type, int nnodes, int nelems1d, int nelems2d, int nelems3d, int *nSubMods1d, int *nSubMods2d, int *nSubMods3d) {
    
    int ie;
    
    // input integrity check
    assert(nnodes > 1); // CJT: Could have 1, 1D element
    assert(nelems1d + nelems2d + nelems3d > 0);
    
    // allocate grid
    (*pgrid) = (SGRID *) tl_alloc(sizeof(SGRID), 1);
    SGRID *grid = (*pgrid);  // alias
    grid->filename = filename;
    grid->type = type;
    grid->nnodes = nnodes;
    grid->nelems1d = nelems1d;
    grid->nelems2d = nelems2d;
    grid->nelems3d = nelems3d;
    grid->max_nnodes = nnodes;
    grid->max_nelems1d = nelems1d;
    grid->max_nelems2d = nelems2d;
    grid->max_nelems3d = nelems3d;
    
    // allocated nodes and elements and physics on each element
    if (grid->nelems1d > 0) {
        selem1d_init_alloc_array(&(grid->elem1d), grid->nelems1d); // allocate elements
        grid->nSubMods1d = nSubMods1d; // pointer assignment
        elem_physics_alloc_init(grid->elem1d_physics,nelems1d,grid->nSubMods1d); // allocate element physics
    }
    if (grid->nelems2d > 0) {
        selem2d_init_alloc_array(&(grid->elem2d), grid->nelems2d); // allocate elements
        grid->nSubMods2d = nSubMods2d; // pointer assignment
        elem_physics_alloc_init(grid->elem2d_physics,nelems2d,grid->nSubMods2d); // allocate element physics
    }
    if (grid->nelems3d > 0) {
        selem3d_init_alloc_array(&(grid->elem3d), grid->nelems3d); // allocate elements
        grid->nSubMods3d = nSubMods3d; // pointer assignment
        elem_physics_alloc_init(grid->elem3d_physics,nelems3d,grid->nSubMods3d); // allocate element physics
    }
    snode_init_alloc_array(&(grid->node), grid->nnodes);

    // calculate total grid volume and 3D bounds
    //grid->total_volume =
    //grid->x_min =
    //grid->x_max =
    //grid->y_min =
    //grid->y_max =
    //grid->z_min =
    //grid->z_max =
    
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Frees and AdH grid
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
void sgrid_free(SGRID *grid) {
    int ie, inode;
    
    if (grid->elem1d != NULL && grid->nelems1d > 0) {
        elem_physics_free(grid->elem1d_physics,grid->nelems1d,grid->nSubMods1d); // free element physics
        selem1d_free_array(grid->elem1d, grid->max_nelems1d);
    }
    if (grid->elem2d != NULL && grid->nelems2d > 0) {
        elem_physics_free(grid->elem2d_physics,grid->nelems2d,grid->nSubMods2d); // free element physics
        selem2d_free_array(grid->elem2d, grid->max_nelems2d);
    }
    if (grid->elem3d != NULL && grid->nelems3d > 0) {
        elem_physics_free(grid->elem3d_physics,grid->nelems3d,grid->nSubMods3d); // free element physics
        selem3d_free_array(grid->elem3d, grid->max_nelems3d);
    }
    if (grid->nnodes > 0) {
        for (inode=0; inode<grid->max_nnodes; inode++) {
            snode_free(&(grid->node[inode]));
        }
        grid->node = (SNODE *) tl_free(sizeof(SNODE), grid->max_nnodes, grid->node);
    }
    
    grid = (SGRID *) tl_free(sizeof(SGRID), 1, grid);
    
}

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
void sgrid_printScreen(SGRID *grid) {
    
    printf("\n");
    printf("***********************************************************\n");
    printf("***********************************************************\n");
    printf("----- Grid: %s statistics || type: %s \n", grid->filename, grid->type);
    printf("---------- total number of nodes: %d || refinement max: %d\n", grid->nnodes, grid->max_nnodes);
    printf("---------- total number of 1d elements: %d || refinement max: %d\n", grid->nelems1d, grid->max_nelems1d);
    printf("---------- total number of 2d elements: %d || refinement max: %d\n", grid->nelems2d, grid->max_nelems2d);
    printf("---------- total number of 2d elements: %d || refinement max: %d\n", grid->nelems3d, grid->max_nelems3d);
    printf("---------- total grid volumne: %20.10e\n", grid->total_volume);
    printf("---------- grid x-bounds: %20.10f, %20.10f\n", grid->x_min,grid->x_max);
    printf("---------- grid y-bounds: %20.10f, %20.10f\n", grid->y_min,grid->y_max);
    printf("---------- grid z-bounds: %20.10f, %20.10f\n", grid->z_min,grid->z_max);
    printf("***********************************************************\n");
    printf("***********************************************************\n");
    
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads and stores an AdH grid file
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
void sgrid_read(SGRID **pgrid, FILE *fp, char *filename) {

    // Go through grid first to get counts for allocation
    
    // Allocate grid
    //sgrid_alloc_init(pgrid,filename,type,nnodes,nelems1d,nelems2d,nelems3d,nSubMods1d,nSubMods2d,nSubMods3d);
    
    //
    
    // Now read through again to store
    rewind(fp);
    
    
}

