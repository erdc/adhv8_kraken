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
void sgrid_alloc_init(SGRID **pgrid, char *filename, char *type, int nnodes, int nelems1d, int nelems2d, int nelems3d) {
    
    int ie;

    // input integrity check
    assert(nnodes > 1); // CJT: Could have 1, 1D element
    assert(nelems1d + nelems2d + nelems3d > 0);
    
    // allocate grid
    (*pgrid) = (SGRID *) tl_alloc(sizeof(SGRID), 1);
    SGRID *grid = (*pgrid);  // alias
    strcpy(grid->filename,filename);
    strcpy(grid->type,type);
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
        //grid->nSubMods1d = nSubMods1d; // pointer assignment
        //elem_physics_alloc_init(grid->elem1d_physics,nelems1d,grid->nSubMods1d); // allocate element physics
    }
    if (grid->nelems2d > 0) {
        selem2d_init_alloc_array(&(grid->elem2d), grid->nelems2d); // allocate elements
        //grid->nSubMods2d = nSubMods2d; // pointer assignment
        //elem_physics_alloc_init(grid->elem2d_physics,nelems2d,grid->nSubMods2d); // allocate element physics
    }
    if (grid->nelems3d > 0) {
        selem3d_init_alloc_array(&(grid->elem3d), grid->nelems3d); // allocate elements
        //grid->nSubMods3d = nSubMods3d; // pointer assignment
        //elem_physics_alloc_init(grid->elem3d_physics,nelems3d,grid->nSubMods3d); // allocate element physics
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
        //elem_physics_free(grid->elem1d_physics,grid->nelems1d,grid->nSubMods1d); // free element physics
        selem1d_free_array(grid->elem1d, grid->max_nelems1d);
    }
    if (grid->elem2d != NULL && grid->nelems2d > 0) {
        //elem_physics_free(grid->elem2d_physics,grid->nelems2d,grid->nSubMods2d); // free element physics
        selem2d_free_array(grid->elem2d, grid->max_nelems2d);
    }
    if (grid->elem3d != NULL && grid->nelems3d > 0) {
        //elem_physics_free(grid->elem3d_physics,grid->nelems3d,grid->nSubMods3d); // free element physics
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
    printf("----- Grid: \"%s\" statistics || type: %s \n", grid->filename, grid->type);
    printf("---------- total number of nodes: %d || refinement max: %d\n", grid->nnodes, grid->max_nnodes);
    printf("---------- total number of 1d elements: %d || refinement max: %d\n", grid->nelems1d, grid->max_nelems1d);
    printf("---------- total number of 2d elements: %d || refinement max: %d\n", grid->nelems2d, grid->max_nelems2d);
    printf("---------- total number of 2d elements: %d || refinement max: %d\n", grid->nelems3d, grid->max_nelems3d);
    printf("---------- total grid volumne: %20.10e\n", grid->total_volume);
    printf("---------- grid x-bounds: %20.10f, %20.10f\n", grid->x_min,grid->x_max);
    printf("---------- grid y-bounds: %20.10f, %20.10f\n", grid->y_min,grid->y_max);
    printf("---------- grid z-bounds: %20.10f, %20.10f\n", grid->z_min,grid->z_max);
    
//    int i;
//    for (i=0; i<grid->nnodes; i++){
//        printf("node[%d]: {%f, %f, %f} \n",grid->node[i].id,grid->node[i].x,grid->node[i].y, grid->node[i].z);
//        //snode_printScreen(&(grid->node[i]));
//    }
//    for (i=0; i<grid->nelems1d; i++) {
//
//        //selem1d_printScreen(&(grid->elem1d[i]));
//    }
//    for (i=0; i<grid->nelems2d; i++) {
//        //selem2d_printScreen(&(grid->elem2d[i]));
//    }
//    for (i=0; i<grid->nelems3d; i++) {
//        selem3d_printScreen(&(grid->elem3d[i]));
//    }
    
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
    
    int id, idum;
    char *line = NULL, cdum[20], line_save[100];
    size_t len = 0;
    ssize_t read;
    char *token, *subtoken;
    char delim[] = " ";
    
    // open the geometry file
    //open_input_file(geo,"geometry file",0);
    char *fn = "test.geo";
    FILE *ffp = fopen(fn, "r"); // FOR NOW
    

    // Go through grid first to get counts for allocation
    char grid_type[15];
    int nelems3d=0, nelems2d=0, nelems1d=0, nnodes=0;
    while ((read = getline(&line, &len, ffp)) != -1) {
        strcpy(line_save,line);
        //printf("line: %s \n",line_save);
        token = strtok(line,delim); token[strcspn(token, "\n")] = 0;
        //printf("card: %s \n",token);
        if (strcmp(token, "GRID") == 0) {
            subtoken = strtok(NULL, delim);  subtoken[strcspn(subtoken, "\n")] = 0;
            if (strcmp(subtoken,"COLUMNAR")!=0 && strcmp(subtoken,"UNSTRUCTURED")!=0) {
                tl_error(">> input grid must be of type either UNSTRUCTURED or COLUMNAR\n");
            }
            strcpy(grid_type,subtoken);
        }
        if (strcmp(token, "TET") == 0 || strcmp(token, "PRISM") == 0 ) {nelems3d++;}
        if (strcmp(token, "TRI") == 0 || strcmp(token, "QUAD") == 0 ) {nelems2d++;}
        if (strcmp(token, "SEG") == 0) {nelems1d++;}
        if (strcmp(token, "ND") == 0) {nnodes++;}
    }
    
    // Allocate grid
    sgrid_alloc_init(pgrid,fn,grid_type,nnodes,nelems1d,nelems2d,nelems3d);
    SGRID *grid = *pgrid;
    
    // Now read through again to store
    rewind(ffp);
    int count_node=0, count_elem1d=0, count_elem2d=0, count_elem3d=0;
    while ((read = getline(&line, &len, ffp)) != -1) {
        strcpy(line_save,line);
        token = strtok(line,delim); token[strcspn(token, "\n")] = 0;
        if (strcmp(token, "TET") == 0) {
            selem3d_alloc(&(grid->elem3d[count_elem3d]),4);
            grid->elem3d[count_elem3d].id = count_elem3d;
            grid->elem3d[count_elem3d].nnodes = 4;
            sscanf(line_save, "%s %d %d %d %d %d",cdum,&idum,
                   &(grid->elem3d[count_elem3d].nodes[0]),
                   &(grid->elem3d[count_elem3d].nodes[1]),
                   &(grid->elem3d[count_elem3d].nodes[2]),
                   &(grid->elem3d[count_elem3d].nodes[3])
                   );
            count_elem3d++;
        } else if (strcmp(token, "PRISM") == 0) {
            selem3d_alloc(&(grid->elem3d[count_elem3d]),6);
            grid->elem3d[count_elem3d].id = count_elem3d;
            grid->elem3d[count_elem3d].nnodes = 6;
            sscanf(line_save, "%s %d %d %d %d %d %d %d",cdum,&idum,
                   &(grid->elem3d[count_elem3d].nodes[0]),
                   &(grid->elem3d[count_elem3d].nodes[1]),
                   &(grid->elem3d[count_elem3d].nodes[2]),
                   &(grid->elem3d[count_elem3d].nodes[3]),
                   &(grid->elem3d[count_elem3d].nodes[4]),
                   &(grid->elem3d[count_elem3d].nodes[5])
                   );
            count_elem3d++;
        } else if (strcmp(token, "TRI") == 0) {
            selem2d_alloc(&(grid->elem2d[count_elem2d]),3);
            grid->elem2d[count_elem2d].id = count_elem2d;
            grid->elem2d[count_elem2d].nnodes = 3;
            sscanf(line_save, "%s %d %d %d %d",cdum,&idum,
                   &(grid->elem2d[count_elem2d].nodes[0]),
                   &(grid->elem2d[count_elem2d].nodes[1]),
                   &(grid->elem2d[count_elem2d].nodes[2])
                   );
            count_elem2d++;
        } else if (strcmp(token, "QUAD") == 0 ) {
            selem2d_alloc(&(grid->elem2d[count_elem2d]),4);
            grid->elem2d[count_elem2d].id = count_elem2d;
            grid->elem2d[count_elem2d].nnodes = 4;
            sscanf(line_save, "%s %d %d %d %d %d",cdum,&idum,
                   &(grid->elem2d[count_elem2d].nodes[0]),
                   &(grid->elem2d[count_elem2d].nodes[1]),
                   &(grid->elem2d[count_elem2d].nodes[2]),
                   &(grid->elem2d[count_elem2d].nodes[3])
                   );
            count_elem2d++;
        } else if (strcmp(token, "SEG") == 0) {
            grid->elem1d[count_elem1d].id = count_elem1d;
            grid->elem1d[count_elem1d].nnodes = 2;
            sscanf(line_save, "%s %d %d %d",cdum,&idum,
                   &(grid->elem1d[count_elem1d].nodes[0]),
                   &(grid->elem1d[count_elem1d].nodes[1])
                   );
            count_elem1d++;
        } else if (strcmp(token, "ND") == 0) {
            grid->node[count_node].id = count_node;
            sscanf(line_save, "%s %d %lf %lf %lf",cdum,&idum,
                   &(grid->node[count_node].x),
                   &(grid->node[count_node].y),
                   &(grid->node[count_node].z));
            count_node++;
            
        }
    }
    assert(count_node == grid->nnodes);
    assert(count_elem1d == grid->nelems1d);
    assert(count_elem2d == grid->nelems2d);
    assert(count_elem3d == grid->nelems3d);
    
    // Close grid file
    fclose(ffp);
}

