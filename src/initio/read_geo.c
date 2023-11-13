/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This file reads in the finite element grid geometry
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 *  \returns (int) the dimensionality of the grid
 *
 * @param[in] SIO *io :: pointer to an input file data container
 * @param[in] SGRID *grid :: a grid container
 * \note CJT \::  The following grid variables are filled in:
 * \note int macro_nnodes: total number of nodes in the global mesh
 * \note int macro_nelems1d: total number of 1D elements in the global mesh
 * \note int macro_nelems2d: total number of 2D elements in the global mesh
 * \note int macro_nelems3d: total number of 3D elements in the global mesh
 * \note int orig_macro_nelems1d: total number of original 1D elements in the global mesh
 * \note int orig_macro_nelems2d: total number of original 2D elements in the global mesh
 * \note int orig_macro_nelems3d: total number of original 3D elements in the global mesh
 * \note int macro_nelems2d_bed: total number of 2d bed/surface elements across all PEs on the global mesh
 * \note int initial_nnodes_bed;  total number of initial (before adaption) number of ghost + residential bed nnodes
 * \note int nnodes: number of ghost + residential nodes
 * \note int nelems3d: number of ghost + resdiential 3d elements
 * \note int nelems2d: number of ghost + residential 2d elements
 * \note int nelems1d: number of ghost + residential 1d elements
 * \note int nelems1d2d;
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"
#define DIFF_TOL .0000001           // for finding columns
int DEBUG;
int DEBUG_FULL;
int num_ghosts;                     // count of the ghost nodes
NODE_LIST_ITEM *ghost_nodes=NULL;   // linked list of ghost nodes

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_elements_nodes(SIO *io, SGRID *grid, int JUST_COUNT);
void read_element(SIO *io, SGRID *grid, int JUST_COUNT, int elem_dim, int elem_nnodes, char *data);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int read_geo(SIO *io, SGRID *grid) {
    
    grid->nelems1d = 0; grid->nelems2d = 0; grid->nelems3d = 0; grid->nelems1d2d = 0; grid->nnodes = 0;
    grid->macro_nelems1d = 0; grid->macro_nelems2d = 0; grid->macro_nelems3d = 0; grid->macro_nnodes = 0;
    grid->macro_nnodes_sur = 0; grid->my_nnodes_sur = 0; grid->nnodes_sur = 0;
    grid->macro_nnodes_bed = 0; grid->my_nnodes_bed = 0; grid->nnodes_bed = 0;
    
#ifdef _DEBUG
    if (grid->smpi->myid == 0) printf("\n-- Reading geo file\n");
#endif
    
    // FILE DEBUG **************************************************************************
    DEBUG = OFF;
    DEBUG_FULL = OFF;
    
    // allocate ghost node listing *********************************************************
    num_ghosts = 0;
    ghost_nodes = tl_alloc(sizeof(NODE_LIST_ITEM), 1);
    ghost_nodes->local=UNSET_INT;
    ghost_nodes->global=UNSET_INT;
    ghost_nodes->next = NULL;
    
    // routine variables *******************************************************************
    int i, j, k, ipe, id, my_pe_nnodes_surf = 0;
    int myid = grid->smpi->myid;
    int npes = grid->smpi->npes;
    int nodes_per_pe[npes];
    double surf_x, surf_y, x, y, z;
    char line[MAXLINE];                     /* the input line */
    char *data = NULL, *subdata = NULL;     /* the data */
    
    FILE *fp = NULL;
    char *filename = NULL;
    if (io->geo2d.fp != NULL) {
        grid->ndim = 2;
        fp = io->geo2d.fp;
        filename = io->geo2d.filename;
    } else if (io->geo3d.fp != NULL) {
        grid->ndim = 3;
        fp = io->geo3d.fp;
        filename = io->geo3d.filename;
    } else {
        tl_error("FATAL ERROR :: 2d or 3d geo file not found!");
    }
    
    while (fgets(line, MAXLINE, fp) != NULL) {
        io_save_line(io, fp, filename, line);
        if (strip_comments(line) <= 1) {continue;}    /* ignore empty ('\n') or comment lines */
        switch (parse_card(line, &data)) {
            case CARD_GRID:
                switch (parse_card(data, &subdata)) {
                    case CARD_NOCOL:
                        grid->type = UNSTRUCTURED;
                        break;
                    case CARD_COL:
                        grid->type = COLUMNAR;
                        break;
                    default:
                        grid->type = UNSTRUCTURED;
                        break;
                }
                break;
            default:
                break;
        }
    }
    io_save_line(io, NULL, "", "");
    rewind(fp);
    
#ifdef _DEBUG
    assert(grid->smpi->npes >= 1);
    assert(grid->ndim == 2 || grid->ndim == 3);
    assert(grid->type == COLUMNAR || grid->type == UNSTRUCTURED);
#endif
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // now perform a simple partition where the global nodes are divided by the # of PEs
    // to do this, we first calculate start and end node IDs for each PE
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    grid->start_node_ID = (int *) tl_alloc(sizeof(int), grid->smpi->npes);
    grid->end_node_ID = (int *) tl_alloc(sizeof(int), grid->smpi->npes);
    
    
    if (grid->type == COLUMNAR || grid->ndim == 2) {
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // for columnar grids, do decomposition on the surface grid only
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        // read geo file to get the total number of surface nodes
        grid->macro_nnodes = 0;
        grid->macro_nnodes_sur = 0;
        while (fgets(line, MAXLINE, fp) != NULL) {
            io_save_line(io, fp, filename, line);
            if (strip_comments(line) <= 1) {continue;}    /* ignore empty ('\n') or comment lines */
            switch (parse_card(line, &data)) {
                case CARD_ND:
                    grid->macro_nnodes++;
                    id = read_int_field(*io, &data);
                    x = read_dbl_field(*io, &data);
                    y = read_dbl_field(*io, &data);
                    z = read_dbl_field(*io, &data);
                    if (((fabs(fabs(x) - fabs(surf_x)) > DIFF_TOL) || (fabs(fabs(y) - fabs(surf_y)) > DIFF_TOL)) || grid->macro_nnodes_sur == 0) {
                        surf_x=x;
                        surf_y=y;
                        grid->macro_nnodes_sur++;
                    }
                    break;
                default:
                    break;
            }
        }
        io_save_line(io, NULL, "", "");
        rewind(fp);
        
        // the global bed nodes are the same as surface nodes
        grid->macro_nnodes_bed = grid->macro_nnodes_sur;
        
        // now do a simple decomposition based on these surface nodes
        int surface_nodes_per_pe[npes];
        for (i=0; i<npes; i++) {
            surface_nodes_per_pe[i] = grid->macro_nnodes_sur / npes; // not, integer division floored
        }
        surface_nodes_per_pe[npes-1] += grid->macro_nnodes_sur % npes; // add remaining nodes to last PE
        
        // calculate local processor resident surface and bed nodes (same)
        grid->my_nnodes_sur = surface_nodes_per_pe[grid->smpi->myid];
        grid->my_nnodes_bed = surface_nodes_per_pe[grid->smpi->myid];
        
#ifdef _DEBUG
        assert(grid->macro_nnodes_sur % npes == grid->macro_nnodes_sur - (grid->macro_nnodes_sur / npes) * npes);
#endif
        
        if (npes > 1) {
            int ipe=0, k=0, my_pe_nnodes_surf=0;
            while (fgets(line, MAXLINE, fp) != NULL) {
                io_save_line(io, fp, filename, line);
                if (strip_comments(line) <= 1) {    /* ignore empty ('\n') or comment lines */
                    continue;
                }
                switch (parse_card(line, &data)) {
                    case CARD_ND:
                        k++;
                        id = read_int_field(*io, &data);
                        x = read_dbl_field(*io, &data);
                        y = read_dbl_field(*io, &data);
                        z = read_dbl_field(*io, &data);
                        if (ipe == 0 && my_pe_nnodes_surf == 0) {
                            // first surface node
                            my_pe_nnodes_surf = 1;
                            surf_x=x;
                            surf_y=y;
                            grid->start_node_ID[ipe] = 0;
                        } else if (((fabs(fabs(x) - fabs(surf_x)) > DIFF_TOL) || (fabs(fabs(y) - fabs(surf_y)) > DIFF_TOL))) {
                            // new surface node
                            my_pe_nnodes_surf++;
                            surf_x=x;
                            surf_y=y;
                            // start a new PE if need be
                            if (my_pe_nnodes_surf > surface_nodes_per_pe[ipe]) {
                                ipe++;
                                my_pe_nnodes_surf = 1;
                                grid->start_node_ID[ipe] = k-1;
                                grid->end_node_ID[ipe-1] = k-1-1;
                            }
                        }
                        break;
                    default:
                        break;
                }
            }
        } else {
            grid->start_node_ID[0] = 0;
            grid->end_node_ID[0] = grid->macro_nnodes-1;
        }
        
        io_save_line(io, NULL, "", "");
        rewind(fp);
        
    } else {
        //tl_error("are you running NS?  Doesn't work yet -- Corey");
        
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // Fully unstructured simple partition
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // read geo file to get the total number of global nodes
        grid->macro_nnodes = 0;
        while (fgets(line, MAXLINE, fp) != NULL) {
            io_save_line(io, fp, filename, line);
            if (strip_comments(line) <= 1) {continue;}    /* ignore empty ('\n') or comment lines */
            switch (parse_card(line, &data)) {
                case CARD_ND:
                    grid->macro_nnodes++;
                    break;
                default:
                    break;
            }
        }
        io_save_line(io, NULL, "", "");
        rewind(fp);
        
        for (i=0; i<npes; i++) {
            nodes_per_pe[i] = grid->macro_nnodes / npes; // not, integer division floored
        }
        nodes_per_pe[npes-1] += grid->macro_nnodes % npes; // add remaining nodes to last PE
#ifdef _DEBUG
        assert(grid->macro_nnodes % npes == grid->macro_nnodes - (grid->macro_nnodes / npes) * npes);
#endif
        if (npes > 1) {
            grid->start_node_ID[0] = 0;
            grid->end_node_ID[0] = (nodes_per_pe[0]-1);
            for (i=1; i<npes; i++) {
                grid->start_node_ID[i] = grid->end_node_ID[i-1] + 1;
                grid->end_node_ID[i] = grid->start_node_ID[i] + (nodes_per_pe[i]-1);
            }
        } else {
            grid->start_node_ID[0] = 0;
            grid->end_node_ID[0] = grid->macro_nnodes-1;
        }
    }
    grid->end_node_ID[npes-1] = grid->macro_nnodes-1;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // count both the number of global/macro and local nodes and elements, also builds ghost node list!
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    // compute PE residential node count
    grid->my_nnodes =  grid->end_node_ID[myid] - grid->start_node_ID[myid] + 1;
    
    // compute the global/macro and local elements and nodes
    read_elements_nodes(io, grid, YES);
    if (grid->nelems3d > 0) {
        read_faces(io, grid, YES); // YES here set nelems2d = macro_nelems2d temporarily until local 3d elements are read
    }
    
    // compute total number of local nodes = residential (my_nnodes) and ghost nodes
    grid->nnodes = grid->my_nnodes + num_ghosts;
    
    // for adaption allocation purposes
    grid->max_nnodes = grid->nnodes;
    
    // determine grid dimensionality
    grid->ndim = 2; if (grid->macro_nelems3d > 0) grid->ndim = 3;

    /********************************************************************************/
    /* Assert the .3dm and .faces are ok */
    /********************************************************************************/
    assert(grid->nnodes > 0); /* must have > 0 nodes */
    if (grid->ndim == 3) {
        assert(grid->nelems3d > 0); /* some number of 3d elements must be given */
        assert(grid->nelems2d > 0); /* some number of 2d elements must be given */
    }
    if (grid->ndim == 2) assert(grid->nelems2d > 0); /* some number of 3d elements must be given */
    
    /********************************************************************************/
    /* allocate local grid struct arrays  */
    /********************************************************************************/
    snode_init_alloc_array(&(grid->node), grid->nnodes);
    
    /* allocate & initialize grid 1d elements */
    if (grid->nelems1d > 0) selem1d_init_alloc_array(&(grid->elem1d), grid->nelems1d, NDONSEG);
    
    /* allocate & intialize grid 2d elements */
    if (grid->nelems2d > 0) selem2d_init_alloc_array(&(grid->elem2d), grid->nelems2d);
    
    /* allocate & intialize grid 3d elements */
    if (grid->nelems3d > 0) selem3d_init_alloc_array(&(grid->elem3d), grid->nelems3d);
    
    /********************************************************************************/
    /* read grid files and fill in arrays */
    /********************************************************************************/
    // now that grid elements have been allocated, read and store
    read_elements_nodes(io, grid, NO);
    if (grid->nelems3d > 0) {
        read_faces(io, grid, NO); // cjt :: now that 3d elements are in, we can count local faces
        grid->elem2d = (SELEM_2D *) tl_realloc(sizeof(SELEM_2D), grid->nelems2d, grid->macro_nelems2d, grid->elem2d);
    }

#ifdef _DEBUG
    if (DEBUG || DEBUG_FULL) {
        printf("pe: %d || start_node_ID: %d end_node_ID: %d || grid->macro_nnodes: %d grid->my_nnodes: %d ghost_nnodes: %d :: nelems2d: %d nelems3d: %d macro_nelems2d: %d macro_nelems3d: %d \n",
               myid,grid->start_node_ID[myid],grid->end_node_ID[myid],grid->macro_nnodes,grid->my_nnodes,num_ghosts,grid->nelems2d,grid->nelems3d,grid->macro_nelems2d,grid->macro_nelems3d);
#ifdef _MESSG
        fflush(stdout);
        messg_barrier(MPI_COMM_WORLD);
#endif
    }
#endif
    
    /********************************************************************************/
    /* Integrity check grid */
    /********************************************************************************/
    
    // check for gaps in nodes
    for (i = 0; i < grid->nnodes; i++) {
        if (grid->node[i].string == UNSET_INT) {
            printf("nnodes: %d \n",grid->nnodes);
            printf("pe: %d \t local_node: %d\n",myid,i);
            snode_printScreen(grid->node[i]);
            tl_error("Missing node in input geometry file.");
        }
    }
    
    // check for gaps in elements
    for (i = 0; i < grid->nelems2d; i++) {
        if (grid->elem2d[i].mat == UNSET_INT ) {
            printf("pe: %d \t elemid %d \n", i,grid->elem2d[i].id);
            tl_error("Missing 2d element in input geometry file.");
        }
    }
    for (i = 0; i < grid->nelems3d; i++) {
        if (grid->elem3d[i].mat == UNSET_INT) {
            printf("pe: %d \t elemid %d \n", i,grid->elem3d[i].id);
            tl_error("Missing 3d element in input geometry file.");
        }
        if(grid->elem3d[i].nnodes <= 0 || grid->elem3d[i].nnodes > MAX_NNODES_ON_ELEM3D) {
            selem3d_printScreen(&(grid->elem3d[i]));
            printf("nnodes: %d\n",grid->elem3d[i].nnodes);
            tl_error("ERROR: nnodes");
        }
        if(grid->elem3d[i].nnodes_quad <= 0 || grid->elem3d[i].nnodes_quad > MAX_NNODES_ON_ELEM3D_QUAD) {
            selem3d_printScreen(&(grid->elem3d[i]));
            printf("nnodes_quad: %d\n",grid->elem3d[i].nnodes_quad);
            tl_error("ERROR: nnodes_quad");
        }
    }
    
    //    // close the geo and face file
    //    if (io->geo2d.fp != NULL) {
    //        fclose(io->geo2d.fp);
    //    } else if (io->geo3d.fp != NULL) {
    //        fclose(io->geo3d.fp);
    //    }
    //    if (grid->nelems3d > 0) fclose(io->face.fp);
    
#ifdef _DEBUG
    if (grid->smpi->myid == 0) printf("\n-- Finished reading geo file\n");
#endif
    
    // free the ghost node list
    free_node_list(ghost_nodes);
    
    // return grid dimensionality
    return grid->ndim;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_elements_nodes(SIO *io, SGRID *grid, int JUST_COUNT) {
    
    int i, id, local_index, resident_id, owner, surf_id, global_surf_id = 0;
    double xx, yy, zz, surf_x, surf_y;      // nodal coordinates
    char line[MAXLINE];                     // the input line
    char *data = NULL, *subdata = NULL;     // the data
    
    // there grid variable are being counted in this file
    grid->nelems2d = 0; grid->nelems3d = 0;
    grid->macro_nelems2d = 0; grid->macro_nelems3d = 0; grid->macro_nnodes = 0;
    
#ifdef _DEBUG
    if (grid->smpi->myid == 0) {
        printf("\n---- Reading elements and nodes from geo file");
        if (JUST_COUNT==YES) printf(" (for counting purposes only)");
        printf("\n");
    }
#endif
    
    // open geo and face file
    FILE *fp = NULL;
    char *filename = NULL;
    if (io->geo2d.fp != NULL) {
        fp = io->geo2d.fp;
        filename = io->geo2d.filename;
    } else if (io->geo3d.fp != NULL) {
        fp = io->geo3d.fp;
        filename = io->geo3d.filename;
    } else {
        tl_error("FATAL ERROR :: 2d or 3d geo file not found!");
    }
    
    // loops over the lines in the input geometry file
    rewind(fp);
    while (fgets(line, MAXLINE, fp) != NULL) {
        io_save_line(io, fp, filename, line);
        if (strip_comments(line) <= 1) {    /* ignore empty ('\n') or comment lines */
            continue;
        }
        switch (parse_card(line, &data)) {
                /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
                /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
            case CARD_GRID:
                switch (parse_card(data, &subdata)) {
                    case CARD_NOCOL:
                        grid->type = UNSTRUCTURED;
                        break;
                    case CARD_COL:
                        grid->type = COLUMNAR;
                        break;
                    default:
                        grid->type = UNSTRUCTURED;
                        break;
                }
                break;
                /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
                /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
            case CARD_E4T:
            case CARD_TET:
                read_element(io, grid, JUST_COUNT, 3, NDONTET, data);
                break;
                /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
                /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
            case CARD_PRISM:
                read_element(io, grid, JUST_COUNT, 3, NDONPRISM, data);
                break;
                /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
                /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
            case CARD_E3T:
            case CARD_TRI:
                read_element(io, grid, JUST_COUNT, 2, NDONTRI, data);
                break;
                /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
                /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
            case CARD_QUAD:
                read_element(io, grid, JUST_COUNT, 2, NDONQUAD, data);
                break;
                /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
                /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
            case CARD_ND:
                grid->macro_nnodes++; // add to global node count
                
                if(JUST_COUNT==YES) break;
                
                id = read_int_field(*io, &data) - 1;
                xx = read_dbl_field(*io, &data);
                yy = read_dbl_field(*io, &data);
                zz = read_dbl_field(*io, &data);
                
                // create surface/bed node map
                if (grid->type == COLUMNAR) {
                    surf_id = UNSET_INT;
                    if (((fabs(fabs(xx) - fabs(surf_x)) > DIFF_TOL) || (fabs(fabs(yy) - fabs(surf_y)) > DIFF_TOL)) || global_surf_id == 0) {
                        surf_x=xx;
                        surf_y=yy;
                        surf_id = global_surf_id;
                        global_surf_id++;
                    }
                }
                
                // check to see if I own this node
                if (grid->smpi->npes > 1) {
                    for (i = 0; i < grid->smpi->npes; i++) {
                        if((id >= grid->start_node_ID[i]) && (id <= grid->end_node_ID[i])) owner = i;
                    }
                    if (owner == grid->smpi->myid) {
                        local_index = id - grid->start_node_ID[grid->smpi->myid];
                        resident_id = local_index;
                        if (local_index > grid->end_node_ID[grid->smpi->myid]) {
                            printf("myid: %d local_index: %d grid->start_node_ID[grid->smpi->myid]: %d grid->end_node_ID[grid->smpi->myid]: %d",grid->smpi->myid,local_index,grid->start_node_ID[grid->smpi->myid],grid->end_node_ID[grid->smpi->myid]);
                            tl_error("Nodal index appears to be out of range.");
                        }
                    } else if ((local_index = search_ghost_node(ghost_nodes, id)) >= 0) {
                        resident_id = id - grid->start_node_ID[owner];
                    } else {
                        break;
                    }
                } else {
                    local_index = id;
                    owner = 0;
                    resident_id = 0;
                }
                
                grid->node[local_index].id = local_index;
                grid->node[local_index].gid = id;
                grid->node[local_index].original_id = id;
                grid->node[local_index].x = xx;
                grid->node[local_index].y = yy;
                grid->node[local_index].z = zz;
                grid->node[local_index].myid = grid->smpi->myid;
                grid->node[local_index].resident_pe = owner;
                grid->node[local_index].global_surf_id = surf_id;
                grid->node[local_index].resident_id = resident_id;
                grid->node[local_index].string = NORMAL; /* dummy set up for the nodes */
                break;
            default:
                break;
        }
        io_save_line(io, NULL, "", "");
    }
    
#ifdef _DEBUG
    if (grid->smpi->myid == 0) {
        printf("\n---- Finished reading elements and nodes from geo file");
        if (JUST_COUNT == YES) printf(" (for counting purposes only)");
        printf("\n");
    }
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_element(SIO *io, SGRID *grid, int JUST_COUNT, int elem_dim, int elem_nnodes, char *data) {
    
    int i, j, k, ipe, id, ielem, iel_mat, flag, break_flag = FALSE;
    int nd[6]; sarray_init_value_int(nd,6,UNSET_INT);
    int pe_owner[6]; sarray_init_value_int(pe_owner,6,UNSET_INT);
    int local_id[6]; sarray_init_value_int(local_id,6,UNSET_INT);
    
    // add to global element count
    if (elem_dim == 2) {
        grid->macro_nelems2d++;
    } else {
        grid->macro_nelems3d++;
    }
    
    // read element global/macro ID
    id = read_int_field(*io, &data);
#ifdef _DEBUG
    if (grid->smpi->myid == 0 && DEBUG_FULL == ON) printf("mype: %d reading element %d \n",grid->smpi->myid,id);
#endif
    
    // read element nodes
    for (i=0;i<elem_nnodes;i++) {
        nd[i] = read_int_field(*io, &data);
    }
    
    // determine which pe owns each node of element
    for (i=0;i<elem_nnodes;i++) {
        for (ipe = 0; ipe < grid->smpi->npes; ipe++){
            if(((nd[i]-1) >= grid->start_node_ID[ipe]) && ((nd[i]-1) <= grid->end_node_ID[ipe])) pe_owner[i] = ipe;
        }
    }
    
    // if no nodes are resedential to my pe, then bail
    flag = 0;
    for (i=0;i<elem_nnodes;i++) {
        if (pe_owner[i] == grid->smpi->myid) flag = 1;;
    }
    if (flag == 0) return;
    
    // add to global element count
    if (elem_dim == 2) {
        ielem = grid->nelems2d;
        grid->nelems2d++;
    } else {
        ielem = grid->nelems3d;
        grid->nelems3d++;
    }
    
    // count and store ghost nodes on this element given a partition
    // we need ghost node count to allocated subdomain nnode (residential+ghost) arrays
    for (i = 0; i < elem_nnodes; i++) {
        if (pe_owner[i] == grid->smpi->myid) {
            local_id[i] = (nd[i]-1) - grid->start_node_ID[grid->smpi->myid];
        } else {
            local_id[i]=UNSET_INT;
            local_id[i] = search_ghost_node(ghost_nodes, (nd[i]-1));
            if(local_id[i]<0){
                if(JUST_COUNT==NO) {
                    tl_error("Ghost node being added after node count!\n");
                }
                local_id[i] = add_ghost_node(&ghost_nodes, (nd[i]-1), grid->my_nnodes, &num_ghosts);
            }
        }
    }
    
    // if we're just counting, we still need to allocation memorey before below, so call again after allocation
    if(JUST_COUNT == YES) return;
    
    // read element material ID
    iel_mat = read_int_field(*io, &data);
    if (iel_mat <= 0) tl_error("Element material value less than one.");
    if (iel_mat > grid->nmat) grid->nmat = iel_mat;
    
    if (elem_nnodes == NDONTET) {
        grid->elem3d[ielem].nnodes = elem_nnodes;
        selem3d_alloc(&(grid->elem3d[ielem]), elem_nnodes);
        grid->elem3d[ielem].nedges = 6;
        grid->elem3d[ielem].nnodes_quad = 10;
        grid->elem3d[ielem].edges = grid->nd_on_TetEdge;
        grid->elem3d[ielem].mat = iel_mat - 1;
        for (i=0;i<elem_nnodes;i++) grid->elem3d[ielem].nodes[i] = local_id[i];
        grid->elem3d[ielem].gid = id-1;
        grid->elem3d[ielem].id = grid->nelems3d;
        grid->haveTets = TRUE;
    } else if (elem_nnodes == NDONPRISM) {
        grid->elem3d[ielem].nnodes = elem_nnodes;
        selem3d_alloc(&(grid->elem3d[ielem]), elem_nnodes);
        grid->elem3d[ielem].nedges = 9;
        grid->elem3d[ielem].nnodes_quad = 15;
        grid->elem3d[ielem].edges = grid->nd_on_PrismEdge;
        grid->elem3d[ielem].mat = iel_mat - 1;
        for (i=0;i<elem_nnodes;i++) grid->elem3d[ielem].nodes[i] = local_id[i];
        grid->elem3d[ielem].gid = id-1;
        grid->elem3d[ielem].id = grid->nelems3d;
        grid->isGridMixedElement = TRUE;
        grid->havePrisms = TRUE;
    } else if (elem_nnodes == NDONTRI) {
        grid->elem2d[ielem].nnodes = elem_nnodes;
        selem2d_alloc(&(grid->elem2d[ielem]), elem_nnodes);
        grid->elem2d[ielem].nedges = 3;
        grid->elem2d[ielem].nnodes_quad = 6;
        grid->elem2d[ielem].edges = grid->nd_on_TriEdge;
        grid->elem2d[ielem].mat = iel_mat - 1;
        for (i=0;i<elem_nnodes;i++) grid->elem2d[ielem].nodes[i] = local_id[i];
        grid->elem2d[ielem].gid = id-1;
        grid->elem2d[ielem].id = grid->nelems2d;
        grid->haveTris = TRUE;
    } else if (elem_nnodes == NDONQUAD) {
        grid->elem2d[ielem].nnodes = elem_nnodes;
        selem2d_alloc(&(grid->elem2d[ielem]), elem_nnodes);
        grid->elem2d[ielem].nedges = 4;
        grid->elem2d[ielem].nnodes_quad = 8;
        grid->elem2d[ielem].edges = grid->nd_on_QuadEdge;
        grid->elem2d[ielem].mat = iel_mat - 1;
        for (i=0;i<elem_nnodes;i++) grid->elem2d[ielem].nodes[i] = local_id[i];
        grid->elem2d[ielem].gid = id-1;
        grid->elem2d[ielem].id = grid->nelems2d;
        grid->haveQuads = TRUE;
    } else {
        tl_error("ERROR: element not recognized!\n");
    }
    
#ifdef _DEBUG
    if (grid->smpi->myid == 0 && DEBUG_FULL == ON) printf("finished reading TET \n");
#endif
    
}
