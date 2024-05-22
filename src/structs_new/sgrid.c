/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  sgrid.c This file collects methods of the SGRID structure */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "global_header.h"
//#include "local_header.h"

// file scope variables
static int DEBUG = OFF;
static int DEBUG_WITH_PICKETS = OFF;


// CJT -- TAKE THIS OUT LATER
int doesFileExist(const char *fname) {
    FILE *file;
    if ((file = fopen(fname, "r")))
    {
        fclose(file);
        return 1;
    }
    return 0;
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
 * \note CJT: Grid file has nnodes, nelems, etc. on top of file
 * \note CJT: For speed, nodes should go first in the geo file
 * \note CJT: ND 0/1/2 (surface/bed/side) x y z
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void sgrid_read(SGRID **pgrid, char *root_filename
#ifdef _MESSG
                ,MPI_Comm model_comm
#endif
) {
    
    int i, j, k, id, idum, ipe;
    char *line = NULL, cdum[20], line_save[100];
    size_t len = 0;
    ssize_t read;
    char *token, *subtoken;
    char delim[] = " ";
    double fdum;
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // allocate the grid
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    char grid_filename[80];
    strcpy(grid_filename,root_filename);
    strcat(grid_filename, ".geo");
    (*pgrid) = (SGRID *) tl_alloc(sizeof(SGRID), 1);
    SGRID *g = *pgrid; // alias
    strcpy(g->filename,grid_filename);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // allocate MPI on this grid
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    g->smpi = (SMPI *) tl_alloc(sizeof(SMPI), 1);
    g->part_smpi = NULL;
    g->part_map = NULL;
    smpi_defaults(g->smpi);
#ifndef _MESSG
    smpi_init(g->smpi);
#else
    smpi_init(g->smpi, model_comm);
#endif
    int npes = g->smpi->npes; // alias
    int myid = g->smpi->myid; // alias
    //printf("npes: %d\n",npes);
    //printf("myid: %d\n",myid);
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // open the grid file
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef _DEBUG
    if (myid == 0) printf("\n-- Opening Grid File: %s\n",grid_filename);
#endif
    
    // Check to see if grid file exists
    FILE *fp = NULL;
    if (doesFileExist(grid_filename)) {
        fp = fopen(grid_filename, "r");
    } else {
        tl_error("ERROR: Geo file does not exist!\n");
    }
    
    // Does the grid have a columnar component?
    g->columnar = false;
    g->ndim = 2; // default to 2D grid
    while ((read = getline(&line, &len, fp)) != -1) {
        token = strtok(line,delim); token[strcspn(token, "\n")] = 0;
        if (strcmp(token, "GRID") == 0) {
            subtoken = strtok(NULL, delim);  subtoken[strcspn(subtoken, "\n")] = 0;
            if (strcmp(subtoken,"COLUMNAR") == 0) {
                strcpy(g->type,"COLUMNAR");
                g->columnar = true;
            }
        }
    }
    rewind(fp);
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Read grid nodes
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // NOTE: 2D nodes are also counted as surface and bed nodes
    g->macro_nnodes = 0;
    g->macro_nelems3d = 0;
    g->macro_nelems2d = 0;
    g->macro_nelems1d = 0;
    g->macro_nnodes_sur = 0;
    g->macro_nnodes_bed = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        strcpy(line_save,line);
        token = strtok(line,delim); token[strcspn(token, "\n")] = 0;
        if (strcmp(token, "NODES") == 0) {
            sscanf(line_save, "%s %d",cdum,&(g->macro_nnodes));
        } else if (strcmp(token, "NODES_SURFACE") == 0) {
            sscanf(line_save, "%s %d",cdum,&(g->macro_nnodes_sur));
        } else if (strcmp(token, "NODES_BED") == 0) {
            sscanf(line_save, "%s %d",cdum,&(g->macro_nnodes_bed));
        } else if (strcmp(token, "ELEMS1D") == 0) {
            sscanf(line_save, "%s %d",cdum,&(g->macro_nelems1d));
        } else if (strcmp(token, "ELEMS2D") == 0) {
            sscanf(line_save, "%s %d",cdum,&(g->macro_nelems2d));
        } else if (strcmp(token, "ELEMS3D") == 0) {
            sscanf(line_save, "%s %d",cdum,&(g->macro_nelems3d));
        }
    }
    assert(g->macro_nnodes > 0);
    assert(g->macro_nnodes_sur > 0);
    assert(g->macro_nnodes_bed > 0);
    assert(g->macro_nelems1d + g->macro_nelems2d + g->macro_nelems3d > 0);
    rewind(fp);
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Initialize SGRID components
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    g->mesh_length=0.0; g->my_mesh_length=0.0;
    g->mesh_area=0.0; g->my_mesh_area=0.0;
    g->mesh_area_surface=0.0; g->my_mesh_area_surface=0.0;
    g->mesh_area_bed=0.0; g->my_mesh_area_bed=0.0;
    g->mesh_area_sidewalls=0.0; g->my_mesh_area_sidewalls=0.0;
    g->mesh_volume=0.0; g->my_mesh_volume=0.0;
    g->x_min = 9.9e+99; g->x_max = -9.9e+99;
    g->y_min = 9.9e+99; g->y_max = -9.9e+99;
    g->z_min = 9.9e+99; g->z_max = -9.9e+99;
    g->my_x_min = 9.9e+99; g->my_x_max = -9.9e+99;
    g->my_y_min = 9.9e+99; g->my_y_max = -9.9e+99;
    g->my_z_min = 9.9e+99; g->my_z_max = -9.9e+99;
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Do a simple domain decomposition
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Find where start/stop nodes are for each PE
    int *nodes_per_pe = (int *) tl_alloc(sizeof(int), npes);
    int *start_node_id = (int *) tl_alloc(sizeof(int), npes);
    int *end_node_id   = (int *) tl_alloc(sizeof(int), npes);
    sarray_init_value_int(nodes_per_pe, npes, UNSET_INT);
    sarray_init_value_int(start_node_id, npes, UNSET_INT);
    sarray_init_value_int(end_node_id, npes, UNSET_INT);
    tl_check_all_pickets(__FILE__,__LINE__);
    
    if (npes > 1) {
#ifdef _MESSG
        
        // Even divide the total nodes among PEs
        int nodesOnPe = g->macro_nnodes / npes;
        start_node_id[0] = 0;
        for (ipe=1; ipe<npes; ipe++) {
            start_node_id[ipe] = ipe * nodesOnPe;
            end_node_id[ipe-1] = start_node_id[ipe] - 1;
        }
        end_node_id[npes-1] = g->macro_nnodes - 1;
        
        // Debug
        if (myid == 0) for (ipe=0; ipe<npes; ipe++) {
            printf("g->macro_nnodes: %d || pe: %d || start_node_id: %d || end_node_id: %d \n",
                   g->macro_nnodes,ipe,start_node_id[ipe],end_node_id[ipe]);
        }
        tl_check_all_pickets(__FILE__,__LINE__);
        messg_barrier(g->smpi->ADH_COMM);
        
        // Make sure a full column resides on a PE if columnar grid
        if (g->columnar) {
            int node_surf = UNSET_INT, last_node_surf = UNSET_INT;
            ipe = 0;
            while ((read = getline(&line, &len, fp)) != -1) {
                strcpy(line_save,line);
                token = strtok(line,delim); token[strcspn(token, "\n")] = 0;
                if (strcmp(token, "ND") == 0) {
                    sscanf(line_save, "%s %d %lf %lf %lf %d",cdum,&idum,&fdum,&fdum,&fdum,&node_surf);
                    if (idum - 1 > end_node_id[ipe] && node_surf != last_node_surf) {
                        end_node_id[ipe] = idum - 2;
                        start_node_id[ipe+1] = idum - 1;
                        ipe++;
                    }
                    last_node_surf = node_surf;
                }
            }
            rewind(fp);
        }
        
        // Nodes per PE after columnar grid changes
        for (ipe=0; ipe<npes; ipe++) {
            nodes_per_pe[ipe] = end_node_id[ipe] - start_node_id[ipe];
        }
        
        // Columnar grid variables
        if (g->columnar) {
            g->my_nnodes_sur = nodes_per_pe[myid]; // residential surface nodes
            g->my_nnodes_bed = nodes_per_pe[myid]; // residential bed nodes
        }
        
        // Grid residential nodes
        g->my_nnodes = end_node_id[myid] - start_node_id[myid] + 1; // residential nodes
#endif
    } else {
        
        start_node_id[0] = 0;
        end_node_id[0] = g->macro_nnodes - 1;
        
        g->nnodes = g->macro_nnodes;
        g->my_nnodes = g->macro_nnodes;
        g->my_nnodes_sur = g->macro_nnodes_sur;
        g->my_nnodes_bed = g->macro_nnodes_bed;
        g->nelems1d = g->macro_nelems1d;
        g->nelems2d = g->macro_nelems2d;
        g->nelems3d = g->macro_nelems3d;
        assert(g->nnodes > 0);
    }
    
    g->orig_macro_nnodes = g->macro_nnodes;
    g->orig_macro_nnodes_sur = g->macro_nnodes_sur;
    g->orig_macro_nnodes_bed = g->macro_nnodes_bed;
    g->orig_macro_nelems1d = g->macro_nelems1d;
    g->orig_macro_nelems2d = g->macro_nelems2d;
    g->orig_macro_nelems3d = g->macro_nelems3d;
    
//    printf("start_node_id[myid]: %d || end_node_id[myid]: %d || my_nnodes: %d || macro_nnodes: %d\n",start_node_id[myid],end_node_id[myid],g->my_nnodes,g->macro_nnodes);
//    MPI_Finalize();
//    exit(-1);
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Count the # of local ghost nodes and elements on my processor
    // for local allocation
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // allocate ghost node listing
    int num_ghosts = 0;
    NODE_LIST_ITEM *ghost_nodes = tl_alloc(sizeof(NODE_LIST_ITEM), 1);
    ghost_nodes->local=UNSET_INT;
    ghost_nodes->global=UNSET_INT;
    ghost_nodes->next = NULL;
    
    if (npes > 1) {
#ifdef _MESSG
        g->nelems1d = 0; g->nelems2d = 0; g->nelems3d = 0;
        while ((read = getline(&line, &len, fp)) != -1) {
            strcpy(line_save,line);
            token = strtok(line,delim); token[strcspn(token, "\n")] = 0;
            if        (strcmp(token, "SEG") == 0) {
                g->nelems1d += sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,2,1,true);
            } else if (strcmp(token, "TRI") == 0) {
                g->nelems2d += sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,3,2,true);
            } else if (strcmp(token, "QUAD") == 0) {
                g->nelems2d += sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,4,2,true);
            } else if (strcmp(token, "TET") == 0) {
                g->nelems3d += sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,4,3,true);
            } else if (strcmp(token, "PRISM") == 0) {
                g->nelems3d += sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,6,3,true);
            }
        }
        g->nnodes = g->my_nnodes + num_ghosts;
        rewind(fp);
#endif
    }
    assert(g->nnodes > 0);
    //MPI_Finalize();
    //exit(-1);
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Allocate grid
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // allocated nodes and elements and physics on each element
    if (g->nelems1d > 0) {selem1d_init_alloc_array(&(g->elem1d), g->nelems1d);}
    if (g->nelems2d > 0) {selem2d_init_alloc_array(&(g->elem2d), g->nelems2d);}
    if (g->nelems3d > 0) {selem3d_init_alloc_array(&(g->elem3d), g->nelems3d);}
    snode_init_alloc_array(&(g->node), g->nnodes);
    tl_check_all_pickets(__FILE__,__LINE__);
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Read through to store nodes first.  Need this to calculate element jacobians.
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    while ((read = getline(&line, &len, fp)) != -1) {
        strcpy(line_save,line);
        token = strtok(line,delim); token[strcspn(token, "\n")] = 0;
        if (strcmp(token, "ND") == 0) {
            sgrid_read_node(g,line_save,start_node_id,end_node_id,ghost_nodes);
        }
    }
    rewind(fp);
    
    char myfile[20];
    sprintf(myfile, "pe_%d_nodes.txt",myid);
    FILE *ft = fopen(myfile, "w");
    for (i=0; i<g->nnodes; i++) {
        fprintf(ft,"local_id: %d || global_id: %d || residential_pe: %d\n",g->node[i].id+1,g->node[i].gid+1,g->node[i].resident_pe);
    }
    //MPI_Finalize();
    //exit(-1);
    
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Read through again to store elements
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    g->nelems1d = 0; g->nelems2d = 0; g->nelems3d = 0;
    g->haveSegs = false; g->haveTets = false;
    g->havePrisms = false; g->haveQuads = false;
    while ((read = getline(&line, &len, fp)) != -1) {
        strcpy(line_save,line);
        token = strtok(line,delim); token[strcspn(token, "\n")] = 0;
        if        (strcmp(token, "SEG") == 0) {
            g->haveSegs = true;
            g->nelems1d += sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,2,1,false);
        } else if (strcmp(token, "TRI") == 0) {
            g->haveTris = true;
            g->nelems2d += sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,3,2,false);
        } else if (strcmp(token, "QUAD") == 0) {
            g->haveQuads = true;
            g->nelems2d += sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,4,2,false);
        } else if (strcmp(token, "TET") == 0) {
            g->haveTets = true;
            g->nelems3d += sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,4,3,false);
        } else if (strcmp(token, "PRISM") == 0) {
            g->havePrisms = true;
            g->nelems3d += sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,6,3,false);
        }
        //tl_check_all_pickets(__FILE__,__LINE__);
    }
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // find highest grid dimension
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    g->ndim = 1;
    if (g->haveTris || g->haveQuads) g->ndim = 2;
    if (g->havePrisms || g->haveTets) g->ndim = 3;\
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Close grid file
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fclose(fp);
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // free the ghost node list
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    free_node_list(ghost_nodes);
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // calculate quadrature points/weights
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int quad_flag;
    if (g->haveSegs)   quad_flag = squad_segment_alloc_init(&(g->quad_seg));
    if (g->haveTris)   quad_flag = squad_triangle_alloc_init(&(g->quad_tri));
    if (g->haveQuads)  quad_flag = squad_rectangle_alloc_init(&(g->quad_rect));
    if (g->haveTets)   quad_flag = squad_tetrahedron_alloc_init(&(g->quad_tet));
    if (g->havePrisms) quad_flag = squad_triprism_alloc_init(&(g->quad_prism));
    assert(quad_flag == 0);
    
    g->max_nnodes = g->nnodes;
    g->max_nelems1d = g->nelems1d;
    g->max_nelems2d = g->nelems2d;
    g->max_nelems3d = g->nelems3d;
    
    g->nnodes_old = g->nnodes;
    g->nelems1d_old = g->nelems1d;
    g->nelems2d_old = g->nelems2d;
    g->nelems3d_old = g->nelems3d;
    
    //g->my_nnodes_old = g->my_nnodes;
    g->max_nnodes_sur = g->my_nnodes_sur;
    g->max_nnodes_bed = g->my_nnodes_bed;
    g->nnodes_sur_old = g->my_nnodes_sur;
    g->nnodes_bed_old = g->my_nnodes_bed;
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // calculate total grid length, area, and/or volume over all PEs
    // note: this will ** NOT ** include non-zero initial surface/bed displacements
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (npes > 1) {
#ifdef _MESSG
        g->mesh_length = messg_dsum(g->my_mesh_length, g->smpi->ADH_COMM);
        g->mesh_area   = messg_dsum(g->my_mesh_area,   g->smpi->ADH_COMM);
        g->mesh_area_surface     = messg_dsum(g->my_mesh_area_surface,   g->smpi->ADH_COMM);
        g->mesh_area_bed         = messg_dsum(g->my_mesh_area_bed,       g->smpi->ADH_COMM);
        g->mesh_area_sidewalls   = messg_dsum(g->my_mesh_area_sidewalls, g->smpi->ADH_COMM);
        g->mesh_volume = messg_dsum(g->my_mesh_volume, g->smpi->ADH_COMM);
        g->x_min = messg_dsum(g->my_x_min, g->smpi->ADH_COMM);
        g->x_max = messg_dsum(g->my_x_max, g->smpi->ADH_COMM);
        g->y_min = messg_dsum(g->my_y_min, g->smpi->ADH_COMM);
        g->y_max = messg_dsum(g->my_y_max, g->smpi->ADH_COMM);
        g->z_min = messg_dsum(g->my_z_min, g->smpi->ADH_COMM);
        g->z_max = messg_dsum(g->my_z_max, g->smpi->ADH_COMM);
#endif
    } else {
        g->mesh_length = g->my_mesh_length;
        g->mesh_area = g->my_mesh_area;
        g->mesh_area_surface = g->my_mesh_area_surface;
        g->mesh_area_bed = g->my_mesh_area_bed;
        g->mesh_area_sidewalls = g->my_mesh_area_sidewalls;
        g->mesh_volume = g->my_mesh_volume;
        g->x_min = g->my_x_min; g->x_max = g->my_x_max;
        g->y_min = g->my_y_min; g->y_max = g->my_y_max;
        g->z_min = g->my_z_min; g->z_max = g->my_z_max;
    }
    
#ifdef _DEBUG
    if (myid == 0) printf("\n-- Finished reading geo file: %s\n",grid_filename);
#endif
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void sgrid_read_node(SGRID *g, char *line, int *start_node_id, int *end_node_id, NODE_LIST_ITEM *ghost_nodes) {
    
    int i, owner, local_index, resident_id, gid, column_flag;
    int npes = g->smpi->npes; // alias
    int myid = g->smpi->myid; // alias
    double x,y,z;
    char cdum[20];
    
    // read the global nodes IDs
    //printf("sgrid_read_node || line: %s \n",line);
    sscanf(line, "%s %d %lf %lf %lf %d",cdum,&gid,&x,&y,&z,&column_flag);
    gid--;
    
    
    // check to see if I own this node
    if (npes > 1) {
        local_index = search_ghost_node(ghost_nodes, gid); //printf("local_index: %d\n",local_index);
        for (i=0; i<npes; i++) {
            if((gid >= start_node_id[i]) && (gid <= end_node_id[i])) owner = i;
        }
        if (owner == myid) {
            local_index = gid - start_node_id[myid];
            resident_id = local_index;
            if (local_index > end_node_id[myid]) {tl_error("Nodal index appears to be out of range.");}
        } else if (local_index >= 0) {
            resident_id = gid - start_node_id[owner];
        } else {
            return;  // this node is not a part of my subdomain
        }
    } else {
        local_index = gid;
        owner = 0;
        resident_id = 0;
    }
    
    //printf("local_index: %d || npes: %d || gid: %d\n",local_index,npes,gid);
    g->node[local_index].id = local_index;
    g->node[local_index].gid = gid;
    g->node[local_index].original_id = gid;
    g->node[local_index].x = x;
    g->node[local_index].y = y;
    g->node[local_index].z = z;
    g->node[local_index].myid = myid;
    g->node[local_index].resident_pe = owner;
    g->node[local_index].global_surf_id = column_flag;
    g->node[local_index].resident_id = resident_id;
    g->node[local_index].string = NORMAL; /* dummy set up for the nodes */
    
    if (resident_id == myid) {
        if (x > g->my_x_max) g->my_x_max = x;
        if (x < g->my_x_min) g->my_x_min = x;
        if (y > g->my_y_max) g->my_y_max = y;
        if (y < g->my_y_min) g->my_y_min = y;
        if (z > g->my_z_max) g->my_z_max = z;
        if (z < g->my_z_min) g->my_z_min = z;
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// calculates nnodes, nelems1d/2d/3d and creates list of nodes
// CJT :: Must read nodes first before fully loading this!
int sgrid_read_elem(SGRID *g, char *line, int *start_node_id, int *end_node_id, int *num_ghosts, NODE_LIST_ITEM **ghost_nodes, int nnodes_on_elem, int elem_dim, bool JUST_COUNT) {
    
    //printf("sgrid_read_elem || line: %s || JUST COUNT: %d\n",line,JUST_COUNT);
    
    int i, bflag = UNSET_INT, id = UNSET_INT;
    int npes = g->smpi->npes; // alias
    int myid = g->smpi->myid; // alias
    SVECT nds[nnodes_on_elem]; svect_init_array(nds, nnodes_on_elem);
    int nd[nnodes_on_elem]; sarray_init_value_int(nd,nnodes_on_elem,UNSET_INT);
    int doIown[nnodes_on_elem]; sarray_init_value_int(doIown,nnodes_on_elem,NO);
    
    char cdum[20];
    if (nnodes_on_elem == 2) {
        sscanf(line, "%s %d %d %d %d",cdum,&id,&nd[0],&nd[1],&bflag);
    } else if (nnodes_on_elem == 3) {
        sscanf(line, "%s %d %d %d %d %d",cdum,&id,&nd[0],&nd[1],&nd[2],&bflag);
        //printf("element: %s || gid: %d || nodes: %d %d %d || bflag: %d\n",cdum,id,nd[0],nd[1],nd[2],bflag);
    } else if (nnodes_on_elem == 4) {
        sscanf(line, "%s %d %d %d %d %d %d",cdum,&id,&nd[0],&nd[1],&nd[2],&nd[3],&bflag);
        //printf("element: %s || gid: %d || nodes: %d %d %d %d || bflag: %d\n",cdum,id,nd[0],nd[1],nd[2],nd[3],bflag);
    } else if (nnodes_on_elem == 6) {
        sscanf(line, "%s %d %d %d %d %d %d %d %d",cdum,&id,&nd[0],&nd[1],&nd[2],&nd[3],&nd[4],&nd[5],&bflag);
    }
    
    id--;
    for (i = 0; i < nnodes_on_elem; i++) {nd[i]--;}
    
    int elem_resid_on_myid = false;
    int local_id[nnodes_on_elem]; sarray_init_value_int(local_id,nnodes_on_elem,UNSET_INT);
    if (npes > 1) {
        
        // determine which pe owns each node of element
        for (i=0; i<nnodes_on_elem; i++) {
            if((nd[i] >= start_node_id[myid]) && (nd[i] <= end_node_id[myid])) doIown[i] = YES;
        }
        
        // if no nodes are residential to my pe, then bail
        int flag = 0;
        for (i=0; i<nnodes_on_elem; i++) {
            if (doIown[i] == YES) flag += 1;
        }
        if (flag == 0) {
            return 0;           // do not add to local counts
        } else if (flag == nnodes_on_elem) {
            //if (myid==1 && (nnodes_on_elem == 4 || bflag == BODY)) printf("ielem || nnodes: %d || gid: %d\n",nnodes_on_elem, id);
            elem_resid_on_myid = true; // all nodes on this element are residential, so the element is
        }
        //if (myid == 1) printf("flag: %d\n",flag);
        
        // this element is part of myid
        if (JUST_COUNT) {
            // only add ghost nodes if the element is a body element
            if (nnodes_on_elem == 4 || nnodes_on_elem == 6 || bflag == BODY) {
                // add ghost nodes to list
                for (i = 0; i < nnodes_on_elem; i++) {
                    if (doIown[i] != YES) {
                        local_id[i] = UNSET_INT;
                        local_id[i] = search_ghost_node(*ghost_nodes,nd[i]);
                        if(local_id[i] < 0){
                            local_id[i] = add_ghost_node(ghost_nodes,nd[i], g->my_nnodes, num_ghosts);
                            //printf("ghost node being added to list || nd: %d || g->my_nnodes: %d || num_ghosts: %d || local_id: %d \n ",nd[i],g->my_nnodes,*num_ghosts,local_id[i] );
                        }
                    }
                }
            }
        } else {
            // find local node numbers
            for (i = 0; i < nnodes_on_elem; i++) {
                if (doIown[i] == YES) {
                    local_id[i] = nd[i] - start_node_id[myid];
                } else {
                    local_id[i] = search_ghost_node(*ghost_nodes, nd[i]);
                    if(local_id[i] < 0 && JUST_COUNT == true){tl_error("Ghost node being added after node count!\n");}
                }
                nds[i].x = g->node[local_id[i]].x;
                nds[i].y = g->node[local_id[i]].y;
                nds[i].z = g->node[local_id[i]].z;
                
                //printf("JUST COUNT: %d || element gid: %d || NODE[%d] || local_id: %d || position: %f %f %f || nd: %d || start_node_id: %d || end_node_id: %d || flag: %d\n",JUST_COUNT, id+1,i,local_id[i]+1,nds[i].x,nds[i].y,nds[i].z,nd[i],start_node_id[myid],end_node_id[myid],flag);
            }
        }
        
    } else {
        elem_resid_on_myid = true;
        for (i = 0; i < nnodes_on_elem; i++) {
            local_id[i] = nd[i];
            nds[i].x = g->node[nd[i]].x;
            nds[i].y = g->node[nd[i]].y;
            nds[i].z = g->node[nd[i]].z;
        }
    }
    
    // if we're just counting for allocation, bail here
    if(JUST_COUNT) return 1;
    
    // if the element has some ghost nodes, arbitrarily assing its resident element PE to the lowest resident node PE
    int lowest_pe = 99999999;
    for (i=0; i<nnodes_on_elem; i++) {
        if (g->node[local_id[i]].resident_pe < lowest_pe) lowest_pe = g->node[local_id[i]].resident_pe;
    }
    if (myid == lowest_pe) elem_resid_on_myid = true;
    
    // load elements into AdH elemental array and add to subdomain mesh volume if the element is a body element
    if (elem_dim == 1)      {
        selem1d_init(&g->elem1d[g->nelems1d]);
        selem1d_load(&(g->elem1d[g->nelems1d]),id,g->nelems1d,nnodes_on_elem,local_id,bflag,nds);
        if ( (g->elem1d[g->nelems1d].bflag == BODY) && elem_resid_on_myid) {
            g->my_mesh_length += g->elem1d[g->nelems1d].length;
        }
    }
    else if (elem_dim == 2) {
        selem2d_init(&g->elem2d[g->nelems2d]);
        selem2d_load(&(g->elem2d[g->nelems2d]),id,g->nelems2d,nnodes_on_elem,local_id,bflag,nds);
        if (elem_resid_on_myid) g->elem2d[g->nelems2d].resident_pe = myid;
        if ( (g->elem2d[g->nelems2d].bflag == BODY) && elem_resid_on_myid) {
            g->my_mesh_area += g->elem2d[g->nelems2d].area;
        } else if ( (g->elem2d[g->nelems2d].bflag == SURFACE) && elem_resid_on_myid) {
            g->my_mesh_area_surface += g->elem2d[g->nelems2d].area;
        } else if ( (g->elem2d[g->nelems2d].bflag == BED) && elem_resid_on_myid) {
            g->my_mesh_area_bed += g->elem2d[g->nelems2d].area;
        } else if ( (g->elem2d[g->nelems2d].bflag == SIDEWALL) && elem_resid_on_myid) {
            g->my_mesh_area_sidewalls += g->elem2d[g->nelems2d].area;
        }
    }
    else if (elem_dim == 3) {
        selem3d_init(&g->elem3d[g->nelems3d]);
        selem3d_load(&(g->elem3d[g->nelems3d]),id,g->nelems3d,nnodes_on_elem,local_id,bflag,nds);
        if (elem_resid_on_myid) {
            g->elem3d[g->nelems3d].resident_pe = myid;
            g->my_mesh_volume += g->elem3d[g->nelems3d].volume;
        }
    }
    else {tl_error("ERROR: element not recognized!\n");}
    
    return 1;
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
void sgrid_free(SGRID *g) {
    int ie, inode;
    
    if (g->elem1d != NULL && g->nelems1d > 0) {selem1d_free_array(g->elem1d, g->max_nelems1d);}
    if (g->elem2d != NULL && g->nelems2d > 0) {selem2d_free_array(g->elem2d, g->max_nelems2d);}
    if (g->elem3d != NULL && g->nelems3d > 0) {selem3d_free_array(g->elem3d, g->max_nelems3d);}
    if (g->nnodes > 0) {
        for (inode=0; inode<g->max_nnodes; inode++) {snode_free(&(g->node[inode]));}
        g->node = (SNODE *) tl_free(sizeof(SNODE), g->max_nnodes, g->node);
    }
    g = (SGRID *) tl_free(sizeof(SGRID), 1, g);
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
                printf("------------ PE: %d || nnodes: %d || max_nnodes: %d || nnodes_old: %d\n",i,buffer1[i],buffer2[i],buffer3[i]);
            }
        }
        
        // print local nnode values in PE order
        sarray_init_int(buffer1, UNSET_INT); messg_gather_int(0, &g->my_nnodes_sur,  buffer1, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer2, UNSET_INT); messg_gather_int(0, &g->max_nnodes_sur, buffer2, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer3, UNSET_INT); messg_gather_int(0, &g->nnodes_sur_old, buffer3, 1, g->smpi->ADH_COMM);
        if (myid == 0) {
            for (i=0; i<npes; i++) {
                printf("------------ PE: %d || surface || nnodes: %d || max_nnodes: %d || nnodes_old: %d\n",i,buffer1[i],buffer2[i],buffer3[i]);
            }
        }
        
        // print local nnode values in PE order
        sarray_init_int(buffer1, UNSET_INT); messg_gather_int(0, &g->my_nnodes_bed,  buffer1, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer2, UNSET_INT); messg_gather_int(0, &g->max_nnodes_bed, buffer2, 1, g->smpi->ADH_COMM);
        sarray_init_int(buffer3, UNSET_INT); messg_gather_int(0, &g->nnodes_bed_old, buffer3, 1, g->smpi->ADH_COMM);
        if (myid == 0) {
            for (i=0; i<npes; i++) {
                printf("------------ PE: %d || bed || nnodes: %d || max_nnodes: %d || nnodes_old: %d\n",i,buffer1[i],buffer2[i],buffer3[i]);
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
