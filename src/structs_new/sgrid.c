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
    //Mark added
    g->macro_nTets = 0;
    g->macro_nPrisms = 0;
    g->macro_nQuads = 0;
    g->macro_nTris = 0; 
    while ((read = getline(&line, &len, fp)) != -1) {
        strcpy(line_save,line);
        token = strtok(line,delim); token[strcspn(token, "\n")] = 0;
        if (strcmp(token, "NODES") == 0) {
            sscanf(line_save, "%s %d",cdum,&(g->macro_nnodes));
        } else if (strcmp(token, "NODES_SURFACE") == 0) {
            sscanf(line_save, "%s %d",cdum,&(g->macro_nnodes_sur));
        } else if (strcmp(token, "NODES_BED") == 0) {
            sscanf(line_save, "%s %d",cdum,&(g->macro_nnodes_bed));
        } else if (strcmp(token, "ELEMS_SEG") == 0) {
            sscanf(line_save, "%s %d",cdum,&(g->macro_nelems1d));
        } else if (strcmp(token, "ELEMS_TRI") == 0) {
            sscanf(line_save, "%s %d",cdum,&(g->macro_nTris));
        } else if (strcmp(token, "ELEMS_QUAD") == 0) {
            sscanf(line_save, "%s %d",cdum,&(g->macro_nQuads));
        } else if (strcmp(token, "ELEMS_TET") == 0) {
            sscanf(line_save, "%s %d",cdum,&(g->macro_nTets));
        } else if (strcmp(token, "ELEMS_PRISM") == 0) {
            sscanf(line_save, "%s %d",cdum,&(g->macro_nQuads));
        }
    }

    //use element types to compute totals
    g->macro_nelems2d = g->macro_nTris + g->macro_nQuads;
    g->macro_nelems3d = g->macro_nTets + g->macro_nPrisms;

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
        
        // Grid residential nodes
        g->my_nnodes = end_node_id[myid] - start_node_id[myid] + 1; // residential nodes
#endif
    } else {
        
        start_node_id[0] = 0;
        end_node_id[0] = g->macro_nnodes - 1;
        
        g->nnodes = g->macro_nnodes;
        g->my_nnodes = g->macro_nnodes;
        g->nelems1d = g->macro_nelems1d;
        g->nelems2d = g->macro_nTris + g->macro_nQuads;
        g->nelems3d = g->macro_nTets + g->macro_nPrisms;
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
        //mark also add
        g->nTris=0; g->nQuads=0; g->nTets=0; g->nPrisms=0;

        while ((read = getline(&line, &len, fp)) != -1) {
            strcpy(line_save,line);
            token = strtok(line,delim); token[strcspn(token, "\n")] = 0;
            if        (strcmp(token, "SEG") == 0) {
                g->nelems1d += sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,2,1,true);
            } else if (strcmp(token, "TRI") == 0) {
                g->nTris += sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,3,2,true);
            } else if (strcmp(token, "QUAD") == 0) {
                g->nQuads += sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,4,2,true);
            } else if (strcmp(token, "TET") == 0) {
                g->nTets += sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,4,3,true);
            } else if (strcmp(token, "PRISM") == 0) {
                g->nPrisms += sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,6,3,true);
            }
        }
        //Mark add at end to get totals
        g->nelems2d = g->nTris + g->nQuads;
        g->nelems3d = g->nTets + g->nPrisms;
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
    g->my_nelems1d=0; g->my_nelems2d=0; g->my_nelems3d=0;
    g->my_nTris=0; g->my_nQuads=0; g->my_nTets=0; g->my_nPrisms=0;
    g->haveSegs = false; g->haveTets = false;
    g->havePrisms = false; g->haveQuads = false;
    int temp=0;
    while ((read = getline(&line, &len, fp)) != -1) {
        strcpy(line_save,line);
        token = strtok(line,delim); token[strcspn(token, "\n")] = 0;
        if (strcmp(token, "SEG") == 0) {
            g->haveSegs = true;
            temp = sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,2,1,false);
            g->nelems1d += temp;
            //Mark fix later
            if(false){g->my_nelems1d+=temp;}
        } else if (strcmp(token, "TRI") == 0) {
            g->haveTris = true;
            temp = sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,3,2,false);
            g->nelems2d += temp;
            if(g->elem2d[g->nelems2d-1].resident_pe == myid){g->my_nTris+=temp;}
        } else if (strcmp(token, "QUAD") == 0) {
            g->haveQuads = true;
            temp = sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,4,2,false);
            g->nelems2d += temp;
            if(g->elem2d[g->nelems2d-1].resident_pe == myid){g->my_nQuads+=temp;}
        } else if (strcmp(token, "TET") == 0) {
            g->haveTets = true;
            temp = sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,4,3,false);
            g->nelems3d += temp;
            if(g->elem3d[g->nelems3d-1].resident_pe == myid){g->my_nTets+=temp;}
        } else if (strcmp(token, "PRISM") == 0) {
            g->havePrisms = true;
            temp = sgrid_read_elem(g,line_save,start_node_id,end_node_id,&num_ghosts,&ghost_nodes,6,3,false);
            g->nelems3d += temp;
            if(g->elem3d[g->nelems3d-1].resident_pe == myid){g->my_nPrisms+=temp;}
        }
        //tl_check_all_pickets(__FILE__,__LINE__);
    }
    

    //Mark add at end to get totals
    g->my_nelems2d = g->my_nTris + g->my_nQuads;
    g->my_nelems3d = g->my_nTets + g->my_nPrisms;
    
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
    
    g->my_nnodes_sur = 0; g->nnodes_sur = 0;
    g->my_nnodes_bed = 0; g->nnodes_bed = 0;
    for (i=0; i<g->nnodes; i++) {
        if (myid==0) printf("PE: %d || bflag: %d || gid: %d || node PE: %d\n",myid,g->node[i].bflag,g->node[i].gid,g->node[i].resident_pe);
        if (g->node[i].bflag == SURFACE) {
            g->nnodes_sur++;
            if (g->node[i].resident_pe == myid) {
                g->my_nnodes_sur++;
            }
        } else if (g->node[i].bflag == BED) {
            g->nnodes_bed++;
            if (g->node[i].resident_pe == myid) {
                g->my_nnodes_bed++;
            }
        } else if (g->node[i].bflag == BODY2D) {
            g->nnodes_sur++;
            g->nnodes_bed++;
            if (g->node[i].resident_pe == myid) {
                g->my_nnodes_sur++;
                g->my_nnodes_bed++;
            }
        }
    }
    g->max_nnodes_sur = g->nnodes_sur;
    g->max_nnodes_bed = g->nnodes_bed;
    g->nnodes_sur_old = g->nnodes_sur;
    g->nnodes_bed_old = g->nnodes_bed;
    printf("PE: %d || g->my_nnodes_sur: %d || g->my_nnodes_bed: %d \n",myid,g->my_nnodes_sur,g->my_nnodes_bed);
    
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
        //sums go hear
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
    column_flag--;
    
    
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

    // CJT: Some of these will be ambiguous.  Surface/Bed always over-ride
    if(elem_dim == 2){
        if(bflag == BODY){
            for (i = 0; i < nnodes_on_elem; i++){
                if (g->node[local_id[i]].bflag == UNSET_INT) { // CJT: if this node is already defined, leave it alone
                    g->node[local_id[i]].bflag = BODY2D;
                }
            }
        }else if(bflag == SURFACE){
            for (i = 0; i < nnodes_on_elem; i++){g->node[local_id[i]].bflag = SURFACE;} // overide previous bflag
        }else if(bflag == BED){
            for (i = 0; i < nnodes_on_elem; i++){g->node[local_id[i]].bflag = BED;}     // overide previous bflag
        }else if(bflag == SIDEWALL){
            for (i = 0; i < nnodes_on_elem; i++){
                if (g->node[local_id[i]].bflag == UNSET_INT) { // CJT: if this node is already defined, leave it alone
                    g->node[local_id[i]].bflag = SIDEWALL;
                }
            }
        }
    }

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

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Print an AdH grid to XDMF
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] g           (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define MATRIX_RANK 2 // Tensor dimension :: Always be 2
#define VECTOR_RANK 1 // Always be 1
#define SPATIAL_DIM 3 // Always be 3



void sgrid_write_hdf5(SGRID *g){



#ifdef _ADH_HDF5

    hid_t     file_id, plist_id, filespace, dset_id, memspace;
    hid_t     grp1,grp2;
    char      fname[50];
    hsize_t   dims[MATRIX_RANK];
    herr_t    status;
    hsize_t   count[MATRIX_RANK],offset[MATRIX_RANK];
    int nt=0;


    // for adaptive grid writing
    char number[50] = "";
    sprintf(number, "%d", nt);

    // mesh properties
    int *connectivity;

    // create object to store h5 config
    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    #ifdef _MESSG
        MPI_Info info  = MPI_INFO_NULL;
        H5Pset_fapl_mpio(plist_id, g->smpi->ADH_COMM, info);
    #endif


    // open file
    strcpy(fname,g->filename);
    strcat(fname, ".h5");
    file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);

    // group 1 is anything with mesh nodes
    grp1 = H5Gopen(file_id, "/Mesh/XY", H5P_DEFAULT); // a group is like a folder

    //++++++++++++++++++++++++++++++++++++++
    // Writing Nodes
    //++++++++++++++++++++++++++++++++++++++

    // declare global sizes of data set nodes first. This is local data.
    double (*xyz)[SPATIAL_DIM] = malloc(sizeof(double[g->my_nnodes][SPATIAL_DIM]));




    // create a dataspace.  This is a global quantity
    dims[0]  = g->macro_nnodes;
    dims[1]  = SPATIAL_DIM;
    filespace = H5Screate_simple(MATRIX_RANK, dims, NULL);

    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp1, number, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    //create a local dataspace where each process has a few of the rows
    count[0]  = g->my_nnodes;
    count[1]  = dims[1];
    //offset[0] = g->node[0].gid;
    //offset[1] = 0;

    // give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(MATRIX_RANK, count, NULL);

    // determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    
    // the 2 is dimension of our data structure, RANK is spatial dimension of mesh
    hsize_t   coord1[g->my_nnodes * SPATIAL_DIM * MATRIX_RANK];
    
    // the values in coord1 should be the global node numbers, this time it is trivial
    int i,j,k=0,l=0;
    for (i=0; i<g->nnodes; i++){
        //temporarily store xyz if it is a residential node
        if(g->node[i].resident_pe == g->smpi->myid){
            xyz[l][0] = g->node[i].x; 
            xyz[l][1] = g->node[i].y;
            xyz[l][2] = g->node[i].z;
            l+=1;

        for(j=0;j<SPATIAL_DIM;j++) {
            coord1[k] = g->node[i].gid;
            coord1[k+1] = j;
            k += MATRIX_RANK;
        }

    }
    }

    assert(l==g->my_nnodes);
    H5Sselect_elements(filespace, H5S_SELECT_SET,g->my_nnodes * SPATIAL_DIM,(const hsize_t *)&coord1);
    
    // Create property list for collective dataset write.
    // declare collextive data file writing

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #ifdef _MESSG
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    #endif
    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xyz);
    free(xyz);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp1);
    
    // +++++++++++++++++++++++++++++++++++++++
    // Nodal write complete
    // +++++++++++++++++++++++++++++++++++++++


    // +++++++++++++++++++++++++++++++++++++++
    // Now write elemental data
    // +++++++++++++++++++++++++++++++++++++++

    // Connectivity data
    //size of array should be numQuad*(4+1) + numTri*(3+1) + numTet*(4+1) + numPrism*(6+1)
    //+1 is for element code
    grp2 = H5Gopen(file_id, "/Mesh/Elements", H5P_DEFAULT);

    //global attributes
    int nentry = g->macro_nQuads * 5 + g->macro_nTris * 4 + g->macro_nTets * 5 + g->macro_nPrisms * 7;

    //local data, different on each PE
    int nentry_local = g->my_nQuads * 5 + g->my_nTris * 4 + g->my_nTets * 5 + g->my_nPrisms * 7;
    connectivity = (int *)malloc(sizeof(int) * nentry_local);
    count[0] = nentry_local; // local dim[0]
    count[1] = 0; // ignores for elemental since just dim 1 array
    

    int ie;
    k = 0;
    for (ie=0; ie<g->nelems3d; ie++) {
        if (g->elem3d[ie].resident_pe == g->smpi->myid) {

            if      (g->elem3d[ie].nnodes == 4) {connectivity[k] = 6;} // XDMF CODE
            else if (g->elem3d[ie].nnodes == 6) {connectivity[k] = 7;} // XDMF CODE
            k += 1;
            for (i=0; i<g->elem3d[ie].nnodes; i++) {
                connectivity[k] = g->node[g->elem3d[ie].nodes[i]].gid;
                k++;
            }
        }
    }



    for (ie=0; ie<g->nelems2d; ie++) {
        //if (g->elem2d[ie].bflag == BODY) { // only body elements??  // MAKE SURE THEY ARE RESIDENTIAL!!!
        if (g->elem2d[ie].resident_pe == g->smpi->myid) {
            if      (g->elem2d[ie].nnodes == 3) {connectivity[k] = 4;} // XDMF CODE
            else if (g->elem2d[ie].nnodes == 4) {connectivity[k] = 5;} // XDMF CODE
            k += 1;
            for (i=0; i<g->elem2d[ie].nnodes; i++) {
                connectivity[k] = g->node[g->elem2d[ie].nodes[i]].gid;
                k++;
            }
        }
    }



    //add 1d element too
     for (ie=0; ie<g->nelems1d; ie++) {
        //if (g->elem2d[ie].bflag == BODY) { // only body elements??  // MAKE SURE THEY ARE RESIDENTIAL!!!
        //should be g->elem1d[ie].resident_pe == g->smpi->myid
        if (true) {
            if      (g->elem1d[ie].nnodes == 2) {connectivity[k] = 2;} // XDMF CODE
            k += 1;
            for (i=0; i<g->elem1d[ie].nnodes; i++) {
                connectivity[k] = g->node[g->elem1d[ie].nodes[i]].gid;
                k++;
            }
        }
    }


    assert(k == nentry_local);

    // create dataspace and add
    //give global size
    dims[0] = nentry;
    filespace = H5Screate_simple(VECTOR_RANK, dims, NULL);
    
    // create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp2, number, H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);
    
    // give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(VECTOR_RANK, count, NULL);

    //determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    //H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL); // only contiguous arrays
    
    // select elements
    // to do this we need to loop through each element and find out what is before it
    // the 2 is dimension of our data structure, RANK is spatial dimension of mesh
    hsize_t   coord2[nentry_local];
    
    // global element id number always will follow this convention:
    // 3d -> 2d -> 1d (descending order in terms of # nodes per element)
    // and prisms -> Tets
    // Quads -> Tris
    k=0;
    int gid,nnodes;
    for (ie=0; ie<g->nelems3d; ie++){

        if (g->elem3d[ie].resident_pe == g->smpi->myid) {
        gid = g->elem3d[ie].gid;
        nnodes = g->elem3d[ie].nnodes;
        //use gid to calculate global coord numbers

        //Prisms first
        if(nnodes == 6) {
            for(j=0;j<nnodes+1;j++){
                coord2[k] = gid*(nnodes+1) + j;
                k+=1;
            }
        //Tets second
        }else if(nnodes == 4){
            for(j=0;j<nnodes+1;j++){
                coord2[k] = g->macro_nPrisms*(7) + (gid-g->macro_nPrisms)*(nnodes+1) + j;
                k+=1;
            }

        }
        } 

    }



    for (ie=0; ie<g->nelems2d; ie++){
        gid = g->elem2d[ie].gid;
        nnodes = g->elem2d[ie].nnodes;
        //use gid to calculate global coord numbers
        if (g->elem2d[ie].resident_pe == g->smpi->myid) {
        //Quads first
        if(nnodes == 4) {
            for(j=0;j<nnodes+1;j++){
                coord2[k] = g->macro_nPrisms*7 + g->macro_nTets*5 + (gid)*(nnodes+1) + j;
                k+=1;
            }
        //Tris second
        }else if(nnodes == 3){
            for(j=0;j<nnodes+1;j++){
                coord2[k] = g->macro_nPrisms*7 + g->macro_nTets*5 + g->macro_nQuads*5 + (gid-g->macro_nQuads)*(nnodes+1) + j;
                k+=1;
            }

        } 
        }

    }

    //add 1d here
    for (ie=0; ie<g->nelems1d; ie++){
        gid = g->elem1d[ie].gid;
        nnodes = g->elem1d[ie].nnodes;
        //use gid to calculate global coord numbers
        //RN no mpi infor for 1d elements, need to fix
        //should be g->elem1d[ie].resident_pe == g->smpi->myid
        if (true) {
        //line seg only
        if(nnodes == 2) {
            for(j=0;j<nnodes+1;j++){
                coord2[k] = g->macro_nPrisms*7 + g->macro_nTets*5 + g->macro_nQuads*5 + g->macro_nTris*4 + (gid)*(nnodes+1) + j;
                k+=1;
            }
        //Tris second
        }
        }

    }

    //put into appropriate part of global array
    H5Sselect_elements(filespace, H5S_SELECT_SET,nentry_local,(const hsize_t *)&coord2);


    //Create property list for collective dataset write.
    //declare collextive data file weiting
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #ifdef _MESSG
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    #endif    
    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, connectivity);
    free(connectivity);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp2);

    ////////////////////////////////////
    ///Elemental connections complete///

    //close mesh file
    H5Fclose(file_id);
    
#endif

}

void init_hdf5_file(SGRID *g)
{
    #ifdef _ADH_HDF5
    //this function should do the following
    // 1. Create HDf5 File for the simulation
    // 2. Create all necessary data groups for the run

    

    hid_t     file_id, dset_id, memspace, filespace, plist_id;
    hid_t     grp1,grp2,grp3,grp4,grp5,grp6,grp7,grp8,grp9;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    //setting options for parallel

    #ifdef _MESSG
    MPI_Info info  = MPI_INFO_NULL;
    H5Pset_fapl_mpio(plist_id, g->smpi->ADH_COMM, info);
    #endif
    //H5Pset_all_coll_metadata_ops(plist_id, true);
    //H5Pset_coll_metadata_write(plist_id, true);


    //create file
    char fname[50];
    strcpy(fname,g->filename);
    strcat(fname, ".h5");
    file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);


    //1st 2 groups always active even if data is empty

    //group 1 is anything with mesh
    grp1 = H5Gcreate(file_id, "/Mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


    //group 2 contains any info with data
    grp2 = H5Gcreate(file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //Eventually want a routine to handle this based on active models but for
    

    //now just hard code some exmamples
    grp3 = H5Gcreate(grp2, "Time", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    grp4 = H5Gcreate(grp2, "NodalScalar", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    grp5 = H5Gcreate(grp2, "ElementalScalar", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    grp6 = H5Gcreate(grp2, "NodalVector", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    grp7 = H5Gcreate(grp2, "ElementalVector", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //these groups on the other hand will always be there so this can be hard coded
    grp8 = H5Gcreate(grp1, "XY", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    grp9 = H5Gcreate(grp1, "Elements", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);



    H5Gclose(grp1);
    H5Gclose(grp2);
    H5Gclose(grp3);
    H5Gclose(grp4);
    H5Gclose(grp5);
    H5Gclose(grp6);
    H5Gclose(grp7);
    H5Gclose(grp8);
    H5Gclose(grp9);
    H5Fclose(file_id);

    #endif

}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Print an AdH grid to XDMF
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] g           (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void sgrid_write_xdmf(SGRID *g){

if(g->smpi->myid==0){
    FILE *xmf = 0;
    char fname[50];
    float t=0;
    int mesh_no=0;
    strcpy(fname, g->filename);
    strcat(fname, ".xmf");
    xmf = fopen(fname, "w");
    
    //call from xdmf utils
    write_xdmf_header(xmf);
    
    fprintf(xmf, "\t\t\t<Grid Name=\"2D Unstructured Mesh\">\n");
            fprintf(xmf, "\t\t\t<Time Value=\"%f\" />\n",t);
            //Geometry start, this is nodes
            fprintf(xmf, "\t\t\t\t<Geometry GeometryType=\"XYZ\">\n");
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d %d\" Format=\"HDF\">\n", g->macro_nnodes,SPATIAL_DIM);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Mesh/XY/%d\n",g->filename,mesh_no);
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Geometry>\n");
            //Geometry finish, this is nodes

            //Topology start
            fprintf(xmf, "\t\t\t\t<Topology TopologyType=\"Mixed\" NumberOfElements=\"%d\">\n", g->macro_nelems1d+ g->macro_nelems2d+g->macro_nelems3d);
            //Data Item is the mixed connectivity
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", g->macro_nQuads * 5 + g->macro_nTris * 4 + g->macro_nTets * 5 + g->macro_nPrisms * 7);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Mesh/Elements/%d\n",g->filename,mesh_no);    
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Topology>\n");
            //Topology finish, this contains connectivities


            

    //call from xdmf utils        
    write_xdmf_tail(xmf);
    //XDMF File finished
    
    fclose(xmf);

}


}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Print an AdH grid to XDMF
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] g           (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void sgrid_write_xdmf_nodal_pe(SGRID *g){

if(g->smpi->myid==0){
    FILE *xmf = 0;
    char fname[50];
    float t=0;
    int mesh_no=0;
    strcpy(fname, g->filename);
    strcat(fname, ".xmf");
    xmf = fopen(fname, "r+");


  int lines = 0;
  int nchar = 0,nline=0;
  char line[256],c;

  //get line count first
  // Extract characters from file and store in character c
  for (c = getc(xmf); c != EOF; c = getc(xmf))
        if (c == '\n') // Increment count if this character is newline
            nline+= 1;
  //printf("Nline,%d\n",nline);
  rewind(xmf);
  // Count lines and move to the second last line (considering 0-based indexing)
  while (fgets(line, sizeof(line), xmf) != NULL) {
    lines++;
    if(lines==nline-4){break;}
    }
    //get the pointer
    fseek(xmf, 0, SEEK_CUR);

  // Insert text before the current position
  //Any nodal data we can also write
            fprintf(xmf, "\t\t\t\t<Attribute Name=\"NodalData\" AttributeType=\"Scalar\" Center=\"Node\">\n");
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", g->macro_nnodes);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Data/NodalScalar/%d\n",g->filename,mesh_no);
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Attribute>\n");

  // append file tail
  write_xdmf_tail(xmf);



    //XDMF File finished
    
    fclose(xmf);
    //printf("Text inserted successfully!\n");

}


}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Print an AdH grid to XDMF
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] g           (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void sgrid_write_xdmf_elemental_pe(SGRID *g){

if(g->smpi->myid==0){
    FILE *xmf = 0;
    char fname[50];
    float t=0;
    int mesh_no=0;
    strcpy(fname, g->filename);
    strcat(fname, ".xmf");
    xmf = fopen(fname, "r+");


  int lines = 0;
  int nchar = 0,nline=0;
  char line[256],c;

  //get line count first
  // Extract characters from file and store in character c
  for (c = getc(xmf); c != EOF; c = getc(xmf))
        if (c == '\n') // Increment count if this character is newline
            nline+= 1;
  //printf("Nline,%d\n",nline);
  rewind(xmf);
  // Count lines and move to the second last line (considering 0-based indexing)
  while (fgets(line, sizeof(line), xmf) != NULL) {
    lines++;
    if(lines==nline-4){break;}
    }
    //get the pointer
    fseek(xmf, 0, SEEK_CUR);
  

  // Insert text before the current position
  //Any nodal data we can also write
  //Should also be way to write elemental data here too just change Center=\"Node" to Center=\"Cell"
            fprintf(xmf, "\t\t\t\t<Attribute Name=\"ElementalData\" AttributeType=\"Scalar\" Center=\"Cell\">\n");
            fprintf(xmf, "\t\t\t\t\t<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", g->macro_nelems1d+ g->macro_nelems2d+g->macro_nelems3d);
            fprintf(xmf, "\t\t\t\t\t\t%s.h5:/Data/ElementalScalar/%d\n",g->filename,mesh_no);
            fprintf(xmf, "\t\t\t\t\t</DataItem>\n");
            fprintf(xmf, "\t\t\t\t</Attribute>\n");

  // append file tail
  write_xdmf_tail(xmf);

    //XDMF File finished
    
    fclose(xmf);
    

}


}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Print an AdH grid to XDMF
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] g           (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void sgrid_write_nodal_pe(SGRID *g){



#ifdef _ADH_HDF5

    hid_t     file_id, plist_id, filespace, dset_id, memspace;
    hid_t     grp1,grp2;
    char      fname[50];
    hsize_t   dims[MATRIX_RANK];
    herr_t    status;
    hsize_t   count[MATRIX_RANK],offset[MATRIX_RANK];
    int nt=0;


    // for adaptive grid writing
    char number[50] = "";
    sprintf(number, "%d", nt);

    // mesh properties
    int *connectivity;

    // create object to store h5 config
    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    #ifdef _MESSG
        MPI_Info info  = MPI_INFO_NULL;
        H5Pset_fapl_mpio(plist_id, g->smpi->ADH_COMM, info);
    #endif


    // open file
    strcpy(fname,g->filename);
    strcat(fname, ".h5");
    file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);

    // group 1 is anything with mesh nodes
    grp1 = H5Gopen(file_id, "/Data/NodalScalar", H5P_DEFAULT); // a group is like a folder

    //++++++++++++++++++++++++++++++++++++++
    // Writing Nodes
    //++++++++++++++++++++++++++++++++++++++

    // declare global sizes of data set nodes first. This is local data.
    //Nodally based scalar data, maybe outsource this to other routine later

    int *scalardata; 
    scalardata = (int *)malloc(sizeof(int) * g->my_nnodes);




    // create a dataspace.  This is a global quantity
    dims[0]  = g->macro_nnodes;
    dims[1]  = 0;
    filespace = H5Screate_simple(VECTOR_RANK, dims, NULL);

    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp1, number, H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    //create a local dataspace where each process has a few of the rows
    count[0]  = g->my_nnodes;
    count[1]  = dims[1];
    

    // give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(VECTOR_RANK, count, NULL);

    // determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    
    // the 2 is dimension of our data structure
    hsize_t   coord1[g->my_nnodes];
    
    // the values in coord1 should be the global node numbers, this time it is trivial
    int i,j,k=0,l=0;
    for (i=0; i<g->nnodes; i++){
        //temporarily store PE
        if(g->node[i].resident_pe == g->smpi->myid){
            scalardata[l] = g->node[i].resident_pe;
            coord1[l] = g->node[i].gid;
            l+=1;

    }
    }

    assert(l==g->my_nnodes);
    H5Sselect_elements(filespace, H5S_SELECT_SET,g->my_nnodes,(const hsize_t *)&coord1);
    
    // Create property list for collective dataset write.
    // declare collextive data file writing

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #ifdef _MESSG
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    #endif
    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, scalardata);
    free(scalardata);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp1);
    
    // +++++++++++++++++++++++++++++++++++++++
    // Nodal write complete
    // +++++++++++++++++++++++++++++++++++++++

    //close mesh file
    H5Fclose(file_id);
    
#endif

}




/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Print an AdH grid to XDMF
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] g           (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void sgrid_write_elemental_pe(SGRID *g){



#ifdef _ADH_HDF5

    hid_t     file_id, plist_id, filespace, dset_id, memspace;
    hid_t     grp1,grp2;
    char      fname[50];
    hsize_t   dims[MATRIX_RANK];
    herr_t    status;
    hsize_t   count[MATRIX_RANK],offset[MATRIX_RANK];
    int nt=0;


    // for adaptive grid writing
    char number[50] = "";
    sprintf(number, "%d", nt);

    // mesh properties
    int *connectivity;

    // create object to store h5 config
    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    #ifdef _MESSG
        MPI_Info info  = MPI_INFO_NULL;
        H5Pset_fapl_mpio(plist_id, g->smpi->ADH_COMM, info);
    #endif


    // open file
    strcpy(fname,g->filename);
    strcat(fname, ".h5");
    file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);

    // group 1 is anything with mesh nodes
    grp1 = H5Gopen(file_id, "/Data/ElementalScalar", H5P_DEFAULT); // a group is like a folder

    //++++++++++++++++++++++++++++++++++++++
    // Writing Nodes
    //++++++++++++++++++++++++++++++++++++++

    // declare global sizes of data set nodes first. This is local data.
    //Nodally based scalar data, maybe outsource this to other routine later

    int *scalardata; 
    scalardata = (int *)malloc(sizeof(int) * (g->my_nelems3d + g->my_nelems2d + g->my_nelems1d));




    // create a dataspace.  This is a global quantity
    dims[0]  = g->macro_nelems3d + g->macro_nelems2d + g->macro_nelems1d;
    dims[1]  = 0;
    filespace = H5Screate_simple(VECTOR_RANK, dims, NULL);

    //create a parallel dataset object and close filespace
    dset_id = H5Dcreate(grp1, number, H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    //create a local dataspace where each process has a few of the rows
    count[0]  = g->my_nelems3d + g->my_nelems2d + g->my_nelems1d;
    count[1]  = dims[1];
    

    // give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(VECTOR_RANK, count, NULL);

    // determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    
    // the 2 is dimension of our data structure
    hsize_t   coord1[g->my_nelems3d + g->my_nelems2d + g->my_nelems1d];
    
    // the values in coord1 should be the global node numbers, this time it is trivial
    int i,j,k=0,l=0;
    for (i=0; i<g->nelems3d; i++){
        //temporarily store PE
        if(g->elem3d[i].resident_pe == g->smpi->myid){
            scalardata[l] = g->elem3d[i].resident_pe;
            coord1[l] = g->elem3d[i].gid;
            l+=1;
    }
    }

    for (i=0; i<g->nelems2d; i++){
        //temporarily store PE
        if(g->elem2d[i].resident_pe == g->smpi->myid){
            scalardata[l] = g->elem2d[i].resident_pe;
            coord1[l] = g->elem2d[i].gid + g->macro_nelems3d;
            l+=1;
    }
    }

    //currently 1d not available to have PE plotted

    assert(l==g->my_nelems3d + g->my_nelems2d + g->my_nelems1d);
    H5Sselect_elements(filespace, H5S_SELECT_SET,  g->my_nelems3d + g->my_nelems2d + g->my_nelems1d  ,(const hsize_t *)&coord1);
    
    // Create property list for collective dataset write.
    // declare collextive data file writing

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #ifdef _MESSG
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    #endif
    //collective write
    status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, scalardata);
    free(scalardata);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(grp1);
    
    // +++++++++++++++++++++++++++++++++++++++
    // Nodal write complete
    // +++++++++++++++++++++++++++++++++++++++

    //close mesh file
    H5Fclose(file_id);
    
#endif

}



/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Read a nodal attribute
 *  \author    Corey Trahan, Ph.D.
 *  \author    Mark Loveand, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] g           (SGRID *)  pointer to an AdH grid
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void sgrid_read_nodal_attribute(SGRID *g){



    #ifdef _ADH_HDF5

    int i,j=0;
    int nsurf,my_nsurf;
    int lnode3d;

    //hard code one nodal attribute?
    g->nodal_attribute.rain = (double *) tl_alloc(sizeof(double), nsurf);

    
    //Note, need these to match to avoid above loop
    nsurf = g->nnodes_sur;
    //doesnt appear to be neccesary
    //my_nsurf = g->my_nnodes_sur;
    //printf("PE %d || nsurf %d || my nsurf %d || macro_nsurf%d\n",g->smpi->myid,g->nnodes_sur,g->my_nnodes_sur,g->macro_nnodes_sur);
    int local_to_global_sur[nsurf];
    int nodeID_2d_to_3d_sur[nsurf];
    //set up map to go from local to global
    for (i=0;i<g->nnodes;i++){
        //printf("PE %d || global id %d Global surface id %d\n",g->smpi->myid,g->node[i].gid, g->node[i].global_surf_id);
        if (g->node[i].bflag==SURFACE || g->node[i].bflag==BODY2D){
            //not necessary for read but will need this map most likely in residual routines
            //add to grid class?
            nodeID_2d_to_3d_sur[j] = i;
            //necessary for reading surface properly
            local_to_global_sur[j]=g->node[i].global_surf_id;
            j+=1;
        }
    }
    assert(j==nsurf);
    //key attribute here, 
    //g->node[local_index].global_surf_id
    //localSurf2GlobalMap[nnodes_surf]; initialize to UNSET_INT
    //count = 0
    //for (i=0; i<nnodes_sur; i++) {
     //lnode3d = nodeID_2d_to_3d_sur[i]
     //gnode3d = g->node[lnode3d].gid
    //localSurf2GlobalMap[count] = gnode3d;
     //count ++;
    //}


    //attempt to read in H5 file in parallel
    hid_t     file_id, plist_id, filespace, dset_id, memspace;
    char      fname[50];
    hsize_t   dims[MATRIX_RANK];
    herr_t    status;
    hsize_t   count[MATRIX_RANK],offset[MATRIX_RANK];
    int nt=0;
    // for adaptive grid writing
    char number[50] = "";
    sprintf(number, "%d", nt);

    // create object to store h5 config
    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    #ifdef _MESSG
        MPI_Info info  = MPI_INFO_NULL;
        H5Pset_fapl_mpio(plist_id, g->smpi->ADH_COMM, info);
    #endif


    // open file

    strcpy(fname,g->filename);
    strcat(fname, ".att.h5");
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);

    

    //++++++++++++++++++++++++++++++++++++++
    // Writing Nodes
    //++++++++++++++++++++++++++++++++++++++

    // declare global sizes of data set nodes first. This is local data.
    //Nodally based scalar data, maybe outsource this to other routine later

    //double *scalardata; 
    //scalardata = (double *)malloc(sizeof(double) * nsurf);



    //create a parallel dataset object and close filespace
    //For now hard code but this should be an argument into the function
    dset_id = H5Dopen(file_id, "/NodalAttributes/Rain", H5P_DEFAULT);
    

    //create a local dataspace where each process has a few of the rows
    count[0]  = nsurf;
    count[1]  = 0;
    

    // give the part that each processor will give
    // this one assumes each process is the same
    // there is way to specify each count per process
    memspace  = H5Screate_simple(VECTOR_RANK, count, NULL);

    // determine the hyperslabs in the file
    filespace = H5Dget_space(dset_id);
    
    // the 2 is dimension of our data structure
    hsize_t   coord1[nsurf];

    //fill in
    for (i=0; i<nsurf; i++){
        coord1[i] = local_to_global_sur[i];
        //printf("PE %d || local surface no. %d , global surface id %d\n",g->smpi->myid,i,mapping[i]);
    }


    H5Sselect_elements(filespace, H5S_SELECT_SET,nsurf,(const hsize_t *)&coord1);
    
    // Create property list for collective dataset write.
    // declare collextive data file writing

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #ifdef _MESSG
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    #endif
    //collective write
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, g->nodal_attribute.rain);



    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    
    // +++++++++++++++++++++++++++++++++++++++
    // Nodal write complete
    // +++++++++++++++++++++++++++++++++++++++

    //close mesh file
    H5Fclose(file_id);


    

    for (i=0; i<nsurf; i++){
            printf("PE %d, global sur id %d, gid %d, Rain data = %f\n",g->smpi->myid,local_to_global_sur[i],g->node[nodeID_2d_to_3d_sur[i]].gid,g->nodal_attribute.rain[i]);
    }



    #endif

}

