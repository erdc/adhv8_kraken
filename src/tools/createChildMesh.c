#include "global_header.h"

static int DEBUG = OFF;

void createChildMesh(SGRID *grid, SIO *io) {
    
    int i, inode, nodeID_3d, ipe;
    char str[MAXLINE];
    int myid = grid->smpi->myid;
    int npes = grid->smpi->npes;
    
#ifdef _MESSG
    MPI_Barrier(grid->smpi->ADH_COMM);
#endif
    
    assert(grid->ndim == 3);
    
    // write grid header
    FILE *fp_2dm = io_fopen(build_filename(str, MAXLINE, io->proj_name, ".2dm"), "w", TRUE);
    if (grid->smpi->myid == 0) {
        fprintf(fp_2dm, "MESH2D\n");
        fprintf(fp_2dm, "Child mesh of 3d grid: %s\n", io->geo2d.filename);
    }
    fclose(fp_2dm);

    // get global surface grid node counts
    int gnnodes = 0;
#ifdef _MESSG
    MPI_Allreduce(&grid->my_nnodes, &gnnodes, 1, MPI_INT, MPI_SUM, grid->smpi->ADH_COMM);
#else
    gnnodes = grid->nnodes;
#endif
    
    
#ifdef _MESSG
    MPI_Barrier(grid->smpi->ADH_COMM); // all nodes must get here at once
#endif
    
#ifdef _DEBUG
    fflush(stdout);
    if (DEBUG) {
        printf("myid: %d \t my_nnodes: %d\n",myid, grid->my_nnodes);
        printf("CHILD GRID :: global node count: %d\n",gnnodes);
        fflush(stdout);
    }
#endif
    
    // Store all my_nnodes for each PE on all PES
    int *my_nnode_PE = (int *)tl_alloc(sizeof(int), npes);
#ifdef _MESSG
    MPI_Allgather(&grid->my_nnodes, 1, MPI_INT, my_nnode_PE, 1, MPI_INT, grid->smpi->ADH_COMM);
#ifdef _DEBUG
    if (DEBUG) {
        for (i = 0; i < npes; i++) {
            printf("myid: %d my_nnodes[PE=%d]: %d\n",myid, i, my_nnode_PE[i]);
        }
        messg_barrier(grid->smpi->ADH_COMM);
    }
#endif
#else
    my_nnode_PE[0] = grid->my_nnodes;
#endif
    
    // Now find the global map starting points for each set of my_nnodes
    int *start_index = (int *)tl_alloc(sizeof(int), npes);
    start_index[0] = 0;
    for (i = 1; i < npes; i++) {
        start_index[i] = start_index[i-1] + my_nnode_PE[i-1];
    }
    
    // allocate global to local node map
    int *node_map = (int *) tl_alloc(sizeof(int), gnnodes);
    sarray_init_value_int(node_map, gnnodes, UNSET_INT);
    
    // create node and node maps
    int kk=0, gindex = UNSET_INT, ierr = 0;
    for (ipe=0; ipe<npes; ipe++) {
        if (myid == ipe) {
            fp_2dm = io_fopen(build_filename(str, MAXLINE, io->proj_name, ".2dm"), "a", TRUE);
            for (i=0; i<grid->nnodes_sur; i++) {
                nodeID_3d = grid->nodeID_2d_to_3d_sur[i];
                if (grid->node[nodeID_3d].resident_pe == myid) {
                    gindex = start_index[ipe] + nodeID_3d;
#ifdef _DEBUG
                    if (DEBUG) {
                        printf("ND  %d %14.8e %14.8e %14.8e \n", kk+1, grid->node[nodeID_3d].x, grid->node[nodeID_3d].y, grid->node[nodeID_3d].z);
                    }
#endif
                    fprintf(fp_2dm, "ND  %d %14.8e %14.8e %14.8e \n", kk+1, grid->node[nodeID_3d].x, grid->node[nodeID_3d].y, grid->node[nodeID_3d].z);
                    node_map[gindex] = kk;
                    kk++;
                }
            }
#ifdef _DEBUG
            if (DEBUG) {
                printf("myid: %d \t kk: %d \t start_index[ipe]: %d\n",myid, kk,start_index[ipe]);
            }
#endif
            fclose(fp_2dm);
        }
#ifdef _MESSG
        ierr = MPI_Bcast(node_map, gnnodes, MPI_INT, ipe, grid->smpi->ADH_COMM);
        ierr = MPI_Bcast(&kk, 1, MPI_INT, ipe,  grid->smpi->ADH_COMM);
        MPI_Barrier(grid->smpi->ADH_COMM);
#endif
    }
    int ND = kk;
    
    // now write element and conectivities
    int gid_0, gid_1, gid_2;
    int ie, ie2d, pe0, pe1, pe2, nd0_3d, nd1_3d, nd2_3d, count=0;
    kk = 0;
    for (ipe=0; ipe<npes; ipe++) {
        if (myid == ipe) {
            fp_2dm = io_fopen(build_filename(str, MAXLINE, io->proj_name, ".2dm"), "a", TRUE);
            for (ie = 0; ie < grid->nelems2d_sur; ie++) {
                ie2d = grid->elem2d_sur[ie]; // 2d element id in the complete list of 2d elements
                nd0_3d = grid->elem2d[ie2d].nodes[0];
                nd1_3d = grid->elem2d[ie2d].nodes[1];
                nd2_3d = grid->elem2d[ie2d].nodes[2];
                
                // now find the PE each one of these is on - this will help find positions in global map
                pe0 = grid->node[nd0_3d].resident_pe;
                pe1 = grid->node[nd1_3d].resident_pe;
                pe2 = grid->node[nd2_3d].resident_pe;
                
                count = 0;
                if (pe0 == ipe) count++;
                if (pe1 == ipe) count++;
                if (pe2 == ipe) count++;
                
                if (count > 1) {
                    kk++;
                    gid_0 = start_index[grid->node[nd0_3d].resident_pe] + grid->node[nd0_3d].resident_id;
                    gid_1 = start_index[grid->node[nd1_3d].resident_pe] + grid->node[nd1_3d].resident_id;
                    gid_2 = start_index[grid->node[nd2_3d].resident_pe] + grid->node[nd2_3d].resident_id;
                    if (node_map[gid_0] == UNSET_INT) tl_error("something wrong with node0 map inverse");
                    if (node_map[gid_1] == UNSET_INT) tl_error("something wrong with node1 map inverse");
                    if (node_map[gid_2] == UNSET_INT) tl_error("something wrong with node2 map inverse");
                    fprintf(fp_2dm, "E3T  %d %6d %6d %6d %d\n", kk, node_map[gid_0]+1, node_map[gid_1]+1, node_map[gid_2]+1, grid->elem2d[ie].mat + 1);
                }
                
            }
            fclose(fp_2dm);
        }
#ifdef _MESSG
        ierr = MPI_Bcast(&kk, 1, MPI_INT, ipe,  grid->smpi->ADH_COMM);
        MPI_Barrier(grid->smpi->ADH_COMM);
#endif
    }
    int NC = kk;
    
    // write 2D hotstart
    FILE *fp_2d_hot = io_fopen(build_filename(str, MAXLINE, io->proj_name, "_2d.hot"), "w", TRUE);
    if (grid->smpi->myid == 0) {
        fprintf(fp_2d_hot, "DATASET\nOBJTYPE \"mesh2d\"\nBEGSCL\nND %d\nNC %d\nNAME \"ioh\"\nTIMEUNITS seconds\n",ND,NC);
    }
    fclose(fp_2d_hot);
    kk=0; gindex = UNSET_INT; ierr = 0;
    for (ipe=0; ipe<npes; ipe++) {
        if (myid == ipe) {
            fp_2d_hot = io_fopen(build_filename(str, MAXLINE, io->proj_name, "_2d.hot"), "a", TRUE);
            for (i=0; i<grid->nnodes_bed; i++) {
                nodeID_3d = grid->nodeID_2d_to_3d_bed[i];
                if (grid->node[nodeID_3d].resident_pe == myid) {
                    gindex = start_index[ipe] + nodeID_3d;
#ifdef _DEBUG
                    if (DEBUG) {
                        printf("ND  %d %14.8e \n", kk+1, -grid->node[nodeID_3d].z);
                    }
#endif
                    fprintf(fp_2dm, "%14.8e \n", -grid->node[nodeID_3d].z);
                    node_map[gindex] = kk;
                    kk++;
                }
            }
#ifdef _DEBUG
            if (DEBUG) {
                printf("myid: %d \t kk: %d \t start_index[ipe]: %d\n",myid, kk,start_index[ipe]);
            }
#endif
            fclose(fp_2d_hot);
        }
#ifdef _MESSG
        ierr = MPI_Bcast(node_map, gnnodes, MPI_INT, ipe, grid->smpi->ADH_COMM);
        ierr = MPI_Bcast(&kk, 1, MPI_INT, ipe,  grid->smpi->ADH_COMM);
        MPI_Barrier(grid->smpi->ADH_COMM);
#endif
    }
    
    // free arrays
    start_index = (int *) tl_free(sizeof(int), npes, start_index);
    my_nnode_PE = (int *) tl_free(sizeof(int), npes, my_nnode_PE);
    node_map = (int *) tl_free(sizeof(int), gnnodes, node_map);
}

