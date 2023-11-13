#include "global_header.h"

void sgrid_read_adjacency(SGRID *grid, char *file_base);
double sgrid_min(SGRID *, const char);
double sgrid_max(SGRID *, const char);

// ----------------------------------------------------------
// ----------------------------------------------------------
void sgrid_read_adh(SGRID **pgrid, char *file_base) {
    
    //RUN_VERBOSE = true;
    
    SGRID *grid = *pgrid;
    
    FILE *fptr = NULL;
    char filename[100], binary[100];
    char *line = NULL, card[100], dum[100];
    size_t len = 0;
    ssize_t readfile = -1;
    double x,y,z;
    int i,j,k, id, nd1, nd2, nd3, nd4, imat, ie1, ie2;
    
    // first check for binary
    strcpy(binary,file_base);
    strcat(binary, ".bin");
    fptr = fopen(binary,"rb");
    if (fptr == NULL) {
        // binary not found, reading geometry file
        strcpy(filename,file_base);
        strcat(filename, ".3dm");
        if (RUN_VERBOSE) printf("-- READING ADH GEO FILE:  %s \n",filename);
        if ((fptr = fopen(filename,"r+")) == NULL){
            printf("Error! opening geo file. If you are running -t, you must be in verification directory.\n");
            exit(-1);
        }
    } else {
        if (RUN_VERBOSE) printf("-- READING BINARY GEO FILE:  %s \n",binary);
        sgrid_read_binary(grid,fptr);
        if (RUN_VERBOSE) sgrid_printScreen(grid);
        return;
    }
    
    assert(fptr);
 
    // initialize vertical displacement flag
    grid->vertical_dpl_flag = OFF;

    // triangle edge layout
    grid->nd_on_TriEdge[0][0] = 0; grid->nd_on_TriEdge[0][1] = 1;
    grid->nd_on_TriEdge[1][0] = 1; grid->nd_on_TriEdge[1][1] = 2;
    grid->nd_on_TriEdge[2][0] = 2; grid->nd_on_TriEdge[2][1] = 0;
    
    // tetrahedron edge layout
    grid->nd_on_TetEdge[0][0] = 0; grid->nd_on_TetEdge[0][1] = 1;
    grid->nd_on_TetEdge[1][0] = 0; grid->nd_on_TetEdge[1][1] = 2;
    grid->nd_on_TetEdge[2][0] = 0; grid->nd_on_TetEdge[2][1] = 3;
    grid->nd_on_TetEdge[3][0] = 1; grid->nd_on_TetEdge[3][1] = 2;
    grid->nd_on_TetEdge[4][0] = 1; grid->nd_on_TetEdge[4][1] = 3;
    grid->nd_on_TetEdge[5][0] = 2; grid->nd_on_TetEdge[5][1] = 3;
    
    // count nodes and elements
    grid->nelems1d = 0;
    grid->nelems2d = 0;
    grid->nelems3d = 0;
    grid->nnodes = 0;
    grid->ndim = UNSET_INT;
    
    readfile = getline(&line, &len, fptr); // "mesh"
    //printf("%s", line);
    while ((readfile = getline(&line, &len, fptr)) != -1) {
        //printf("%s", line);
        sscanf(line,"%s %s",card,dum);
        //printf("%s %s \n",card,dum);
        if (strcmp(card, "TET") == 0) {grid->nelems3d++; grid->ndim = 3;}
        if (strcmp(card, "TRI") == 0) {grid->nelems2d++; grid->ndim = 2;}
        if (strcmp(card, "E3T") == 0) {grid->nelems2d++; grid->ndim = 2;}
        if (strcmp(card, "ND") == 0)  grid->nnodes++;
        
        line = NULL;
    }
    assert(grid->nnodes > 0);
    assert(grid->nelems2d > 0 || grid->nelems3d > 0);
    
    //printf("nnodes: %d nelems2d: %d nelems3d: %d \n",grid->nnodes,grid->nelems2d,grid->nelems3d);
    
    // allocated grid arrays
    if (grid->nelems2d > 0) grid->elem2d = (SELEM_2D *) calloc(grid->nelems2d, sizeof(SELEM_2D));
    if (grid->nelems3d > 0) grid->elem3d = (SELEM_3D *) calloc(grid->nelems3d, sizeof(SELEM_3D));
    //grid->node_sur =  (int *)  calloc(grid->nnodes, sizeof(int)); // not used atm
    grid->node     = (SVECT *) calloc(grid->nnodes, sizeof(SVECT));
    grid->node_t0  = (SVECT *) calloc(grid->nnodes, sizeof(SVECT));
    grid->node_flag = (int *)  calloc(grid->nnodes, sizeof(int));
    
    selem2d_init_array(grid->elem2d, grid->nelems2d);
    selem3d_init_array(grid->elem3d, grid->nelems3d);
    
    
    // read and store
    rewind(fptr);
    line = NULL; len = 0;
    readfile = getline(&line, &len, fptr); // "mesh"
    while ((readfile = getline(&line, &len, fptr)) != -1) {
        sscanf(line,"%s",card);
        if (strcmp(card, "TET") == 0) {
            sscanf(line,"%s %d %d %d %d %d",dum,&id,&nd1,&nd2,&nd3,&nd4);
            selem3d_alloc(&grid->elem3d[id-1], 4);
            grid->elem3d[id-1].nodes[0] = nd1-1;
            grid->elem3d[id-1].nodes[1] = nd2-1;
            grid->elem3d[id-1].nodes[2] = nd3-1;
            grid->elem3d[id-1].nodes[3] = nd4-1;
        }
        if ((strcmp(card, "TRI") == 0) || (strcmp(card, "E3T") == 0)) {
            sscanf(line,"%s %d %d %d %d %d",dum,&id,&nd1,&nd2,&nd3,&imat);
            selem2d_alloc(&grid->elem2d[id-1], 3);
            //printf("TRI id: %d \n",id-1);
            grid->elem2d[id-1].nodes[0] = nd1-1;
            grid->elem2d[id-1].nodes[1] = nd2-1;
            grid->elem2d[id-1].nodes[2] = nd3-1;
        }
        if (strcmp(card, "ND") == 0)  {
            sscanf(line,"%s %d %lf %lf %lf",dum,&id,&x,&y,&z);
            grid->node[id-1].x = x;
            grid->node[id-1].y = y;
            grid->node[id-1].z = z;
            grid->node_t0[id-1].x = x;
            grid->node_t0[id-1].y = y;
            grid->node_t0[id-1].z = z;
            //printf("NODE id: %d %f %f %f \n",id-1,grid->node[id-1].x,grid->node[id-1].y,grid->node[id-1].z);
        }
        
        line = NULL;
    }
    
    fclose(fptr);

    // compute grid min/max and extent
    grid->xmin = sgrid_min(grid, 'x');
    grid->ymin = sgrid_min(grid, 'y');

    grid->xmax = sgrid_max(grid, 'x');
    grid->ymax = sgrid_max(grid, 'y');

    grid->xL = grid->xmax - grid->xmin;
    grid->yL = grid->ymax - grid->ymin;

    grid->xLinv = 1.0 / grid->xL;
    grid->yLinv = 1.0 / grid->yL;

    if (RUN_VERBOSE) printf("---- FINISHED READING ADH GEO FILE \n");
    
        // ---------------------------------------------------
        // read face file if the grid is 3D
        // ---------------------------------------------------
    //    if (grid->ndim == 3) {
    //        FILE *fptr_face = NULL;
    //        char filename_face[100];
    //        strcpy(filename_face,file_base);
    //        strcat(filename_face, ".faces");
    //
    //
    //        printf("-- READING ADH GEO FACE FILE:  %s \n",filename_face);
    //
    //        if ((fptr_face = fopen(filename_face,"r+")) == NULL){
    //            printf("Error! opening geo file");
    //            exit(-1);
    //        }
    //        assert(fptr_face);
    //
    //        line = NULL; len = 0;
    //        readfile = getline(&line, &len, fptr_face); // "mesh"
    //        grid->nelems2d=0;
    //        while ((readfile = getline(&line, &len, fptr_face)) != -1) {
    //            sscanf(line,"%s",card);
    //            if (strcmp(card, "TRI") == 0) grid->nelems2d++;
    //            line = NULL;
    //        }
    //        printf("---- # of 3d element external faces: %d\n",grid->nelems2d);
    //        if (grid->nelems2d > 0) grid->elem2d = (SELEM_2D *) malloc(grid->nelems2d);
    //
    //
    //        rewind(fptr_face);
    //        int id_3d=UNSET_INT,ibflag=UNSET_INT;
    //        line = NULL; len = 0;
    //        readfile = getline(&line, &len, fptr_face); // "mesh"
    //        k=0;
    //        while ((readfile = getline(&line, &len, fptr_face)) != -1) {
    //            sscanf(line,"%s",card);
    //            if (strcmp(card, "TRI") == 0) {
    //                sscanf(line,"%s %d %d %d %d %d",dum,&id_3d,&nd1,&nd2,&nd3,&ibflag,imat);
    //                //printf("ielem2d: %d id_3d: %d \n",k,id_3d);
    //                selem2d_alloc(&grid->elem2d[k], 3);
    //                grid->elem2d[k].id_3d = id_3d-1;
    //                grid->elem2d[k].nodes[0] = nd1;
    //                grid->elem2d[k].nodes[1] = nd2;
    //                grid->elem2d[k].nodes[2] = nd3;
    //                grid->elem2d[k].bflag = ibflag;
    //                k++;
    //            }
    //            line = NULL;
    //        }
    //    }
    //    printf("---- FINISHED READING ADH GEO FACE FILE \n");
    
    // ---------------------------------------------------
    // creating node connection look up for all grids
    // cjt :: note :: this is pretty fast
    // ---------------------------------------------------
    if (RUN_VERBOSE) printf("-- CREATING NODE CONNECTION TABLE \n");
    grid->nc_nelems = (int *) calloc(grid->nnodes, sizeof(int));
    grid->nc_elems = (int **) calloc(grid->nnodes, sizeof(int *));
    
    if (grid->ndim == 2) {

        for (i=0; i<grid->nnodes; i++) {
            grid->nc_nelems[i] = 0;
            for (ie1=0; ie1<grid->nelems2d; ie1++) {
                for (j=0; j<3; j++) {
                    if (grid->elem2d[ie1].nodes[j] == i) {
                        grid->nc_nelems[i]++;
                        break;
                    }
                }
            }
            grid->nc_elems[i] = (int *) calloc(grid->nc_nelems[i], sizeof(int));
            k=0;
            for (ie1=0; ie1<grid->nelems2d; ie1++) {
                for (j=0; j<3; j++) {
                    if (grid->elem2d[ie1].nodes[j] == i) {
                        grid->nc_elems[i][k] = ie1;
                        k++;
                        break;
                        
                    }
                }
            }
        }
    } else if (grid->ndim == 3) {
        for (i=0; i<grid->nnodes; i++) {
            grid->nc_nelems[i] = 0;
            for (ie1=0; ie1<grid->nelems3d; ie1++) {
                for (j=0; j<4; j++) {
                    if (grid->elem3d[ie1].nodes[j] == i) {
                        grid->nc_nelems[i]++;
                        break;
                    }
                }
            }
            grid->nc_elems[i] = (int *) calloc(grid->nc_nelems[i], sizeof(int));
            k=0;
            for (ie1=0; ie1<grid->nelems3d; ie1++) {
                for (j=0; j<4; j++) {
                    if (grid->elem3d[ie1].nodes[j] == i) {
                        grid->nc_elems[i][k] = ie1;
                        k++;
                        break;
                        
                    }
                }
            }
        }
    }
    if (RUN_VERBOSE) printf("---- FINISHED CREATING NODE CONNECTION TABLE \n");
    
    // ---------------------------------------------------
    // creating face adjacency look up for 3d grids
    // cjt :: note :: slow as christmas
    // ---------------------------------------------------
    int elem3d_lface_flag[grid->nelems3d][4];
    int iface1, iface2, nd1_ID, nd2_ID, nd3_ID, nd1_2_ID, nd2_2_ID, nd3_2_ID, match;
    if (grid->ndim == 3) {
        if (RUN_VERBOSE) printf("-- CREATING FACE ADJACENCY TABLE \n");
        for (i=0; i<grid->nnodes; i++) {
            grid->node_flag[i] =  UNSET_INT;
            
            for (j=0; j<grid->nc_nelems[i]; j++) {
                ie1 = grid->nc_elems[i][j];
                for (iface1=0; iface1<4; iface1++) {
                    if (grid->elem3d[ie1].face_flag[iface1] != UNSET_INT) continue;
                    grid->elem3d[ie1].face_flag[iface1] = UNSET_INT;
                    elem3d_lface_flag[ie1][iface1] = UNSET_INT;
                    nd1_ID = grid->elem3d[ie1].nodes[nd_on_fc[iface1][0]];
                    nd2_ID = grid->elem3d[ie1].nodes[nd_on_fc[iface1][1]];
                    nd3_ID = grid->elem3d[ie1].nodes[nd_on_fc[iface1][2]];
                    
                    match = 0;
                    for (k=0; k<grid->nc_nelems[i]; k++) {
                        if (j==k) continue; // same element
                        ie2 = grid->nc_elems[i][k];
                        for (iface2=0; iface2<4; iface2++) {
                            nd1_2_ID = grid->elem3d[ie2].nodes[nd_on_fc[iface2][0]];
                            nd2_2_ID = grid->elem3d[ie2].nodes[nd_on_fc[iface2][1]];
                            nd3_2_ID = grid->elem3d[ie2].nodes[nd_on_fc[iface2][2]];
                            
                            if ((nd1_ID == nd1_2_ID && nd2_ID == nd2_2_ID && nd3_ID == nd3_2_ID) ||
                                (nd1_ID == nd1_2_ID && nd2_ID == nd3_2_ID && nd3_ID == nd2_2_ID) ||
                                (nd1_ID == nd2_2_ID && nd2_ID == nd1_2_ID && nd3_ID == nd3_2_ID) ||
                                (nd1_ID == nd2_2_ID && nd2_ID == nd3_2_ID && nd3_ID == nd1_2_ID) ||
                                (nd1_ID == nd3_2_ID && nd2_ID == nd2_2_ID && nd3_ID == nd1_2_ID) ||
                                (nd1_ID == nd3_2_ID && nd2_ID == nd1_2_ID && nd3_ID == nd2_2_ID)) {
                                
                                grid->elem3d[ie1].face_flag[iface1] = ie2; // found a match
                                elem3d_lface_flag[ie1][iface1] = iface2;
                                match = 1;
                            }
                            if (match == 1) break;
                        }
                        if (match == 1) break;
                    }
                }
            }
        }
        
        // use the above to flag nodes and edges as external boundaries
        int nd[3];
        for (ie1=0; ie1<grid->nelems3d; ie1++) {
            for (iface1=0; iface1<4; iface1++) {
                if (grid->elem3d[ie1].face_flag[iface1] == UNSET_INT) {
                    // face nodes
                    nd[0] = grid->elem3d[ie1].nodes[nd_on_fc[iface1][0]];
                    nd[1] = grid->elem3d[ie1].nodes[nd_on_fc[iface1][1]];
                    nd[2] = grid->elem3d[ie1].nodes[nd_on_fc[iface1][2]];
                    // nodes
                    for (j=0; j<3; j++) grid->node_flag[nd[j]] = 1;
                    // edges
                    for (j=0; j<6; j++) {
                        // edge nodes
                        nd1_ID = grid->elem3d[ie1].nodes[grid->nd_on_TetEdge[j][0]];
                        nd2_ID = grid->elem3d[ie1].nodes[grid->nd_on_TetEdge[j][1]];
                        // is edge on this face
                        grid->elem3d[ie1].edge_flag[j] = UNSET_INT;
                        if ((nd1_ID == nd[0] && nd2_ID == nd[1]) ||
                            (nd1_ID == nd[1] && nd2_ID == nd[0]) ||
                            (nd1_ID == nd[0] && nd2_ID == nd[2]) ||
                            (nd1_ID == nd[2] && nd2_ID == nd[0]) ||
                            (nd1_ID == nd[1] && nd2_ID == nd[2]) ||
                            (nd1_ID == nd[2] && nd2_ID == nd[1]) ) {
                            grid->elem3d[ie1].edge_flag[j] = 1;
                        }
                    }
                }
            }
        }
        
//        int nt = 230;
//        for (i=0; i<grid->nc_nelems[nt]; i++) {
//            printf("node[%d] {%f %f %f} :: elem_connect[%d]: %d \n",nt,grid->node[nt].x,grid->node[nt].y,grid->node[nt].z,
//                   i,grid->nc_elems[nt][i]);
//        }
        if (RUN_VERBOSE) printf("---- FINISHED CREATING FACE ADJACENCY TABLE \n");
    }
        
    // ---------------------------------------------------
    // creating edge adjacency look up for 2d grids
    // ---------------------------------------------------
    if (RUN_VERBOSE) printf("-- CREATING EDGE ADJACENCY TABLE \n");
    int iedge1, iedge2, nd1_1, nd1_2, nd2_1, nd2_2;
    if (grid->ndim == 2) {
        for (ie1=0; ie1<grid->nelems2d; ie1++) {
            for (iedge1=0; iedge1<3; iedge1++) {
                grid->elem2d[ie1].edge_flag[iedge1] = UNSET_INT;
                nd1_1 = grid->elem2d[ie1].nodes[ grid->nd_on_TriEdge[iedge1][0] ];
                nd2_1 = grid->elem2d[ie1].nodes[ grid->nd_on_TriEdge[iedge1][1] ];
                for (ie2=0; ie2<grid->nelems2d; ie2++) {
                    if (ie1==ie2) continue;
                    
                    for (iedge2=0; iedge2<3; iedge2++) {
                        nd1_2 = grid->elem2d[ie2].nodes[ grid->nd_on_TriEdge[iedge2][0] ];
                        nd2_2 = grid->elem2d[ie2].nodes[ grid->nd_on_TriEdge[iedge2][1] ];
                        
                        if ((nd1_1 == nd1_2 && nd2_1 == nd2_2) || (nd1_1 == nd2_2 && nd2_1 == nd1_2) ) {
                            grid->elem2d[ie1].edge_flag[iedge1] = ie2; // found a match
                            break;
                        }
                    }
                    if (grid->elem2d[ie1].edge_flag[iedge1] != UNSET_INT) break;
                }
            }
        }
        
        for (ie1=0; ie1<grid->nelems2d; ie1++) {
            for (iedge1=0; iedge1<3; iedge1++) {
                if (grid->elem2d[ie1].edge_flag[iedge1] == UNSET_INT) {
                    // edge is on boundary
                    nd1_1 = grid->elem2d[ie1].nodes[ grid->nd_on_TriEdge[iedge1][0] ];
                    nd2_1 = grid->elem2d[ie1].nodes[ grid->nd_on_TriEdge[iedge1][1] ];
                    grid->node_flag[nd1_1] = 1;
                    grid->node_flag[nd2_1] = 1;
                }
            }
        }
    } else if (grid->ndim == 3) {
        
        // first count
        for (ie1=0; ie1<grid->nelems3d; ie1++) {
            for (iedge1=0; iedge1<6; iedge1++) {
                grid->elem3d[ie1].nedge2elem[iedge1] = 0;
                nd1_1 = grid->elem3d[ie1].nodes[ grid->nd_on_TetEdge[iedge1][0] ];
                nd2_1 = grid->elem3d[ie1].nodes[ grid->nd_on_TetEdge[iedge1][1] ];
                for (j=0; j<grid->nc_nelems[nd1_1]; j++) {
                    ie2 = grid->nc_elems[nd1_1][j];
                    if (ie1 == ie2) continue;
                    for (iedge2=0; iedge2<6; iedge2++) {
                        nd1_2 = grid->elem3d[ie2].nodes[ grid->nd_on_TetEdge[iedge2][0] ];
                        nd2_2 = grid->elem3d[ie2].nodes[ grid->nd_on_TetEdge[iedge2][1] ];
                        if ((nd1_1 == nd1_2 && nd2_1 == nd2_2) || (nd1_1 == nd2_2 && nd2_1 == nd1_2) ) {
                            grid->elem3d[ie1].nedge2elem[iedge1]++; // found a match
                            break;
                        }
                    }
                }
            }
        }
        
        // allocate edge to element map
        for (ie1=0; ie1<grid->nelems3d; ie1++) {
            grid->elem3d[ie1].edge2elem  = (int **) calloc(6, sizeof(int *));
            for (iedge1=0; iedge1<6; iedge1++) {
                grid->elem3d[ie1].edge2elem[iedge1] = (int *) calloc(grid->elem3d[ie1].nedge2elem[iedge1], sizeof(int));
            }
        }

        // now find map
        int count = 0;
        for (ie1=0; ie1<grid->nelems3d; ie1++) {
            for (iedge1=0; iedge1<6; iedge1++) {
                count = 0;
                nd1_1 = grid->elem3d[ie1].nodes[ grid->nd_on_TetEdge[iedge1][0] ];
                nd2_1 = grid->elem3d[ie1].nodes[ grid->nd_on_TetEdge[iedge1][1] ];
                for (j=0; j<grid->nc_nelems[nd1_1]; j++) {
                    ie2 = grid->nc_elems[nd1_1][j];
                    if (ie1 == ie2) continue;
                    for (iedge2=0; iedge2<6; iedge2++) {
                        nd1_2 = grid->elem3d[ie2].nodes[ grid->nd_on_TetEdge[iedge2][0] ];
                        nd2_2 = grid->elem3d[ie2].nodes[ grid->nd_on_TetEdge[iedge2][1] ];
                        if ((nd1_1 == nd1_2 && nd2_1 == nd2_2) || (nd1_1 == nd2_2 && nd2_1 == nd1_2) ) {
                            grid->elem3d[ie1].edge2elem[iedge1][count] = ie2; // found a match
                            count++;
                            break;
                        }
                    }
                }
            }
        }

        // print
//        for (ie1=0; ie1<grid->nelems3d; ie1++) {
//            printf("elem3d[%d]\n",ie1);
//            for (iedge1=0; iedge1<6; iedge1++) {
//                nd1_1 = grid->elem3d[ie1].nodes[ grid->nd_on_TetEdge[iedge1][0] ];
//                nd2_1 = grid->elem3d[ie1].nodes[ grid->nd_on_TetEdge[iedge1][1] ];
//                printf("---- local edge: %d || node: %d %d || nedge2elem: %d || elements:",iedge1,nd1_1,nd2_1,grid->elem3d[ie1].nedge2elem[iedge1]);
//                for (k=0;k<grid->elem3d[ie1].nedge2elem[iedge1];k++) {
//                    printf(" %d ",grid->elem3d[ie1].edge2elem[iedge1][k]);
//                }
//                printf("\n");
//            }
//        }
        
    }
    if (RUN_VERBOSE) printf("---- FINISHED CREATING EDGE ADJACENCY TABLE \n");
    

    
    
    //sgrid_read_adjacency(grid,file_base);
    
    // write Adh basic grid info to binary
    fptr = fopen(binary,"wb");
    sgrid_write_binary(grid,fptr);
    fclose(fptr);
    
    if (RUN_VERBOSE) sgrid_printScreen(grid);
    
}
// ----------------------------------------------------------
// ----------------------------------------------------------
void sgrid_write_binary(SGRID *grid, FILE *fptr) {
    int i,j,k;
    
    //fwrite(grid, sizeof(SGRID), 1, fptr);
    fwrite(&grid->vertical_dpl_flag, sizeof(int), 1, fptr);
    fwrite(&grid->ndim, sizeof(int), 1, fptr);
    fwrite(&grid->isGridMixedElement, sizeof(int), 1, fptr);
    fwrite(&grid->haveTets, sizeof(int), 1, fptr);
    fwrite(&grid->haveTris, sizeof(int), 1, fptr);
    fwrite(&grid->haveQuads, sizeof(int), 1, fptr);
    fwrite(&grid->type, sizeof(int), 1, fptr);
    fwrite(&grid->nnodes, sizeof(int), 1, fptr);
    fwrite(&grid->nelems1d, sizeof(int), 1, fptr);
    fwrite(&grid->nelems2d, sizeof(int), 1, fptr);
    fwrite(&grid->nelems3d, sizeof(int), 1, fptr);
    fwrite(&grid->nedges, sizeof(int), 1, fptr);
    fwrite(&grid->nmat, sizeof(int), 1, fptr);
    for (i=0; i<grid->nelems3d; i++) {
        fwrite(&grid->elem3d[i].id, sizeof(int), 1, fptr);
        fwrite(&grid->elem3d[i].nnodes, sizeof(int), 1, fptr);
        fwrite(grid->elem3d[i].nodes, sizeof(int), grid->elem3d[i].nnodes, fptr);
        fwrite(&grid->elem3d[i].djac, sizeof(double), 1, fptr);
        fwrite(&grid->elem3d[i].mat, sizeof(int), 1, fptr);
        //fwrite(&grid->elem3d[i].nadjc_elems, sizeof(int), 1, fptr);
        //fwrite(&grid->elem3d[i].adjc_elems, sizeof(int), grid->elem3d[i].nadjc_elems, fptr);
        fwrite(grid->elem3d[i].face_flag, sizeof(int), 4, fptr);
        fwrite(grid->elem3d[i].elem2d_face, sizeof(int), 4, fptr);
        fwrite(grid->elem3d[i].edge_flag, sizeof(int), 6, fptr);
        fwrite(grid->elem3d[i].nedge2elem, sizeof(int), 6, fptr);
        for (j=0; j<6; j++) {
            for (k=0; k<grid->elem3d[i].nedge2elem[j]; k++) {
                fwrite(&grid->elem3d[i].edge2elem[j][k], sizeof(int), 1, fptr);
            }
        }
        //exit(-1);
    }
    for (i=0; i<grid->nelems2d; i++) {
        fwrite(&grid->elem2d[i].id, sizeof(int), 1, fptr);
        fwrite(&grid->elem2d[i].id_3d, sizeof(int), 1, fptr);
        fwrite(&grid->elem2d[i].nnodes, sizeof(int), 1, fptr);
        fwrite(grid->elem2d[i].nodes, sizeof(int), grid->elem2d[i].nnodes, fptr);
        fwrite(&grid->elem2d[i].djac, sizeof(double), 1, fptr);
        fwrite(&grid->elem2d[i].djac3d, sizeof(double), 1, fptr);
        fwrite(&grid->elem2d[i].nrml.x, sizeof(double), 1, fptr);
        fwrite(&grid->elem2d[i].nrml.y, sizeof(double), 1, fptr);
        fwrite(&grid->elem2d[i].nrml.z, sizeof(double), 1, fptr);
        fwrite(&grid->elem2d[i].mat, sizeof(int), 1, fptr);
        fwrite(&grid->elem2d[i].bflag, sizeof(int), 1, fptr);
        //fwrite(&grid->elem2d[i].nadjc_elems, sizeof(int), 1, fptr);
        //fwrite(&grid->elem2d[i].adjc_elems, sizeof(int), grid->elem2d[i].nadjc_elems, fptr);
        fwrite(grid->elem2d[i].edge_flag, sizeof(int), 3, fptr);
    }
    for (i=0; i<grid->nnodes; i++) {
        fwrite(&grid->node[i].x, sizeof(double), 1, fptr);
        fwrite(&grid->node[i].y, sizeof(double), 1, fptr);
        fwrite(&grid->node[i].z, sizeof(double), 1, fptr);
        fwrite(&grid->node_t0[i].x, sizeof(double), 1, fptr);
        fwrite(&grid->node_t0[i].y, sizeof(double), 1, fptr);
        fwrite(&grid->node_t0[i].z, sizeof(double), 1, fptr);
        fwrite(&grid->node_flag[i], sizeof(int), 1, fptr);
        fwrite(&grid->nc_nelems[i], sizeof(int), 1, fptr);
        for (j=0; j<grid->nc_nelems[i]; j++) fwrite(&grid->nc_elems[i][j], sizeof(int), 1, fptr);
    }
    fwrite(&grid->xmin, sizeof(double), 1, fptr);
    fwrite(&grid->xmax, sizeof(double), 1, fptr);
    fwrite(&grid->ymin, sizeof(double), 1, fptr);
    fwrite(&grid->ymax, sizeof(double), 1, fptr);
    fwrite(&grid->xL, sizeof(double), 1, fptr);
    fwrite(&grid->yL, sizeof(double), 1, fptr);
    fwrite(&grid->xLinv, sizeof(double), 1, fptr);
    fwrite(&grid->yLinv, sizeof(double), 1, fptr);
    fwrite(&grid->mesh_volume, sizeof(double), 1, fptr);
    for (i=0; i<3; i++) {
        for (j=0; j<2; j++) fwrite(&grid->nd_on_TriEdge[i][j], sizeof(int), 1, fptr);
    }
    for (i=0; i<6; i++) {
        for (j=0; j<2; j++) fwrite(&grid->nd_on_TetEdge[i][j], sizeof(int), 1, fptr);
    }
}

// ----------------------------------------------------------
// ----------------------------------------------------------
void sgrid_read_binary(SGRID *grid, FILE *fptr) {
    int i,j,k;
    
    fread(&grid->vertical_dpl_flag, sizeof(int), 1, fptr);
    fread(&grid->ndim, sizeof(int), 1, fptr);
    fread(&grid->isGridMixedElement, sizeof(int), 1, fptr);
    fread(&grid->haveTets, sizeof(int), 1, fptr);
    fread(&grid->haveTris, sizeof(int), 1, fptr);
    fread(&grid->haveQuads, sizeof(int), 1, fptr);
    fread(&grid->type, sizeof(int), 1, fptr);
    fread(&grid->nnodes, sizeof(int), 1, fptr);
    fread(&grid->nelems1d, sizeof(int), 1, fptr);
    fread(&grid->nelems2d, sizeof(int), 1, fptr);
    fread(&grid->nelems3d, sizeof(int), 1, fptr);
    fread(&grid->nedges, sizeof(int), 1, fptr);
    fread(&grid->nmat, sizeof(int), 1, fptr);

    
    if (grid->nelems2d > 0) {
        grid->elem2d = (SELEM_2D *) calloc(grid->nelems2d, sizeof(SELEM_2D));
        selem2d_init_array(grid->elem2d, grid->nelems2d);
    }
    if (grid->nelems3d > 0) {
        grid->elem3d = (SELEM_3D *) calloc(grid->nelems3d, sizeof(SELEM_3D));
        selem3d_init_array(grid->elem3d, grid->nelems3d);
    }
    grid->node =    (SVECT *) calloc(grid->nnodes, sizeof(SVECT));
    grid->node_t0 = (SVECT *) calloc(grid->nnodes, sizeof(SVECT));
    grid->node_flag = (int *)  calloc(grid->nnodes, sizeof(int));
    grid->nc_nelems = (int *)  calloc(grid->nnodes, sizeof(int));
    grid->nc_elems  = (int **) calloc(grid->nnodes, sizeof(int *));
    
    for (i=0; i<grid->nelems3d; i++) {
        fread(&grid->elem3d[i].id, sizeof(int), 1, fptr);
        selem3d_alloc(&grid->elem3d[i],4);
        fread(&grid->elem3d[i].nnodes, sizeof(int), 1, fptr);
        fread(grid->elem3d[i].nodes, sizeof(int), grid->elem3d[i].nnodes, fptr);
        fread(&grid->elem3d[i].djac, sizeof(double), 1, fptr);
        fread(&grid->elem3d[i].mat, sizeof(int), 1, fptr);
//        fread(&grid->elem3d[i].nadjc_elems, sizeof(int), 1, fptr);
//        fread(&grid->elem3d[i].adjc_elems, sizeof(int), grid->elem3d[i].nadjc_elems, fptr);
        fread(grid->elem3d[i].face_flag, sizeof(int), 4, fptr);
        fread(grid->elem3d[i].elem2d_face, sizeof(int), 4, fptr);
        fread(grid->elem3d[i].edge_flag, sizeof(int), 6, fptr);
        fread(grid->elem3d[i].nedge2elem, sizeof(int), 6, fptr);
        grid->elem3d[i].edge2elem  = (int **) calloc(6, sizeof(int *));
        for (j=0; j<6; j++) {
            grid->elem3d[i].edge2elem[j] = (int *) calloc(grid->elem3d[i].nedge2elem[j], sizeof(int));
            for (k=0; k<grid->elem3d[i].nedge2elem[j]; k++) {
                fread(&grid->elem3d[i].edge2elem[j][k], sizeof(int), 1, fptr);
            }
        }
    }
    for (i=0; i<grid->nelems2d; i++) {
        fread(&grid->elem2d[i].id, sizeof(int), 1, fptr);
        selem2d_alloc(&grid->elem2d[i],3);
        fread(&grid->elem2d[i].id_3d, sizeof(int), 1, fptr);
        fread(&grid->elem2d[i].nnodes, sizeof(int), 1, fptr);
        fread(grid->elem2d[i].nodes, sizeof(int), grid->elem2d[i].nnodes, fptr);
        fread(&grid->elem2d[i].djac, sizeof(double), 1, fptr);
        fread(&grid->elem2d[i].djac3d, sizeof(double), 1, fptr);
        fread(&grid->elem2d[i].nrml.x, sizeof(double), 1, fptr);
        fread(&grid->elem2d[i].nrml.y, sizeof(double), 1, fptr);
        fread(&grid->elem2d[i].nrml.z, sizeof(double), 1, fptr);
        fread(&grid->elem2d[i].mat, sizeof(int), 1, fptr);
        fread(&grid->elem2d[i].bflag, sizeof(int), 1, fptr);
//        fread(&grid->elem2d[i].nadjc_elems, sizeof(int), 1, fptr);
//        fread(&grid->elem2d[i].adjc_elems, sizeof(int), grid->elem2d[i].nadjc_elems, fptr);
        fread(grid->elem2d[i].edge_flag, sizeof(int), 3, fptr);
    }

    for (i=0; i<grid->nnodes; i++) {
        fread(&grid->node[i].x, sizeof(double), 1, fptr);
        fread(&grid->node[i].y, sizeof(double), 1, fptr);
        fread(&grid->node[i].z, sizeof(double), 1, fptr);
        fread(&grid->node_t0[i].x, sizeof(double), 1, fptr);
        fread(&grid->node_t0[i].y, sizeof(double), 1, fptr);
        fread(&grid->node_t0[i].z, sizeof(double), 1, fptr);
        fread(&grid->node_flag[i], sizeof(int), 1, fptr);
        fread(&grid->nc_nelems[i], sizeof(int), 1, fptr);
        grid->nc_elems[i] = (int *) calloc(grid->nc_nelems[i], sizeof(int));
        for (j=0; j<grid->nc_nelems[i]; j++) fread(&grid->nc_elems[i][j], sizeof(int), 1, fptr);
    }

    fread(&grid->xmin, sizeof(double), 1, fptr);
    fread(&grid->xmax, sizeof(double), 1, fptr);
    fread(&grid->ymin, sizeof(double), 1, fptr);
    fread(&grid->ymax, sizeof(double), 1, fptr);
    fread(&grid->xL, sizeof(double), 1, fptr);
    fread(&grid->yL, sizeof(double), 1, fptr);
    fread(&grid->xLinv, sizeof(double), 1, fptr);
    fread(&grid->yLinv, sizeof(double), 1, fptr);
    fread(&grid->mesh_volume, sizeof(double), 1, fptr);
    for (i=0; i<3; i++) {
        for (j=0; j<2; j++) fread(&grid->nd_on_TriEdge[i][j], sizeof(int), 1, fptr);
    }
    for (i=0; i<6; i++) {
        for (j=0; j<2; j++) fread(&grid->nd_on_TetEdge[i][j], sizeof(int), 1, fptr);
    }
}


// ----------------------------------------------------------
// ----------------------------------------------------------

void sgrid_read_adjacency(SGRID *grid, char *file_base) {
    
    FILE *fptr;
    char filename[100];
    char *line=NULL, card[100], dum[100];
    size_t len = 0;
    ssize_t readfile = -1;
    double x,y,z;
    int i, id, nadjc_elems,n,ielem, found;
    
    strcpy(filename,file_base);
    strcat(filename, ".adj");
    
    fptr = fopen(filename,"r");
    
    if(fptr != NULL) {
        // file exists
        printf("-- READING ADH GEO ADJACENCY FILE:  %s \n",filename);
        
        readfile = getline(&line, &len, fptr); // "mesh"
        //printf("%s", line);
        while ((readfile = getline(&line, &len, fptr)) != -1) {
            //printf("%s", line);
            sscanf(line,"%s %d %d %n",card,&id,&nadjc_elems,&n);
            //printf("card: %s id: %d nadcj: %d\n",card,id,nadjc_elems);
            line += n;
            
            if (strcmp(card, "TET") == 0) {
                assert(nadjc_elems >= 4 && nadjc_elems < 30);
                grid->elem3d[id-1].nadjc_elems = nadjc_elems;
                grid->elem3d[id-1].adjc_elems = (int *) calloc(nadjc_elems, sizeof(int));
                
                for (i=0; i<nadjc_elems; i++) {
                    sscanf(line,"%d %n",&ielem,&n);
                    grid->elem3d[id-1].adjc_elems[i] = ielem-1;
                    line += n;
                }
                
            }
            if (strcmp(card, "TRI") == 0) {
                assert(nadjc_elems >= 3 && nadjc_elems < 30);
                grid->elem2d[id-1].nadjc_elems = nadjc_elems;
                grid->elem2d[id-1].adjc_elems = (int *) calloc(nadjc_elems, sizeof(int));
                
                for (i=0; i<nadjc_elems; i++) {
                    //printf("i: %d\n",i);
                    sscanf(line,"%d %n",&ielem,&n);
                    //printf("TRI line: %s\n",line);
                    grid->elem2d[id-1].adjc_elems[i] = ielem-1;
                    line += n;
                }
            }
            line = NULL;
        }
        printf("-- FINISHED READING ADH GEO ADJACENCY FILE:  %s \n",filename);
        
    } else {
        fclose(fptr);
        fptr = fopen(filename,"w");
        
        // file doesn't exist :: create adjacency info
        //  find 2d adjacent elements
        int ie, ie2, i, j, match, count;
        
        if (grid->nelems2d > 0) {
            
            for (ie=0; ie<grid->nelems2d; ie++) {
                grid->elem2d[ie].nadjc_elems = 0;
                for (ie2=0; ie2<grid->nelems2d; ie2++) {
                    if (ie == ie2) continue;
                    found = 0;
                    for (i=0; i<3; i++) {
                        for (j=0; j<3; j++) {
                            if (grid->elem2d[ie].nodes[i] == grid->elem2d[ie2].nodes[j]) {
                                grid->elem2d[ie].nadjc_elems++;
                                found = 1;
                                break;
                            }
                        }
                        if (found == 1) break;
                    }
                }
            }
            
            for (ie=0; ie<grid->nelems2d; ie++) {
                printf("ie %d adj_elems: %d\n",ie,grid->elem2d[ie].nadjc_elems);
                grid->elem2d[ie].adjc_elems = (int *) calloc(grid->elem2d[ie].nadjc_elems, sizeof(int));
            }
            
            for (ie=0; ie<grid->nelems2d; ie++) {
                count = 0;
                for (ie2=0; ie2<grid->nelems2d; ie2++) {
                    if (ie == ie2) continue;
                    found = 0;
                    for (i=0; i<3; i++) {
                        for (j=0; j<3; j++) {
                            if (grid->elem2d[ie].nodes[i] == grid->elem2d[ie2].nodes[j]) {
                                grid->elem2d[ie].adjc_elems[count] = ie2;
                                found = 1;
                                count++;
                                break;
                            }
                        }
                        if (found == 1) break;
                    }
                }
            }
        }
        
        if (grid->nelems3d > 0) {
            
            for (ie=0; ie<grid->nelems3d; ie++) {
                grid->elem3d[ie].adjc_elems = (int *) calloc(grid->elem3d[ie].nadjc_elems, sizeof(int));
            }
            
            for (ie=0; ie<grid->nelems3d; ie++) {
                count = 0;
                for (ie2=0; ie2<grid->nelems3d; ie2++) {
                    if (ie == ie2) continue;
                    match = 0;
                    for (i=0; i<4; i++) {
                        for (j=0; j<4; j++) {
                            if (grid->elem3d[ie].nodes[i] == grid->elem3d[ie2].nodes[j]) match++;
                            break;
                        }
                        if (match == 3) {
                            grid->elem3d[ie].adjc_elems[count] = ie2;
                            count++;
                            break;
                        }
                    }
                }
            }
        }
        
        // print to file for future use
        fprintf(fptr,"MESH%dD\n",grid->ndim);
        if (grid->nelems2d > 0) {
            for (ie=0; ie<grid->nelems2d; ie++) {
                //printf("ie: %d\n",ie);
                fprintf(fptr,"TRI %d %d ",ie+1,grid->elem2d[ie].nadjc_elems);
                for (i=0; i<grid->elem2d[ie].nadjc_elems; i++) {
                    fprintf(fptr,"%d ",grid->elem2d[ie].adjc_elems[i]+1);
                }
                fprintf(fptr,"\n");
            }
        }
        if (grid->nelems3d > 0) {
            for (ie=0; ie<grid->nelems3d; ie++) {
                fprintf(fptr,"TET %d %d ",ie+1,grid->elem3d[ie].nadjc_elems);
                for (i=0; i<grid->elem3d[ie].nadjc_elems; i++) {
                    fprintf(fptr,"%d ",grid->elem3d[ie].adjc_elems[i]+1);
                }
                fprintf(fptr,"\n");
            }
        }
    }
    
    fclose(fptr);
}

// ----------------------------------------------------------
// ----------------------------------------------------------

void sgrid_free(SGRID *grid) {
    int inode, ie, i, k;
    
    if (grid != NULL) {
        
        //if (grid->elem1d != NULL && grid->nelems1d > 0) {
        //    selem1d_free_array(grid->elem1d, grid->nelems1d);
        //}
        if (grid->elem2d != NULL && grid->nelems2d > 0) {
            selem2d_free_array(grid->elem2d, grid->nelems2d);
        }
        if (grid->elem3d != NULL && grid->nelems3d > 0) {
            selem3d_free_array(grid->elem3d, grid->nelems3d);
        }
        
        if (grid->node != NULL &&grid->nnodes > 0) {
            free(grid->node);
        }
        
        free(grid);
    }
}

// ----------------------------------------------------------
// ----------------------------------------------------------

void sgrid_printScreen(SGRID *grid) {
    
    FILE *fp = stdout;
    
    //fprintf(fp,"\n");
    fprintf(fp,"***********************************************************\n");
    fprintf(fp,"***********************************************************\n");
    fprintf(fp,"---------- nnodes: %d \n", grid->nnodes);
    fprintf(fp,"---------- nelems1d: %d \n", grid->nelems1d);
    fprintf(fp,"---------- nelems2d: %d \n", grid->nelems2d);
    fprintf(fp,"---------- nelems3d: %d \n", grid->nelems3d);
    fprintf(fp,"---------- total number of materials: %d\n", grid->nmat);
    fprintf(fp,"---------- total mesh area: %14.6e\n", grid->mesh_volume);
    fprintf(fp,"***********************************************************\n");
    fprintf(fp,"***********************************************************\n");
    
}

// ----------------------------------------------------------
// ----------------------------------------------------------

void sgrid_create2d(SGRID *grid, double xmin, double xmax, double ymin, double ymax, int npx, int npy, double theta, double dz, double za, double zb, double zc) {
    
    if (RUN_VERBOSE) printf("-- CREATING ADH 2D GRID\n");
    
    // grid verticle displacement
    grid->vertical_dpl_flag = OFF;

    // triangle edge layout
    grid->nd_on_TriEdge[0][0] = 0; grid->nd_on_TriEdge[0][1] = 1;
    grid->nd_on_TriEdge[1][0] = 1; grid->nd_on_TriEdge[1][1] = 2;
    grid->nd_on_TriEdge[2][0] = 2; grid->nd_on_TriEdge[2][1] = 0;
    
    // count nodes and elements
    grid->nelems1d = 0;
    grid->nelems2d = 0;
    grid->nelems3d = 0;
    grid->nnodes = 0;
    grid->ndim = 2;
    
    SVECT ref_vector;
    ref_vector.x = 0 - 0;
    ref_vector.y = 0 - 0;
    ref_vector.z = -1 - 0;
    
    assert(npx > 1);
    assert(npy > 1);
    assert(dz > 0);
    if (npx > 2) assert(npx % 2 != 0);
    if (npy > 2) assert(npy % 2 != 0);
    
    assert(theta > -1e-6);
    theta *= 3.141592653589793 / 180.;
    
    double dx = (xmax - xmin)/(double)(npx-1);
    double dy = (ymax - ymin)/(double)(npy-1);
    
    grid->nnodes = npx * npy;
    grid->node = (SVECT *) calloc(grid->nnodes, sizeof(SVECT));
    grid->node_flag = (int *) calloc(grid->nnodes, sizeof(int));
    
    int i, j, k, ie1, ie2;
    double x = xmin, xr, zr;
    double y = ymin, yr;
    int npoints = 0;
    for (i=0; i<npx; i++) {
        for (j=0; j<npy; j++) {
            x = xmin + i*dx;
            y = ymax - j*dy; //ymin + j*dy;
            
            // totate 45 degrees
            xr = x*cos(theta) - y*sin(theta);
            yr = x*sin(theta) + y*cos(theta);
            zr = -(za + zb * x + zc * x * x);
            //printf("ND %d %20.10f %20.10f %20.10f\n",npoints+1,xr,yr,zr);
            
            grid->node[npoints].x = xr;
            grid->node[npoints].y = yr;
            grid->node[npoints].z = -(za + zb * x + zc * x * x);
            
            npoints++;
        }
    }
    
    SVECT cross;
    SVECT side1, side2;
    
    int igrp = 0, inode=0, nd1, nd2, nd3;
    grid->nelems2d=0;
    for (i=0; i<npx-1; i++) {
        for (igrp=0; igrp<npy-1; igrp++) {
            grid->nelems2d++;
            grid->nelems2d++;
        }
    }
    assert(grid->nelems2d > 0);
    grid->elem2d = (SELEM_2D *) calloc(grid->nelems2d, sizeof(SELEM_2D));
    
    k=0;
    for (i=0; i<npx-1; i++) {
        for (igrp=0; igrp<npy-1; igrp++) {
            inode = i*npy + igrp;
            
            // make sure node numbering is counter clock-wise
            nd1 = inode; nd2 = inode+1; nd3 = inode+npy+1;
            side1.x =  grid->node[nd2].x -  grid->node[nd1].x;
            side1.y =  grid->node[nd2].y -  grid->node[nd1].y;
            side1.z =  grid->node[nd2].z -  grid->node[nd1].z;
            side2.x =  grid->node[nd3].x -  grid->node[nd1].x;
            side2.y =  grid->node[nd3].y -  grid->node[nd1].y;
            side2.z =  grid->node[nd3].z -  grid->node[nd1].z;
            cross = svect_cross(side1,side2);
            
            if (svect_dotp(cross,ref_vector) > 0) {
                //if (cross.z > 0) {
                nd1 = inode; nd2 = inode+npy+1; nd3 = inode+1;
                side1.x =  grid->node[nd2].x -  grid->node[nd1].x;
                side1.y =  grid->node[nd2].y -  grid->node[nd1].y;
                side1.z =  grid->node[nd2].z -  grid->node[nd1].z;
                side2.x =  grid->node[nd3].x -  grid->node[nd1].x;
                side2.y =  grid->node[nd3].y -  grid->node[nd1].y;
                side2.z =  grid->node[nd3].z -  grid->node[nd1].z;
                cross = svect_cross(side1,side2);
                assert(svect_dotp(cross,ref_vector) < 0);
                //assert(cross.z < 0);
            }
            //printf("TRI %d nodes: %d %d %d \t cross: %20.10f %20.10f %20.10f \n",k+1,nd1+1,nd2+1,nd3+1,cross.x,cross.y,cross.z);
            
            selem2d_alloc(&grid->elem2d[k], 3);
            grid->elem2d[k].id = k;
            grid->elem2d[k].nodes[0] = nd1;
            grid->elem2d[k].nodes[1] = nd2;
            grid->elem2d[k].nodes[2] = nd3;
            grid->elem2d[k].mat = 0;
            k++;
            
            
            // make sure node numbering is counter clock-wise
            nd1 = inode; nd2 = inode+npy+1; nd3 = inode+npy;
            side1.x =  grid->node[nd2].x -  grid->node[nd1].x;
            side1.y =  grid->node[nd2].y -  grid->node[nd1].y;
            side1.z =  grid->node[nd2].z -  grid->node[nd1].z;
            side2.x =  grid->node[nd3].x -  grid->node[nd1].x;
            side2.y =  grid->node[nd3].y -  grid->node[nd1].y;
            side2.z =  grid->node[nd3].z -  grid->node[nd1].z;
            cross = svect_cross(side1,side2);
            if (svect_dotp(cross,ref_vector) > 0) {
                //if (cross.z > 0) {
                nd1 = inode; nd2 = inode+npy; nd3 = inode+npy+1;
                side1.x =  grid->node[nd2].x -  grid->node[nd1].x;
                side1.y =  grid->node[nd2].y -  grid->node[nd1].y;
                side1.z =  grid->node[nd2].z -  grid->node[nd1].z;
                side2.x =  grid->node[nd3].x -  grid->node[nd1].x;
                side2.y =  grid->node[nd3].y -  grid->node[nd1].y;
                side2.z =  grid->node[nd3].z -  grid->node[nd1].z;
                cross = svect_cross(side1,side2);
                //assert(cross.z < 0);
                assert(svect_dotp(cross,ref_vector) < 0);
            }
            
            //printf("TRI %d nodes: %d %d %d \t cross: %20.10f %20.10f %20.10f \n",k+1,nd1+1,nd2+1,nd3+1,cross.x,cross.y,cross.z);
            
            selem2d_alloc(&grid->elem2d[k], 3);
            grid->elem2d[k].id = k;
            grid->elem2d[k].nodes[0] = nd1;
            grid->elem2d[k].nodes[1] = nd2;
            grid->elem2d[k].nodes[2] = nd3;
            grid->elem2d[k].mat = 0;
            k++;
        }
    }
    if (RUN_VERBOSE) printf("---- FINISHED CREATING ADH 2D GRID\n");
    
    
    // ---------------------------------------------------
    // creating edge adjacency look up for 2d grids
    // ---------------------------------------------------
    if (RUN_VERBOSE) printf("-- CREATING EDGE ADJACENCY TABLE \n");
    int iedge1, iedge2, nd1_1, nd1_2, nd2_1, nd2_2;
    if (grid->ndim == 2) {
        for (ie1=0; ie1<grid->nelems2d; ie1++) {
            for (iedge1=0; iedge1<3; iedge1++) {
                grid->elem2d[ie1].edge_flag[iedge1] = UNSET_INT;
                nd1_1 = grid->elem2d[ie1].nodes[ grid->nd_on_TriEdge[iedge1][0] ];
                nd2_1 = grid->elem2d[ie1].nodes[ grid->nd_on_TriEdge[iedge1][1] ];
                for (ie2=0; ie2<grid->nelems2d; ie2++) {
                    if (ie1==ie2) continue;
                    
                    for (iedge2=0; iedge2<3; iedge2++) {
                        nd1_2 = grid->elem2d[ie2].nodes[ grid->nd_on_TriEdge[iedge2][0] ];
                        nd2_2 = grid->elem2d[ie2].nodes[ grid->nd_on_TriEdge[iedge2][1] ];
                        
                        if ((nd1_1 == nd1_2 && nd2_1 == nd2_2) || (nd1_1 == nd2_2 && nd2_1 == nd1_2) ) {
                            grid->elem2d[ie1].edge_flag[iedge1] = ie2; // found a match
                            break;
                        }
                    }
                }
            }
        }
        
        for (ie1=0; ie1<grid->nelems2d; ie1++) {
            for (iedge1=0; iedge1<3; iedge1++) {
                if (grid->elem2d[ie1].edge_flag[iedge1] == UNSET_INT) {
                    // edge is on boundary
                    nd1_1 = grid->elem2d[ie1].nodes[ grid->nd_on_TriEdge[iedge1][0] ];
                    nd2_1 = grid->elem2d[ie1].nodes[ grid->nd_on_TriEdge[iedge1][1] ];
                    grid->node_flag[nd1_1] = 1;
                    grid->node_flag[nd2_1] = 1;
                }
            }
        }
    }
    if (RUN_VERBOSE) printf("---- FINISHED CREATING EDGE ADJACENCY TABLE \n");
    
    // ---------------------------------------------------
    // creating node connection look up for all grids
    // ---------------------------------------------------
    if (RUN_VERBOSE) printf("-- CREATING NODE CONNECTION TABLE \n");
    grid->nc_nelems = (int *) calloc(grid->nelems2d, sizeof(int));
    grid->nc_elems = (int **) calloc(grid->nelems2d, sizeof(int *));
    for (i=0; i<grid->nnodes; i++) {
        grid->nc_nelems[i] = 0;
        for (ie1=0; ie1<grid->nelems2d; ie1++) {
            for (j=0; j<3; j++) {
                if (grid->elem2d[ie1].nodes[j] == i) {
                    grid->nc_nelems[i]++;
                    break;
                }
            }
        }
        grid->nc_elems[i] = (int *) calloc(grid->nc_nelems[i], sizeof(int));
        k=0;
        for (ie1=0; ie1<grid->nelems2d; ie1++) {
            for (j=0; j<3; j++) {
                if (grid->elem2d[ie1].nodes[j] == i) {
                    grid->nc_elems[i][k] = ie1;
                    k++;
                    break;
                    
                }
            }
        }
        //        printf("global node: %d || nelems: %d || elems: ",i,grid->nc_nelems[i]);
        //        for (j=0; j<grid->nc_nelems[i]; j++) {
        //            printf("%d ",grid->nc_elems[i][j]);
        //        }
        //        printf("\n");
        //        exit(-1);
    }
    if (RUN_VERBOSE) printf("---- FINISHED CREATING NODE CONNECTION TABLE \n");
}

// ----------------------------------------------------------
// ----------------------------------------------------------

void sgrid_create3d(SGRID *grid, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int npx, int npy, int npz, double theta) {
    
    printf("-- CREATING ADH 3D GRID\n");
    
    // grid verticle displacement
    grid->vertical_dpl_flag = OFF;

    // triangle edge layout
    grid->nd_on_TriEdge[0][0] = 0; grid->nd_on_TriEdge[0][1] = 1;
    grid->nd_on_TriEdge[1][0] = 1; grid->nd_on_TriEdge[1][1] = 2;
    grid->nd_on_TriEdge[2][0] = 2; grid->nd_on_TriEdge[2][1] = 0;
    
    // tetrahedron edge layout
    grid->nd_on_TetEdge[0][0] = 0; grid->nd_on_TetEdge[0][1] = 1;
    grid->nd_on_TetEdge[1][0] = 0; grid->nd_on_TetEdge[1][1] = 2;
    grid->nd_on_TetEdge[2][0] = 0; grid->nd_on_TetEdge[2][1] = 3;
    grid->nd_on_TetEdge[3][0] = 1; grid->nd_on_TetEdge[3][1] = 2;
    grid->nd_on_TetEdge[4][0] = 1; grid->nd_on_TetEdge[4][1] = 3;
    grid->nd_on_TetEdge[5][0] = 2; grid->nd_on_TetEdge[5][1] = 3;
    
    // count nodes and elements
    grid->nelems1d = 0;
    grid->nelems2d = 0;
    grid->nelems3d = 0;
    grid->nnodes = 0;
    grid->ndim = 2;
    
    SVECT ref_vector;
    ref_vector.x = 0 - 0;
    ref_vector.y = 0 - 0;
    ref_vector.z = -1 - 0;
    
    assert(npx > 1);
    assert(npy > 1);
    assert(npz > 1);
    if (npx > 2) assert(npx % 2 != 0);
    if (npy > 2) assert(npy % 2 != 0);
    if (npz > 2) assert(npz % 2 != 0);
    
    assert(theta > -1e-6);
    theta *= 3.141592653589793 / 180.;
    
    double dx = (xmax - xmin)/(double)(npx-1);
    double dy = (ymax - ymin)/(double)(npy-1);
    double dz = (zmax - zmin)/(double)(npz-1);
    
    grid->nnodes = npx * npy * npz;
    grid->node = (SVECT *) calloc(grid->nnodes, sizeof(SVECT));
    grid->node_flag = (int *) calloc(grid->nnodes, sizeof(int));
    
    int i, j, k, ie1, ie2, npoints = 0;
    double x = xmin, xr, zr;
    double y = ymin, yr;
    double z = xmin;
    for (i=0; i<npx; i++) {
        for (j=0; j<npy; j++) {
            for (k=0; k<npz; k++) {
                x = xmin + i*dx;
                y = ymax - j*dy;
                z = zmax - k*dz;
            
                // totate 45 degrees
                xr = x*cos(theta) - y*sin(theta);
                yr = x*sin(theta) + y*cos(theta);
            
                grid->node[npoints].x = xr;
                grid->node[npoints].y = yr;
                grid->node[npoints].z = z;
                
                printf("ND %d %20.10f %20.10f %20.10f\n",npoints+1,xr,yr,z);
            
                npoints++;
            }
        }
    }
    
//    for (i=0; i<npx; i++) {
//        for (j=0; j<npy; j++) {
//            for (k=0; k<npz; k++) {
//                nd0 = ip + (npz-1) ;
//
//    int nd0,nd1,nd2,nd3,nd4,nd5,nd6,nd7;
//    int ncubes = (npx-1) * (npy-1) * (npz - 1);
//    k=0;
//    for (i=0; i<ncubes; i++) {
//        nd0 =
//
//
//        grid->elem3d[k].nodes[0] = nd0;
//        grid->elem3d[k].nodes[1] = nd4;
//        grid->elem3d[k].nodes[2] = nd7;
//        grid->elem3d[k].nodes[3] = nd6;
//
//        k++;
//        grid->elem3d[k].nodes[0] = nd0;
//        grid->elem3d[k].nodes[1] = nd3;
//        grid->elem3d[k].nodes[2] = nd7;
//        grid->elem3d[k].nodes[3] = nd6;
//
//        k++
//        grid->elem3d[k].nodes[0] = nd0;
//        grid->elem3d[k].nodes[1] = nd4;
//        grid->elem3d[k].nodes[2] = nd5;
//        grid->elem3d[k].nodes[3] = nd6;
//
//        k++
//        grid->elem3d[k].nodes[0] = nd0;
//        grid->elem3d[k].nodes[1] = nd1;
//        grid->elem3d[k].nodes[2] = nd5;
//        grid->elem3d[k].nodes[3] = nd6;
//
//        k++;
//        grid->elem3d[k].nodes[0] = nd0;
//        grid->elem3d[k].nodes[1] = nd5;
//        grid->elem3d[k].nodes[2] = nd2;
//        grid->elem3d[k].nodes[3] = nd6;
//
//        k++;
//        grid->elem3d[k].nodes[0] = nd0;
//        grid->elem3d[k].nodes[1] = nd1;
//        grid->elem3d[k].nodes[2] = nd2;
//        grid->elem3d[k].nodes[3] = nd6;
//
//        k++;
//    }
    
//    SVECT cross;
//    SVECT side1, side2;
//
//    int igrp = 0, inode=0, nd1, nd2, nd3;
//    grid->nelems2d=0;
//    for (i=0; i<npx-1; i++) {
//        for (igrp=0; igrp<npy-1; igrp++) {
//            grid->nelems2d++;
//            grid->nelems2d++;
//        }
//    }
//    assert(grid->nelems2d > 0);
//    grid->elem2d = (SELEM_2D *) calloc(grid->nelems2d, sizeof(SELEM_2D));
//
//    k=0;
//    for (i=0; i<npx-1; i++) {
//        for (igrp=0; igrp<npy-1; igrp++) {
//            inode = i*npy + igrp;
//
//            // make sure node numbering is counter clock-wise
//            nd1 = inode; nd2 = inode+1; nd3 = inode+npy+1;
//            side1.x =  grid->node[nd2].x -  grid->node[nd1].x;
//            side1.y =  grid->node[nd2].y -  grid->node[nd1].y;
//            side1.z =  grid->node[nd2].z -  grid->node[nd1].z;
//            side2.x =  grid->node[nd3].x -  grid->node[nd1].x;
//            side2.y =  grid->node[nd3].y -  grid->node[nd1].y;
//            side2.z =  grid->node[nd3].z -  grid->node[nd1].z;
//            cross = svect_cross(side1,side2);
//
//            if (svect_dotp(cross,ref_vector) > 0) {
//                //if (cross.z > 0) {
//                nd1 = inode; nd2 = inode+npy+1; nd3 = inode+1;
//                side1.x =  grid->node[nd2].x -  grid->node[nd1].x;
//                side1.y =  grid->node[nd2].y -  grid->node[nd1].y;
//                side1.z =  grid->node[nd2].z -  grid->node[nd1].z;
//                side2.x =  grid->node[nd3].x -  grid->node[nd1].x;
//                side2.y =  grid->node[nd3].y -  grid->node[nd1].y;
//                side2.z =  grid->node[nd3].z -  grid->node[nd1].z;
//                cross = svect_cross(side1,side2);
//                assert(svect_dotp(cross,ref_vector) < 0);
//                //assert(cross.z < 0);
//            }
//            //printf("TRI %d nodes: %d %d %d \t cross: %20.10f %20.10f %20.10f \n",k+1,nd1+1,nd2+1,nd3+1,cross.x,cross.y,cross.z);
//
//            selem2d_alloc(&grid->elem2d[k], 3);
//            grid->elem2d[k].id = k;
//            grid->elem2d[k].nodes[0] = nd1;
//            grid->elem2d[k].nodes[1] = nd2;
//            grid->elem2d[k].nodes[2] = nd3;
//            grid->elem2d[k].mat = 0;
//            k++;
//
//
//            // make sure node numbering is counter clock-wise
//            nd1 = inode; nd2 = inode+npy+1; nd3 = inode+npy;
//            side1.x =  grid->node[nd2].x -  grid->node[nd1].x;
//            side1.y =  grid->node[nd2].y -  grid->node[nd1].y;
//            side1.z =  grid->node[nd2].z -  grid->node[nd1].z;
//            side2.x =  grid->node[nd3].x -  grid->node[nd1].x;
//            side2.y =  grid->node[nd3].y -  grid->node[nd1].y;
//            side2.z =  grid->node[nd3].z -  grid->node[nd1].z;
//            cross = svect_cross(side1,side2);
//            if (svect_dotp(cross,ref_vector) > 0) {
//                //if (cross.z > 0) {
//                nd1 = inode; nd2 = inode+npy; nd3 = inode+npy+1;
//                side1.x =  grid->node[nd2].x -  grid->node[nd1].x;
//                side1.y =  grid->node[nd2].y -  grid->node[nd1].y;
//                side1.z =  grid->node[nd2].z -  grid->node[nd1].z;
//                side2.x =  grid->node[nd3].x -  grid->node[nd1].x;
//                side2.y =  grid->node[nd3].y -  grid->node[nd1].y;
//                side2.z =  grid->node[nd3].z -  grid->node[nd1].z;
//                cross = svect_cross(side1,side2);
//                //assert(cross.z < 0);
//                assert(svect_dotp(cross,ref_vector) < 0);
//            }
//
//            //printf("TRI %d nodes: %d %d %d \t cross: %20.10f %20.10f %20.10f \n",k+1,nd1+1,nd2+1,nd3+1,cross.x,cross.y,cross.z);
//
//            selem2d_alloc(&grid->elem2d[k], 3);
//            grid->elem2d[k].id = k;
//            grid->elem2d[k].nodes[0] = nd1;
//            grid->elem2d[k].nodes[1] = nd2;
//            grid->elem2d[k].nodes[2] = nd3;
//            grid->elem2d[k].mat = 0;
//            k++;
//        }
//    }
//    printf("---- FINISHED CREATING ADH 2D GRID\n");
//
//
//    // ---------------------------------------------------
//    // creating edge adjacency look up for 2d grids
//    // ---------------------------------------------------
//    printf("-- CREATING EDGE ADJACENCY TABLE \n");
//    int iedge1, iedge2, nd1_1, nd1_2, nd2_1, nd2_2;
//    if (grid->ndim == 2) {
//        for (ie1=0; ie1<grid->nelems2d; ie1++) {
//            for (iedge1=0; iedge1<3; iedge1++) {
//                grid->elem2d[ie1].edge_flag[iedge1] = UNSET_INT;
//                nd1_1 = grid->elem2d[ie1].nodes[ grid->nd_on_TriEdge[iedge1][0] ];
//                nd2_1 = grid->elem2d[ie1].nodes[ grid->nd_on_TriEdge[iedge1][1] ];
//                for (ie2=0; ie2<grid->nelems2d; ie2++) {
//                    if (ie1==ie2) continue;
//
//                    for (iedge2=0; iedge2<3; iedge2++) {
//                        nd1_2 = grid->elem2d[ie2].nodes[ grid->nd_on_TriEdge[iedge2][0] ];
//                        nd2_2 = grid->elem2d[ie2].nodes[ grid->nd_on_TriEdge[iedge2][1] ];
//
//                        if ((nd1_1 == nd1_2 && nd2_1 == nd2_2) || (nd1_1 == nd2_2 && nd2_1 == nd1_2) ) {
//                            grid->elem2d[ie1].edge_flag[iedge1] = ie2; // found a match
//                            break;
//                        }
//                    }
//                }
//            }
//        }
//
//        for (ie1=0; ie1<grid->nelems2d; ie1++) {
//            for (iedge1=0; iedge1<3; iedge1++) {
//                if (grid->elem2d[ie1].edge_flag[iedge1] == UNSET_INT) {
//                    // edge is on boundary
//                    nd1_1 = grid->elem2d[ie1].nodes[ grid->nd_on_TriEdge[iedge1][0] ];
//                    nd2_1 = grid->elem2d[ie1].nodes[ grid->nd_on_TriEdge[iedge1][1] ];
//                    grid->node_flag[nd1_1] = 1;
//                    grid->node_flag[nd2_1] = 1;
//                }
//            }
//        }
//    }
//    printf("---- FINISHED CREATING EDGE ADJACENCY TABLE \n");
//
//    // ---------------------------------------------------
//    // creating node connection look up for all grids
//    // ---------------------------------------------------
//    printf("-- CREATING NODE CONNECTION TABLE \n");
//    grid->nc_nelems = (int *) calloc(grid->nelems2d, sizeof(int));
//    grid->nc_elems = (int **) calloc(grid->nelems2d, sizeof(int *));
//    for (i=0; i<grid->nnodes; i++) {
//        grid->nc_nelems[i] = 0;
//        for (ie1=0; ie1<grid->nelems2d; ie1++) {
//            for (j=0; j<3; j++) {
//                if (grid->elem2d[ie1].nodes[j] == i) {
//                    grid->nc_nelems[i]++;
//                    break;
//                }
//            }
//        }
//        grid->nc_elems[i] = (int *) calloc(grid->nc_nelems[i], sizeof(int));
//        k=0;
//        for (ie1=0; ie1<grid->nelems2d; ie1++) {
//            for (j=0; j<3; j++) {
//                if (grid->elem2d[ie1].nodes[j] == i) {
//                    grid->nc_elems[i][k] = ie1;
//                    k++;
//                    break;
//
//                }
//            }
//        }
//        //        printf("global node: %d || nelems: %d || elems: ",i,grid->nc_nelems[i]);
//        //        for (j=0; j<grid->nc_nelems[i]; j++) {
//        //            printf("%d ",grid->nc_elems[i][j]);
//        //        }
//        //        printf("\n");
//        //        exit(-1);
//    }
//    printf("---- FINISHED CREATING NODE CONNECTION TABLE \n");
}

// updates vertical displacement of nodes
void sgrid_dpl_update(SGRID *grid, double *dpl, double tL_dpl, double tR_dpl, double *dpl_tL, double *dpl_tR, double t) {

    int i;
    double max_dpl = -1e6;
    double max_z = -1e6;
    for (i=0; i<grid->nnodes; i++) {
        dpl[i] = interpolate1D(t,tL_dpl,tR_dpl,dpl_tL[i],dpl_tR[i]);
        //if (dpl[i] > max_dpl) {
            max_z = grid->node_t0[2420].z;
            max_dpl = dpl[2420];
        //}
        grid->node[i].z = grid->node_t0[i].z + dpl[i];
    }
    printf("sgrid_dpl_update: time: %f  tL_dpl: %f :: tR_dpl %f :: max_z: %f :: max_dpl: %f\n",t,tL_dpl,tR_dpl,max_z,max_dpl);
    
}

#define max(a,b) \
	({ __typeof__ (a) _a = (a); \
	 __typeof__ (b) _b = (b); \
	 _a > _b ? _a : _b; })
#define min(a,b) \
	({ __typeof__ (a) _a = (a); \
	 __typeof__ (b) _b = (b); \
	 _a < _b ? _a : _b; })

double sgrid_min(SGRID * grid, const char xy) {
    assert((xy == 'x') || (xy == 'y'));
    double themin = (xy == 'x') ? grid->node[0].x : grid->node[0].y;
    if (xy == 'x') {
		for (int i=1; i<grid->nnodes; i++)
			themin = min(themin, grid->node[i].x);
    } else {
		// Identical to the above loop except "x" replaced with "y"
		for (int i=1; i<grid->nnodes; i++)
			themin = min(themin, grid->node[i].y);
    }
    return themin;
}

// identical to above except "min" replaced with "max"
double sgrid_max(SGRID * grid, const char xy) {
    assert((xy == 'x') || (xy == 'y'));
    double themax = (xy == 'x') ? grid->node[0].x : grid->node[0].y;
    if (xy == 'x') {
		for (int i=1; i<grid->nnodes; i++)
			themax = max(themax, grid->node[i].x);
    } else {
		// Identical to the above loop except "x" replaced with "y"
		for (int i=1; i<grid->nnodes; i++)
			themax = max(themax, grid->node[i].y);
    }
    return themax;
}

#undef max
#undef min

void svect_normalize_inplace(bool normalize, SGRID * grid, SVECT * v) {
	if (normalize) {
		v->x = (v->x - grid->xmin) * grid->xLinv;
		v->y = (v->y - grid->ymin) * grid->yLinv;
	}
}

void svect_denormalize_inplace(bool normalize, SGRID * grid, SVECT * v) {
	if (normalize) {
		v->x = fma(grid->xL, v->x, grid->xmin);
		v->y = fma(grid->yL, v->y, grid->ymin);
	}
}

SVECT svect_normalize(bool normalize, SGRID * grid, SVECT in) {
	if (normalize) {
		SVECT v;
		v.x = (in.x - grid->xmin) * grid->xLinv;
		v.y = (in.y - grid->ymin) * grid->yLinv;
		v.z = in.z;
		return v;
	} else return in;
}

SVECT svect_denormalize(bool normalize, SGRID * grid, SVECT in) {
	if (normalize) {
		SVECT v;
		v.x = fma(grid->xL, in.x, grid->xmin);
		v.y = fma(grid->yL, in.y, grid->ymin);
		v.z = in.z;
		return v;
	} else return in;
}

void sgrid_normalize(bool normalize, SGRID * g) {

	if (normalize == false)
		return;

	for (int i=0; i<g->nnodes; i++)
		svect_normalize_inplace(true, g, &g->node[i]);

    return;
}
