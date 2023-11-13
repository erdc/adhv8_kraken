/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads external boundary/face mesh data.
 *  \author    Charlie Berger, Ph.D.
 *  \author    Corey Trahan, Ph.D.
 *  \author    Gaurav Savant, Ph.D.
 *  \author    Gary Brown, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout]  grid                 (SGRID *) a pointer to the model grid structure
 * @param[in]  io                       (SIO *)  a pointer to input/output file structure
 * @param[in]  JUST_COUNT     (int)  a flag to indicate whether the file is only for counting faces
 * @param[in]  elem_nnodes   (int)  the number of nodes on an element face
 * @param[in]  data                  (int)  a line from the face file
 *
 * \note CJT \:: The order that data is read here is IMPORTANT!
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"
static int DEBUG = OFF;
static int DEBUG_CHECK_PICKETS = OFF;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_a_face(SIO *io, SGRID *grid, int JUST_COUNT, int elem_nnodes, char *data) {
    
    int i, j, k, ielem = UNSET_INT, bflag = UNSET_INT, istring = UNSET_INT, ielem3d_local = 0;
    int nd[MAX_NNODES_ON_ELEM2D]; sarray_init_value_int(nd,MAX_NNODES_ON_ELEM2D,UNSET_INT);
    
    // the global/macro 2d element number
    grid->macro_nelems2d++;
    
    if (JUST_COUNT==YES) {
        grid->nelems2d = grid->macro_nelems2d; // cjt :: temporary assignment // max allocated for now
        return;
    }
    
    // get 3d element the 2d face belongs to
    ielem = read_int_field(*io, &data);
    if (ielem > grid->macro_nelems3d) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        io_read_error(*io, "Bad element number.", TRUE);
    }
    ielem--;
    
    // read in macro/global node IDs
    for (i=0;i<elem_nnodes;i++) {
        nd[i] = read_int_field(*io, &data) - 1;
    }
    
    // read boundary flag
    bflag = read_int_field(*io, &data);
    if (bflag < 0 || bflag > 2) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        io_read_error(*io, "Bad boundary flag for face.", TRUE);
    }
    
    // read material string
    istring = read_int_field_custom(*io, &data, NULL, "boundary string ID", UNSET_INT, TRUE);
    istring--;
    
    // store local 2d element if the parent 3d element is local
    for (ielem3d_local=0;ielem3d_local<grid->nelems3d; ielem3d_local++) {
        if(grid->elem3d[ielem3d_local].gid==ielem){
            grid->elem2d[grid->nelems2d].nnodes = elem_nnodes;
            selem2d_alloc(&(grid->elem2d[grid->nelems2d]), elem_nnodes);
            // these need to be local node ids!
            k = 0;
            for (i=0; i<grid->elem3d[ielem3d_local].nnodes; i++) {
                for (j=0; j<grid->elem2d[grid->nelems2d].nnodes; j++) {
                    if (nd[j] == grid->node[grid->elem3d[ielem3d_local].nodes[i]].gid){
                        grid->elem2d[grid->nelems2d].nodes[j] = grid->elem3d[ielem3d_local].nodes[i];
                        k++;
                    }
                }
            }
            if (k != grid->elem2d[grid->nelems2d].nnodes) {
                printf("\nk: %d\n",k);
                tl_error("\nThere seems to be something wrong with the face file input.");
            }
#ifdef _DEBUG
           if (DEBUG){
                printf("\nelem3d[%i].nnodes = %i, elem2d[%i].nnodes = %i",
                        ielem3d_local, grid->elem3d[ielem3d_local].nnodes,
                        grid->nelems2d, grid->elem2d[grid->nelems2d].nnodes);
                printf("\n    elem3d.nodes :- [%i(%i), %i(%i), %i(%i), %i(%i)]",
                        grid->elem3d[ielem3d_local].nodes[0], grid->node[grid->elem3d[ielem3d_local].nodes[0]].gid,
                        grid->elem3d[ielem3d_local].nodes[1], grid->node[grid->elem3d[ielem3d_local].nodes[1]].gid,
                        grid->elem3d[ielem3d_local].nodes[2], grid->node[grid->elem3d[ielem3d_local].nodes[2]].gid,
                        grid->elem3d[ielem3d_local].nodes[3], grid->node[grid->elem3d[ielem3d_local].nodes[3]].gid);
                printf("\n    elem2d.nodes :- [%i, %i, %i]; need gid: [%i, %i, %i]",
                        grid->elem2d[grid->nelems2d].nodes[0],
                        grid->elem2d[grid->nelems2d].nodes[1],
                        grid->elem2d[grid->nelems2d].nodes[2],
                        nd[0], nd[1], nd[2]);
          }
            for (i=0; i<grid->elem2d[grid->nelems2d].nnodes; i++){
                if (grid->elem2d[grid->nelems2d].nodes[i] == UNSET_INT){
                    tl_error("\nThere seems to be something wrong with the face file input.");
                }
            }
#endif
            
            
            grid->elem2d[grid->nelems2d].mat = grid->elem3d[ielem3d_local].mat;
            grid->elem2d[grid->nelems2d].id_3d = ielem3d_local;
            grid->elem2d[grid->nelems2d].string = istring;
            grid->elem2d[grid->nelems2d].bflag = bflag;
            grid->elem2d[grid->nelems2d].gid = grid->macro_nelems2d-1;
            grid->elem2d[grid->nelems2d].id = grid->nelems2d;
            
            if (elem_nnodes == NDONTRI) {
                grid->elem2d[grid->nelems2d].nnodes_quad = 6;
                grid->elem2d[grid->nelems2d].nedges = 3;
                grid->elem2d[grid->nelems2d].edges = grid->nd_on_TriEdge;
                grid->haveTris = TRUE;
            } else if (elem_nnodes == NDONQUAD) {
                grid->elem2d[grid->nelems2d].nnodes_quad = 8;
                grid->elem2d[grid->nelems2d].nedges = 4;
                grid->elem2d[grid->nelems2d].edges = grid->nd_on_QuadEdge;
                grid->haveQuads = TRUE;
            } else {
                tl_error("ERROR: Face card unrecognized!");
            }
            grid->nelems2d++;
            return;
        }
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_faces(SIO *io, SGRID *grid, int JUST_COUNT) {
    
    char line[MAXLINE];           /* the input line */
    char *data = NULL;            /* the data */
    int i, ielem=UNSET_INT, ifc=UNSET_INT, istring=UNSET_INT, idum=UNSET_INT;
    int bflag=UNSET_INT, icount = 0;
    int ielem3d_local;
    int gid=0, nd1, nd2, nd3, nd4, nd1_local, nd2_local, nd3_local, nd4_local;
    SVECT side1, side2, cross;
    
#ifdef _DEBUG
    if (grid->smpi->myid == 0) {
        printf("\n---- Reading face file");
        if (JUST_COUNT==YES) {printf(" (for counting purposes only)");}
    }
#endif
//    int super = strlen(io->sup.filename);
//    open_input_file( &(io->face), "boundary face file", super);
//    if (io->face.fp == NULL) {
//        printf("WARNING :: NO FACE FILE IS FOUND!\n");
//        return;
//    }
    
    // +++++++++++++++++++++++++++++++++
    // +++++++++++++++++++++++++++++++++
    grid->nelems2d=0;
    grid->macro_nelems2d=0;
    // +++++++++++++++++++++++++++++++++
    // +++++++++++++++++++++++++++++++++
    rewind(io->face.fp);
    while (fgets(line, MAXLINE, io->face.fp) != NULL) {
        io_save_line(io, io->face.fp, io->face.filename, line);
        if (strip_comments(line) <= 1) {continue;}    /* ignore empty ('\n') or comment lines */
        switch (parse_card(line, &data)) {
            case CARD_TRI:
                read_a_face(io, grid, JUST_COUNT, NDONTRI, data);
                break;
                /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
                /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
            case CARD_QUAD:
                read_a_face(io, grid, JUST_COUNT, NDONQUAD, data);
                break;
                /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
                /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
            case CARD_FCS: // FOR BACKWARD'S COMPATIBILITY
                
                // the global/macro 2d element number
                grid->macro_nelems2d++;
                
                if (JUST_COUNT==YES) {
                    grid->nelems2d=grid->macro_nelems2d; // cjt :: temporary assignment
                    break;
                }
                
                /* the nodes on a face */
                int nd_on_fc[NFCPRELM][NDPRFC] = { {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1} };
                
                // get 3d element the 2d face belongs to
                ielem = read_int_field(*io, &data);
                if (ielem > grid->macro_nelems3d) {
                    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                    io_read_error(*io, "Bad element number.", TRUE);
                }
                ielem--;
                
                int ifc = read_int_field(*io, &data);
                if (ifc > NFCPRELM) {
                    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                    io_read_error(*io, "Bad face number.", TRUE);
                }
                ifc--;
                
                bflag = read_int_field(*io, &data);
                if (bflag < 0 || bflag > 2) {
                    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                    io_read_error(*io, "Bad boundary flag for face.", TRUE);
                }
                
                istring = read_int_field_custom(*io, &data, NULL, "boundary string ID", UNSET_INT, TRUE);
                istring--;
                
                
                for (ielem3d_local=0;ielem3d_local<grid->nelems3d; ielem3d_local++) {
                    if(grid->elem3d[ielem3d_local].gid==ielem){
                        grid->elem2d[grid->nelems2d].nnodes = 3;
                        selem2d_alloc(&(grid->elem2d[grid->nelems2d]), NDONTRI);
                        grid->elem2d[grid->nelems2d].nodes[0] = grid->elem3d[ielem3d_local].nodes[nd_on_fc[ifc][0]];
                        grid->elem2d[grid->nelems2d].nodes[1] = grid->elem3d[ielem3d_local].nodes[nd_on_fc[ifc][1]];
                        grid->elem2d[grid->nelems2d].nodes[2] = grid->elem3d[ielem3d_local].nodes[nd_on_fc[ifc][2]];
                        grid->elem2d[grid->nelems2d].mat = grid->elem3d[ielem3d_local].mat;
                        grid->elem2d[grid->nelems2d].id_3d = ielem3d_local;
                        grid->elem2d[grid->nelems2d].string = istring;
                        grid->elem2d[grid->nelems2d].bflag = bflag;
                        grid->elem2d[grid->nelems2d].gid = grid->macro_nelems2d-1;
                        grid->elem2d[grid->nelems2d].id = grid->nelems2d;
                        grid->elem2d[grid->nelems2d].nnodes_quad = 6;
                        grid->elem2d[grid->nelems2d].nedges = 3;
                        grid->elem2d[grid->nelems2d].edges = grid->nd_on_TriEdge;
                        grid->haveTris = TRUE;
                        grid->nelems2d++;
                        break;
                    }
                }
                break;
            default:
                break;
        }
    }
    io_save_line(io, NULL, "", "");
    
    /* sanity checks */
    if (grid->nelems2d <= 0) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> no 2d elements listed in face file\n");
    }
    
    //fclose(io->face.fp);
    rewind(io->face.fp);
    
#ifdef _DEBUG
    if (grid->smpi->myid == 0) {
        printf("\n---- Finished reading face file");
        if (JUST_COUNT == YES) printf(" (for counting purposes only)");
        printf("\n");
    }
    
    if (DEBUG_CHECK_PICKETS) tl_check_all_pickets(__FILE__, __LINE__);
#endif
}
