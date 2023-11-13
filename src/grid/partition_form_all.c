/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  partition_main.c This file partitions the grid for HPC usage.   */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"
#include "parmetis.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Paritions a model grid domain for distributed computing.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout]  g        (SGRID*) a pointer to an AdH model grid to be partitioned
 * @param[in]  flag           (int) a flag to indicate = 1 no division of processors necessary or = 0 processors get divided in superfile read
 *
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static int DEBUG = 0;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// initialize edge hash table
void edge_hashtab_init(EDGE_LIST_ITEM **edge) {
    int i;  for (i = 0; i < HASHSIZE; i++) edge[i] = NULL;
}

// given two nodes on an edge, add it to the hashtable
void edge_hashtab_add(EDGE_LIST_ITEM **edge, int nd1, int nd2, int *global_nd) {
    int nd_tmp;
    if (nd1 > nd2) {
        nd_tmp = nd1;
        nd1 = nd2;
        nd2 = nd_tmp;
    }
    if (global_nd != NULL) {
        if((global_nd[nd1]!=UNSET_INT) && (global_nd[nd2]!=UNSET_INT)) {
            edge_hash_add_entry(nd1, nd2, edge);
        }
    } else{
        edge_hash_add_entry(nd1, nd2, edge);
    }
}

// load edge hash table for a 3D columnar grid
void edge_hashtab_load_columnar(EDGE_LIST_ITEM **edge, SGRID *g, int *global_node) {
    int i, icol, ie;
#ifdef _DEBUG
    assert(g->ndim == 3);
    assert(g->type == COLUMNAR);
#endif
    
    for (icol=0; icol<g->ncolumns; icol++) {
        ie = g->elem2d_sur[icol];
        for (i = 0; i < NEDGEPRFC; i++) {
            edge_hashtab_add(edge, g->elem2d[ie].nodes[g->nd_on_TriEdge[i][0]], g->elem2d[ie].nodes[g->nd_on_TriEdge[i][1]],global_node);
        }
    }
}

// load edge hash table for a 3D unstructured (CJT :: do we need to add 2d external boundary edges here too?)
void edge_hashtab_load_unstructured(EDGE_LIST_ITEM **edge, SGRID *g, int *global_node) {
    int i, ie;
#ifdef _DEBUG
    assert(g->ndim == 3);
    assert(g->type == UNSTRUCTURED);
#endif
    for (ie=0; ie<g->nelems3d; ie++) {
        if (g->elem3d[ie].nnodes == 4) { // tet
            for (i=0; i<6; i++) {
                edge_hashtab_add(edge, g->elem3d[ie].nodes[g->nd_on_TetEdge[i][0]], g->elem3d[ie].nodes[g->nd_on_TetEdge[i][1]],global_node);
            }
        } else if (g->elem3d[ie].nnodes == 6) { // prism
            for (i=0; i<9; i++) {
                edge_hashtab_add(edge, g->elem3d[ie].nodes[g->nd_on_PrismEdge[i][0]], g->elem3d[ie].nodes[g->nd_on_PrismEdge[i][1]],global_node);
            }
        }
    }
}

// load edge hash table for a 2D grid
void edge_hashtab_load_2D(EDGE_LIST_ITEM **edge, SGRID *g, int *global_node) {
    int i, ie, icol, nd1, nd2, nd_tmp;
#ifdef _DEBUG
    assert(g->ndim == 2);
#endif
    // 2D element edges
    for (ie=0; ie<g->nelems2d; ie++) {
        if (g->elem2d[ie].nnodes == 3) { // triangle
            for (i=0; i<3; i++) {
                edge_hashtab_add(edge, g->elem2d[ie].nodes[g->nd_on_TriEdge[i][0]], g->elem2d[ie].nodes[g->nd_on_TriEdge[i][1]],global_node);
            }
        } else if (g->elem2d[ie].nnodes == 4) { // quadrilateral
            for (i=0; i<4; i++) {
                edge_hashtab_add(edge, g->elem2d[ie].nodes[g->nd_on_QuadEdge[i][0]], g->elem2d[ie].nodes[g->nd_on_QuadEdge[i][1]],global_node);
            }
        }
    }
    
    // 1D element edges
    for (ie=0; ie<g->nelems1d; ie++){
        edge_hashtab_add(edge, g->elem1d[ie].nodes[0], g->elem1d[ie].nodes[1], global_node);
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void increment_edge_array(SGRID *g, int nd, int surface_nd, int local_nnode, int *nedges) {
    if (g->node[nd].resident_pe == g->smpi->myid) {
        if (surface_nd < 0 || surface_nd > local_nnode - 1) {
            printf("ERROR :: PE: %d node1 < 0 :: node: %d :: surface node: %d\n", g->smpi->myid, nd, surface_nd);
            tl_error("ERROR: partitioning!\n");
        }
        nedges[surface_nd]++;
    }
}

// calculate the total number of edges connected to each node
double count_node_edges(SGRID *g, int local_nnode, int *map2surface, EDGE_LIST_ITEM **edge_hashtab, int *nedges, bool column_flag) {
    int i, node1, node2, snode1, snode2;
    double delx,dely,delz,dist = 0., dist_max = 0.;
    EDGE_LIST_ITEM *ep; // points to an edge list item
    int myid = g->smpi->myid;
    
    for (i = 0; i < local_nnode; i++) nedges[i] = 0;
    
    for (i = 0; i < HASHSIZE; i++) {
        ep = edge_hashtab[i];
        while (ep != NULL) {
            node1 = ep->nd1; if (column_flag) {snode1 = map2surface[node1];} else {snode1=node1;}
            node2 = ep->nd2; if (column_flag) {snode2 = map2surface[node2];} else {snode2=node2;}
            delx = g->node[node1].x - g->node[node2].x;
            dely = g->node[node1].y - g->node[node2].y;
            delz = g->node[node1].z - g->node[node2].z;
            dist = sqrt(delx * delx + dely * dely + delz * delz);
            if (dist > dist_max) dist_max = dist;
            increment_edge_array(g,node1,snode1,local_nnode,nedges);
            increment_edge_array(g,node2,snode2,local_nnode,nedges);
            ep = ep->next;
        }
    }
    dist_max = messg_dmax(dist_max, g->smpi->ADH_COMM);
    return dist_max;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void weight_graph_nodes(SGRID *g, int ncolumns, idx_t *vwgt) {
    int icol, nd, icount = 0;
    ID_LIST_ITEM *ptr;
    for (icol = 0; icol < ncolumns; icol++) {
        ptr = g->vertical_list[icol];
        nd = ptr->id;
        if (g->node[nd].resident_pe == g->smpi->myid) {
            /* Count the number of nodes in this vertical line to get a weight */
            vwgt[icount] = 1;
            while (ptr->next != NULL) {
                vwgt[icount]++;
                ptr = ptr->next;
            }
            icount++;
        }
    }
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

// store edge/node connectivity
void store_node_edges(SGRID *g, int local_nnode, double dist_max, EDGE_LIST_ITEM **edge_hashtab, int *nedges, idx_t *adjncy, idx_t *adjwgt, idx_t *xadj, int column_flag, int *global_nd) {
    int i, node1, node2, snode1, snode2, iedge1, iedge2, nd_org1=UNSET_INT, nd_org2=UNSET_INT;
    double delx,dely,delz,dist = 0.;
    EDGE_LIST_ITEM *ep; // points to an edge list item
    int myid = g->smpi->myid;
    
    for (i = 0; i < local_nnode; i++) nedges[i] = 0;
    
    for (i = 0; i < HASHSIZE; i++) {
        ep = edge_hashtab[i];
        while (ep != NULL) {
            node1 = ep->nd1; nd_org1 = node1; if (column_flag) {snode1 = g->nodeID_3d_to_2d_sur[node1];} else {snode1=node1;}
            node2 = ep->nd2; nd_org2 = node2; if (column_flag) {snode2 = g->nodeID_3d_to_2d_sur[node2];} else {snode2=node2;}
            delx = g->node[node1].x - g->node[node2].x;
            dely = g->node[node1].y - g->node[node2].y;
            delz = g->node[node1].z - g->node[node2].z;
            dist = sqrt(delx * delx + dely * dely + delz * delz);
            if (g->node[node1].resident_pe == myid) {
                if (snode1 < 0) {
                    printf("ERROR: partition_form.c :: PE: %d node1 < 0 :: node1: %d \t snode1: %d \n", myid, node1, snode1);
                    tl_error("ERROR: partitioning!\n");
                }
                iedge1 = xadj[snode1] + nedges[snode1];
                adjwgt[iedge1] = (idx_t) (dist_max / dist);
                nedges[snode1]++;
                if (global_nd != NULL) {
                    adjncy[iedge1] = global_nd[node2];
                } else {
                    if (column_flag) {
                        adjncy[iedge1] = g->node[node2].global_surf_id;
                    } else {
                        adjncy[iedge1] = g->node[node2].gid;
                    }
                }
            }
            if (g->node[node2].resident_pe == myid) {
                if (column_flag) node2 = g->nodeID_3d_to_2d_sur[node2];
                if (snode2 < 0) {
                    printf("ERROR: partition_form.c :: PE: %d node2 < 0 :: node2: %d \t snode2: %d \n", myid, node2, snode2);
                    tl_error("ERROR: partitioning!\n");
                }
                iedge2 = xadj[snode2] + nedges[snode2];
                adjwgt[iedge2] = (idx_t) (dist_max / dist);
                nedges[snode2]++;
                if (global_nd != NULL) {
                    adjncy[iedge2] = global_nd[node1];
                } else {
                    if (column_flag) {
                        adjncy[iedge2] = g->node[node1].global_surf_id;
                    } else {
                        adjncy[iedge2] = g->node[node1].gid;
                    }
                }
            }
            ep = ep->next;
        }
    }
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void partition_form_final(SGRID *g) {
    int i, icol, ierr_code = UNSET_INT;
    EDGE_LIST_ITEM *edge_hashtab[HASHSIZE];   // hash table for the edges
    int npes = g->smpi->npes;
    int myid = g->smpi->myid;
    
    bool column_flag = FALSE;
    if (g->ndim == 3 && g->type == COLUMNAR) column_flag = TRUE;
    //assert(g->type == COLUMNAR);
    //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
    
    // Count the number of owned nodes
    int local_nnode = g->my_nnodes;
    if (column_flag) {
        local_nnode = g->my_nnodes_sur;
    }
    
    // flag to indicate graph weighting
    idx_t wgtflag = 1;
    if (column_flag) wgtflag = 3;

    // allocate the memory for the global node numbers
    int *nnode_pe = (int *) tl_alloc(sizeof(int), npes);
    idx_t *vtxdist = (idx_t *) tl_alloc(sizeof(idx_t), npes + 1); // the number of nodes belonging to each pe
    
    ierr_code = MPI_Allgather(&(local_nnode), 1, MPI_INT, nnode_pe, 1, MPI_INT, g->smpi->ADH_COMM);
    if (ierr_code != MPI_SUCCESS) messg_err(ierr_code);
    
    // gather number of actual nodes on each processor to all processors
    vtxdist[0] = 0;
    for (i = 1; i <= npes; i++){
        vtxdist[i] = vtxdist[i - 1] + nnode_pe[i - 1];
    }
    
    // graph node weights
    idx_t *vwgt = NULL;
    if (column_flag) {
        vwgt = (idx_t *) tl_alloc(sizeof(idx_t), local_nnode);
        weight_graph_nodes(g, local_nnode, vwgt);
    }
    
    idx_t ncon=1; //number of weights per vertex
    real_t *tpwgts = (real_t *) tl_alloc(sizeof(real_t), npes);  // fraction of vertex weight to each part ncon x npart array
    real_t *ubvec  = (real_t *) tl_alloc(sizeof(real_t), ncon);  // imbalance tolerance for each ncon weight
    ubvec[0] = 1.05;
    for (i = 0; i < npes; i++) tpwgts[i] = 1.0/npes;

    // initialize the hash table
    //tag(MPI_COMM_WORLD);
    edge_hashtab_init(edge_hashtab);
    //tag(MPI_COMM_WORLD);
    
    // load nodal connections in hash table to remove redundancies in element connections
    if (g->ndim == 2) { // 2D
        edge_hashtab_load_2D(edge_hashtab,g,NULL);
    } else { // 3D
        if (column_flag) {
            edge_hashtab_load_columnar(edge_hashtab,g,NULL);
        } else {
            edge_hashtab_load_unstructured(edge_hashtab,g,NULL);
        }
    }
    //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
                                           
    // count total number of connections and number of connections per node
    int *nedges = (int *) tl_alloc(sizeof(int), local_nnode);     // the number of nodal connections for each node
    double dist_max = 0;
    int nedge_total = 0;
    dist_max = count_node_edges(g,local_nnode,g->nodeID_3d_to_2d_sur,edge_hashtab,nedges,column_flag);
    
    //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
    
    idx_t *xadj = (idx_t *) tl_alloc(sizeof(idx_t), local_nnode + 1); // the beginning of each nodes list of edges
    for (i = 0, nedge_total = 0, xadj[0] = 0; i < local_nnode; i++) {
        nedge_total += nedges[i];
        xadj[i + 1] = nedge_total;
    }
    
    if (DEBUG == 1) printf("partition_form: pe %d: # edges %d\n", myid, nedge_total);
    
    //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
    
    // load the connections in the adjacency list array
    idx_t *adjncy = (idx_t *) tl_alloc(sizeof(idx_t), nedge_total); // the nodal connection table
    idx_t *adjwgt = (idx_t *) tl_alloc(sizeof(idx_t), nedge_total); // edge weights for the graph
    store_node_edges(g,local_nnode,dist_max,edge_hashtab,nedges,adjncy,adjwgt,xadj,column_flag,NULL);
    
    //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
    
    // partitioning
    idx_t nparts = npes;
    idx_t *metis_part = (idx_t *) tl_alloc(sizeof(idx_t), local_nnode); // the partition returned from metis
    MPI_Comm comm;
    MPI_Comm_dup(g->smpi->ADH_COMM, &comm);
    
    //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
    
    /* CJT :: debug
     char fn[30+1];
     snprintf(fn, 30, "partition_info_pe%d.dat", g->smpi->myid);
     FILE *fp = fopen(fn, "w");
     fprintf(fp,"local_nnode: %d \t nedge_total: %d\n",local_nnode,nedge_total);
     for (i=0; i<nedge_total; i++) {fprintf(fp,"edge: %d \t adjncy: %d \t adjwgt: %d\n",i,adjncy[i],adjwgt[i]);}
     for (i=0; i<local_nnode; i++) {fprintf(fp,"local_nnode: %d \t vwgt: %d \t xadj: %d\n",i,vwgt[i],xadj[i]);}
     for (i=0; i<npes+1; i++) {fprintf(fp,"pe: %d \t vtxdist: %d \n",i,vtxdist[i]);}
     for (i=0; i<npes; i++) {fprintf(fp,"pe: %d \t nnode_pe: %d\n",i,nnode_pe[i]);}
     fflush(fp);
     fclose(fp);
     tl_error("for now");
     */
    
    idx_t edgecut;                       // number of edge cuts (returned)
    idx_t numflag = 0;                   // flag to indicate numbering scheme
    idx_t options[4] = { 0, 2, 150, 0 }; // the options for metis
    ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, options, &edgecut, metis_part, &comm);
    
    //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Comm_free(&comm);
    
    
    if (DEBUG) printf("partition_form: nparts = %d, edgecut = %d\n", nparts, edgecut);
    
    // cast the metis partition into part
    int k=0, nd=UNSET_INT;
    ID_LIST_ITEM *ptr;
    for (i = 0; i < local_nnode; i++) {
        if(column_flag){
            ptr = g->vertical_list[i];
            g->smpi->surface_partition_info[i] = (int)metis_part[i];
            while (ptr->next != NULL) {
                nd=ptr->id;
                g->smpi->partition_info[nd] = (int)metis_part[i];
                ptr = ptr->next;
            }
            if( g->smpi->surface_partition_info[i] != g->smpi->myid) k++;
        }
        else{
            g->smpi->partition_info[i] = (int)metis_part[i];
            if( g->smpi->partition_info[i] != g->smpi->myid) {k++;}
        }
    }
    
    metis_part = (idx_t *) tl_free(sizeof(idx_t), local_nnode, metis_part);
    nnode_pe =     (int *) tl_free(sizeof(int), npes, nnode_pe);
    vtxdist =    (idx_t *) tl_free(sizeof(idx_t), npes + 1, vtxdist);
    xadj =       (idx_t *) tl_free(sizeof(idx_t), local_nnode + 1, xadj);
    nedges =     (idx_t *) tl_free(sizeof(idx_t), local_nnode, nedges);
    adjncy =     (idx_t *) tl_free(sizeof(idx_t), nedge_total, adjncy);
    adjwgt =     (idx_t *) tl_free(sizeof(idx_t), nedge_total, adjwgt);
    tpwgts =    (real_t *) tl_free(sizeof(real_t), npes, tpwgts);
    ubvec =     (real_t *) tl_free(sizeof(real_t), ncon, ubvec);
    
    if (vwgt != NULL) vwgt = (idx_t *) tl_free(sizeof(idx_t), local_nnode, vwgt);
    tl_list_free_all(EDGE_LIST);
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void partition_adpt_final(SGRID *g) {
    
    int i,j,k,nd;
    EDGE_LIST_ITEM *ep;                       // points to an edge list item
    EDGE_LIST_ITEM *edge_hashtab[HASHSIZE];   // hash table for the edges
    
    int *nnode_pe = NULL;     // the number of nodes owned by each pe
    idx_t *adjncy = NULL;     // the nodal connection table
    idx_t *adjwgt = NULL;     // edge weights for the graph
    idx_t *vtxdist = NULL;    // the number of nodes belonging to each pe
    real_t *tpwgts = NULL;    // ordinarily = 1/npes
    real_t *ubvec = NULL;     // not used
    idx_t *metis_part = NULL; // the partition returned from metis
    idx_t *vwgt = NULL;       // node weights for the graph
    idx_t wgtflag = 0;        // flag to indicate graph weighting
    idx_t numflag = 0;        // flag to indicate numbering scheme
    idx_t ncon = 1;           // number of balance constraints (>= 1)
    idx_t edgecut;            // number of edge cuts
    idx_t nparts;             // number of partitions to be made
    idx_t *xadj;              // the beginning of each nodes list of edges
    idx_t *vsize = NULL;      // size of vertices wrt distribution costs
    real_t itr = 1000.0;      // ratio of interproc communication to distribution time
    int ierr_code;            // the error code from an mpi call
    int *nedges;              // number of nodal connections for each node
    int nedge_total = 0;

    bool column_flag = false;
    int nnode = g->nnodes;
    int my_nnode = g->my_nnodes;
    if (g->ndim == 3 && g->type == COLUMNAR) {
        column_flag = true;
        nnode = g->nnodes_sur;
        my_nnode = g->my_nnodes_sur;
    }
    
    int options[4] = { 0, 1, 15, 0 };    // the options for metis
    if (column_flag) {
        wgtflag=3;
        options[1]=2;
        options[2]=150;
    }
    if(DEBUG == 1) {
        options[0] = 1;        /* assigning options[0] == 1 means user must define options [1,2] */
        options[1] = 1;        /* more debug output can be had by increasing this value -- see ParMetis source. */
        options[2] = 15;       /* default value of random seed */
    }
    
    SMPI *smpi; /* alias */
    if(g->part_smpi != NULL) {
        smpi=g->part_smpi;
    }else{
        smpi=g->smpi;
    }
    int npes = smpi->npes;
    int myid = smpi->myid;
    
    if(g->interface == 0){
        // allocate memory
        nparts = smpi->npes;
        nnode_pe = (int *)tl_alloc(sizeof(int), nparts);
        
        vtxdist = (idx_t *) tl_alloc(sizeof(idx_t), nparts + 1);
        tpwgts = (real_t *)tl_alloc(sizeof(real_t), nparts);
        // fprintf(fp,"tpwgts\n");
        for(i = 0; i < nparts; i++){
            tpwgts[i] = 1.0 / nparts;//fprintf(fp,"%f\n",tpwgts[i]);
        }
        ubvec = (real_t *)tl_alloc(sizeof(real_t), ncon);
        ubvec[0] = 1.05;
        metis_part = (idx_t *) tl_alloc(sizeof(idx_t), my_nnode);
        
        // gather number of actual nodes on each processor to all processors
        vwgt = (idx_t *) tl_alloc(sizeof(idx_t), my_nnode);
        
        ierr_code = MPI_Allgather(&my_nnode, 1, MPI_INT, nnode_pe, 1, MPI_INT, smpi->ADH_COMM);
        if(ierr_code != MPI_SUCCESS) messg_err(ierr_code);
        
        /* calculate global node numbers */
        vtxdist[0] = 0;
        for(i = 1, vtxdist[0] = 0; i <= smpi->npes; i++) vtxdist[i] = vtxdist[i - 1] + nnode_pe[i - 1];
    }
    
    ID_LIST_ITEM *ptr;
    int icount=0,icol;
    int *global_nd = (int *) tl_alloc(sizeof(int), g->nnodes);
    int *local_map = (int *) tl_alloc(sizeof(int), g->nnodes);
    sarray_init_value_int(global_nd, g->nnodes, UNSET_INT);
    sarray_init_value_int(local_map, g->nnodes, UNSET_INT);
    if(g->interface == 0){
        if(column_flag){
            for (icol = 0; icol < nnode; icol++) {
                ptr = g->vertical_list[icol];
                nd = ptr->id;
                if (g->node[nd].resident_pe == myid) {
                    /* Count the number of nodes in this vertical line to get a weight */
                    global_nd[nd] = vtxdist[smpi->myid] + icount;
                    local_map[nd] = icount;
                    vwgt[icount] = 1;
                    while (ptr->next != NULL) {
                        vwgt[icount]++;
                        ptr = ptr->next;
                    }
                    icount++;
                }
            }
        }else{
            for(i = 0; i < my_nnode; i++){
                global_nd[i] = vtxdist[smpi->myid] + i;
                local_map[i] = icount++;
            }
        }
    }
    comm_update_int(global_nd, 1,smpi);
    
    if(g->interface == 0){
        
        // load nodal connections in hash table to remove redundancies in element connections
        if (g->ndim == 2) { // 2D
            edge_hashtab_load_2D(edge_hashtab,g,global_nd);
        } else { // 3D
            if (column_flag) {
                edge_hashtab_load_columnar(edge_hashtab,g,global_nd);
            } else {
                edge_hashtab_load_unstructured(edge_hashtab,g,global_nd);
            }
        }
        
        // count edges and find maximum edge length
        int node1,node2,snode1,snode2;
        double delx,dely,delz,dist;
        xadj = (idx_t *) tl_alloc(sizeof(idx_t), my_nnode + 1);
        nedges = (int *)tl_alloc(sizeof(int), my_nnode);
        double dist_max = 0.;
        for (i = 0; i < my_nnode; i++) nedges[i] = 0;
        for (i = 0; i < HASHSIZE; i++) {
            ep = edge_hashtab[i];
            while (ep != NULL) {
                node1 = ep->nd1; if (column_flag) {snode1 = local_map[node1];} else {snode1=node1;}
                node2 = ep->nd2; if (column_flag) {snode2 = local_map[node2];} else {snode2=node2;}
                if (column_flag) {
                    delx = g->node[node1].x - g->node[node2].x;
                    dely = g->node[node1].y - g->node[node2].y;
                    delz = g->node[node1].z - g->node[node2].z;
                    dist = sqrt(delx * delx + dely * dely + delz * delz);
                    if (dist > dist_max) dist_max = dist;
                    increment_edge_array(g,node1,snode1,my_nnode,nedges);
                    increment_edge_array(g,node2,snode2,my_nnode,nedges);
                } else {
                    if(node1 < my_nnode) nedges[node1]++;
                    if(node2 < my_nnode) nedges[node2]++;
                }
                ep = ep->next;
            }
        }
        if (column_flag) dist_max = messg_dmax(dist_max, g->smpi->ADH_COMM);
        
        nedge_total = 0;
        int iedge1, iedge2;
        for(i = 0, nedge_total = 0, xadj[0] = 0; i < my_nnode; i++) {
            nedge_total += nedges[i];
            xadj[i + 1] = nedge_total; //fprintf(fp,"%d\n",xadj[i+1]);
        }
        //  fprintf(fp,"nedge_total %d\n", nedge_total);
        
        // load edges into the adjacency list array
        adjncy = (idx_t *) tl_alloc(sizeof(idx_t), nedge_total); // the nodal connection table
        adjwgt = (idx_t *) tl_alloc(sizeof(idx_t), nedge_total); // edge weights for the graph
        int surf_node1, surf_node2;
        for(i = 0; i < my_nnode; i++) nedges[i] = 0;
        for (i = 0; i < HASHSIZE; i++) {
            ep = edge_hashtab[i];
            while (ep != NULL) {
                node1 = ep->nd1;
                node2 = ep->nd2;
                if(g->ndim==3){
                    delx = g->node[node1].x - g->node[node2].x;
                    dely = g->node[node1].y - g->node[node2].y;
                    delz = g->node[node1].z - g->node[node2].z;
                    dist = sqrt(delx * delx + dely * dely + delz * delz);
                }
                
                if (g->node[node1].resident_pe == myid) {
                    if(column_flag) {surf_node1 = local_map[node1];}
                    else{surf_node1 = node1;}
                    
                    if (surf_node1 < 0) {
                        printf("ERROR: partition_form.c :: PE: %d surf_node1 < 0 :: node1: %d \t surf_node1: %d \n", myid, node1, surf_node1);
                        exit(EXIT_FAILURE);
                    }
                    iedge1 = xadj[surf_node1] + nedges[surf_node1];
                    adjncy[iedge1] = global_nd[node2];
                    if(column_flag){
                        adjwgt[iedge1] = (idx_t) (dist_max / dist);
                    }else{
                        adjwgt[iedge1] = 1;
                    }
                    
                    nedges[surf_node1]++;
                }
                
                if (g->node[node2].resident_pe == myid) {
                    if(column_flag){surf_node2 = local_map[node2];}
                    else{surf_node2 = node2;}
                    
                    if (surf_node2 < 0) {
                        printf("ERROR: partition_form.c :: PE: %d surf_node2 < 0 :: node2: %d \t surf_node2: %d \n", myid, node2, surf_node2);
                        exit(EXIT_FAILURE);
                    }
                    iedge2 = xadj[surf_node2] + nedges[surf_node2];
                    adjncy[iedge2] = global_nd[node1];
                    if(column_flag){
                        adjwgt[iedge2] = (idx_t) (dist_max / dist);
                    }else{
                        adjwgt[iedge2] = 1;
                    }
                    
                    nedges[surf_node2]++;
                }
                
                ep = ep->next;
            }
        }
        // deallocate memory consumed by hash table
        tl_list_free_all(EDGE_LIST);
        
        // partitioning
        for(i = 0; i < my_nnode; i++) {
            metis_part[i] = smpi->myid;
            if(column_flag == false) vwgt[i] = 1;
        }
        
        MPI_Comm comm;
        MPI_Comm_dup(smpi->ADH_COMM, &comm);
        ParMETIS_V3_AdaptiveRepart(vtxdist, xadj, adjncy, vwgt, vsize, adjwgt, &wgtflag, &numflag,
                                   &ncon, &nparts, tpwgts, ubvec, &itr, options, &edgecut,
                                   metis_part, &comm);
        MPI_Comm_free(&comm);
        
        k=0;
        for (i = 0; i < my_nnode; i++) {
            if(column_flag) {
                ptr = g->vertical_list[i];
                if(g->part_map != NULL){
                    smpi->surface_partition_info[i] = g->part_map[(int)metis_part[i]];
                }else{
                    smpi->surface_partition_info[i] = (int)metis_part[i];
                }
                
                while (ptr->next != NULL) {
                    nd=ptr->id;
                    if(g->part_map != NULL){
                        smpi->partition_info[nd] = g->part_map[(int)metis_part[i]];
                    }else{
                        smpi->partition_info[nd] = (int)metis_part[i];
                    }
                    ptr = ptr->next;
                }
                if( smpi->surface_partition_info[i] != smpi->myid) k++;
            }else{
                if(g->part_map != NULL){
                    smpi->partition_info[i] = g->part_map[(int)metis_part[i]];
                }else{
                    smpi->partition_info[i] = (int)metis_part[i];
                }
                if( smpi->partition_info[i] != smpi->myid) {k++;}
            }
        }
        
    } else {
        for (i = 0; i < my_nnode; i++) {
            if (column_flag) {
                ptr = g->vertical_list[i];
                smpi->surface_partition_info[i] = smpi->myid;
                while (ptr->next != NULL) {
                    nd=ptr->id;
                    smpi->partition_info[nd] = smpi->myid;
                    ptr = ptr->next;
                }
                if( smpi->surface_partition_info[i] != smpi->myid) k++;
            } else {
                smpi->partition_info[i] = smpi->myid;
                if( smpi->partition_info[i] != smpi->myid) {k++;}
            }
        }
    }

    if(DEBUG == 2) { /* write partition vectors to files */
        char filename[MAXLINE];
        FILE *partfile;
        sprintf(filename, "metis.part%d", myid);
        partfile = io_fopen(filename, "w", TRUE);
        fclose(partfile);
    }
    
    
    /* free memory */
    if(g->interface==0){
        metis_part = (idx_t *) tl_free(sizeof(idx_t), my_nnode, metis_part);
        nnode_pe = (int *)tl_free(sizeof(int), nparts, nnode_pe);
        vtxdist = (idx_t *) tl_free(sizeof(idx_t), nparts + 1, vtxdist);
        xadj = (idx_t *) tl_free(sizeof(idx_t), my_nnode + 1, xadj);
        nedges = (idx_t *)tl_free(sizeof(idx_t), my_nnode, nedges);
        adjncy = (idx_t *) tl_free(sizeof(idx_t), nedge_total, adjncy);
        adjwgt = (idx_t *) tl_free(sizeof(idx_t), nedge_total, adjwgt);
        vwgt = (idx_t *) tl_free(sizeof(idx_t), my_nnode, vwgt);
        tpwgts = (real_t *)tl_free(sizeof(real_t), nparts, tpwgts);
        ubvec = (real_t *)tl_free(sizeof(real_t), 1, ubvec);
        tl_list_free_all(EDGE_LIST);
    }
    
    global_nd = (int *) tl_free(sizeof(int), g->nnodes, global_nd);
    local_map = (int *) tl_free(sizeof(int), g->nnodes, local_map);
    
}
