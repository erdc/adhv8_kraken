/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  partition_form.c This file partitions the grid for HPC usage.   Parmetis wrapper. */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"
#include "parmetis.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Paritions a model grid domain for distributed computing.  Wraps and calls Parmetis.
 *  \author    Corey Trahan, Ph.D.
 *  \author    Jeff Hensely, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout]  mod        (SMODEL *) a pointer to an AdH model structure
 * @param[in]  flag           (int) a flag to indicate = 1 no division of processors necessary or = 0 processors get divided in superfile read
 *
 * \note CJT\:: if there is a bug, it's all Jeff! Henslsey!
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static int DEBUG;

void partition_form(SGRID *g) {
    
    idx_t *vtxdist;             /* the number of nodes belonging to each pe */
    idx_t *adjncy;              /* the nodal connection table */
    idx_t *xadj;                /* the beginning of each nodes list of edges */
    idx_t *vwgt = NULL;         /* node weights for the graph */
    idx_t *adjwgt = NULL;       /* edge weights for the graph */
    idx_t wgtflag = 1;              /* flag to indicate graph weighting */
    idx_t numflag = 0;              /* flag to indicate numbering scheme */
    idx_t edgecut;                  /* number of edge cuts (returned) */
    idx_t nparts;                   /* number of partitions to be made */
    idx_t options[4] = { 0, 2, 150, 0 };    /* the options for metis */
    idx_t *metis_part;          /* the partition returned from metis */
    idx_t ncon=1;               /* number of weights per vertex */
    real_t *tpwgts;             /* fraction of vertex weight to each part ncon x npart array */
    real_t *ubvec;              /* imbalance tolerance for each ncon weight */
    MPI_Comm comm;
    
    
    EDGE_LIST_ITEM *ep;           /* points to an edge list item */
    EDGE_LIST_ITEM *edge_hashtab[HASHSIZE];   /* the hash table for the edges */
    
    int *global_nd;               /* the global node numbers of the nodes */
    int *nedges;                  /* the number of nodal connections for each node */
    int *nnode_pe;                /* the number of nodes owned by each pe */
    
    int ncolumns;
    
    int i, k, iproc, ierr, ierr2, ierr3;                     /* loop counters */
    int ie;                       /* loop counter over the elements */
    int ierr_code;                /* the error code from an mpi call */
    int nedge_total;              /* the total number of edge entries */
    int nd1, nd2, nd_tmp;         /* the nodes on an edge */
    int node1, node2;             /* the nodes on an edge */
    int iedge1, iedge2;           /* the current edge (nd1, nd2) and (nd2, nd1) */
    double delx, dely, delz, dist, dist_max;
    ID_LIST_ITEM *ptr;
    int icol;
    int nd;
    int icount = 0;
    int verbose = 0;              // cjt :: put back to 0
    int surf_node1, surf_node2;   /* cjt */
    int local_nnode=0;             /* cjt for sanity checking */
    int npes = g->smpi->npes;
    int myid = g->smpi->myid;
    
    
    /* Count the number of owned surface nodes */
    if(g->ndim == 2){
        local_nnode = g->my_nnodes;
        ncolumns = g->nelems2d;
    }
    else {
        wgtflag = 3;
        ncolumns = g->ncolumns;
        local_nnode=g->my_nnodes_sur;
    }
    /* allocate the memory for the global node numbers */
    nnode_pe = (int *) tl_alloc(sizeof(int), npes);
    vtxdist = (idx_t *) tl_alloc(sizeof(idx_t), npes + 1);
    
    ierr_code = MPI_Allgather(&(local_nnode), 1, MPI_INT, nnode_pe, 1, MPI_INT, g->smpi->ADH_COMM);
    if (ierr_code != MPI_SUCCESS) messg_err(ierr_code);
    
    /* gather number of actual nodes on each processor to all processors */
    vtxdist[0] = 0;
    for (i = 1; i <= npes; i++){
        vtxdist[i] = vtxdist[i - 1] + nnode_pe[i - 1];
    }
    
    if (g->ndim == 3) vwgt = (idx_t *) tl_alloc(sizeof(idx_t), local_nnode);
    tpwgts = (real_t *) tl_alloc(sizeof(real_t), npes);
    ubvec = (real_t *) tl_alloc(sizeof(real_t), ncon);
    ubvec[0] = 1.05;
    
    for (i = 0; i < npes; i++) tpwgts[i] = 1.0/npes;
    
    /* Number the local surface nodes */
    for (icol = 0; icol < local_nnode; icol++) {
        if(g->ndim ==3){
            ptr = g->vertical_list[icol];
            nd = ptr->id;
            if (g->node[nd].resident_pe == myid) {
                /* Count the number of nodes in this vertical line to get a weight */
                vwgt[icount] = 1;
                while (ptr->next != NULL) {
                    vwgt[icount]++;
                    ptr = ptr->next;
                }
                icount++;
            }
        }
        else { /*2D node based weighting can go here */
            //vwgt[icol] = 1;
            
        }
    }
    
    /* COMPUTE THE NODAL CONNECTIVITY INFORMATION */
    
    /* initialize the hash table */
    for (i = 0; i < HASHSIZE; i++)
        edge_hashtab[i] = NULL;
    
    /* load nodal connections in hash table to remove redundancies in element connections */
    for (icol = 0; icol < ncolumns; icol++) {
        if(g->ndim == 2){
            ie=icol;
        }
        else {
            ie = g->elem2d_sur[icol];
        }
        for (i = 0; i < NEDGEPRFC; i++) {
            nd1 = g->elem2d[ie].nodes[g->nd_on_TriEdge[i][0]];
            nd2 = g->elem2d[ie].nodes[g->nd_on_TriEdge[i][1]];
            
            if (nd1 > nd2) {
                nd_tmp = nd1;
                nd1 = nd2;
                nd2 = nd_tmp;
            }
            edge_hash_add_entry(nd1, nd2, edge_hashtab);
        }
    }
    
    for (ie = 0;ie < g->nelems1d; ie++){
        nd1 = g->elem1d[ie].nodes[0];
        nd2 = g->elem1d[ie].nodes[1];
        
        if (nd1 > nd2) {
            nd_tmp = nd1;
            nd1 = nd2;
            nd2 = nd_tmp;
        }
        edge_hash_add_entry(nd1, nd2, edge_hashtab);
    }
    
    /* count total number of connections and number of connections per node */
    xadj = (idx_t *) tl_alloc(sizeof(idx_t), local_nnode + 1);
    nedges = (int *) tl_alloc(sizeof(int), local_nnode);
    dist_max = 0;
    
    for (i = 0; i < local_nnode; i++)
        nedges[i] = 0;
    
    for (i = 0; i < HASHSIZE; i++) {
        ep = edge_hashtab[i];
        
        while (ep != NULL) {
            node1 = ep->nd1;
            node2 = ep->nd2;
            delx = g->node[node1].x - g->node[node2].x;
            dely = g->node[node1].y - g->node[node2].y;
            delz = g->node[node1].z - g->node[node2].z;
            dist = sqrt(delx * delx + dely * dely + delz * delz);
            if (dist > dist_max)
                dist_max = dist;
            
            if (g->node[node1].resident_pe == myid) {
                if(g->ndim ==3){surf_node1 = g->nodeID_3d_to_2d_sur[node1];}
                else{surf_node1 = node1;}
                
                if (surf_node1 < 0 || surf_node1 > local_nnode - 1) {
                    printf("ERROR: partition_form.c :: PE: %d surf_node1 < 0 :: node1: %d \t surf_node1: %d \n", myid, node1, surf_node1);
                    exit(EXIT_FAILURE);
                }
                nedges[surf_node1]++;
            }
            
            if (g->node[node2].resident_pe == myid) {
                if(g->ndim ==3){surf_node2 = g->nodeID_3d_to_2d_sur[node2];}
                else{surf_node2 = node2;}
                
                if (surf_node2 < 0 || surf_node2 > local_nnode - 1) {
                    printf("ERROR: partition_form.c :: PE: %d surf_node2 < 0 :: node2: %d \t surf_node2: %d \n", myid, node2, surf_node2);
                    exit(EXIT_FAILURE);
                }
                nedges[surf_node2]++;
            }
            
            ep = ep->next;
        }
    }
    
    dist_max = messg_dmax(dist_max, g->smpi->ADH_COMM);
    
    for (i = 0, nedge_total = 0, xadj[0] = 0; i < local_nnode; i++) {
        nedge_total += nedges[i];
        xadj[i + 1] = nedge_total;
    }
    
    
    if (verbose == 1) printf("partition_form: pe %d: # edges %d\n", myid, nedge_total);
    
    /* load the connections in the adjacency list array */
    adjncy = (idx_t *) tl_alloc(sizeof(idx_t), nedge_total);
    adjwgt = (idx_t *) tl_alloc(sizeof(idx_t), nedge_total);
    
    for (i = 0; i < local_nnode; i++)
        nedges[i] = 0;
    
    for (i = 0; i < HASHSIZE; i++) {
        ep = edge_hashtab[i];
        while (ep != NULL) {
            node1 = ep->nd1;
            node2 = ep->nd2;
            delx = g->node[node1].x - g->node[node2].x;
            dely = g->node[node1].y - g->node[node2].y;
            delz = g->node[node1].z - g->node[node2].z;
            dist = sqrt(delx * delx + dely * dely + delz * delz);
            
            if (g->node[node1].resident_pe == myid) {
                if(g->ndim ==3){surf_node1 = g->nodeID_3d_to_2d_sur[node1];}
                else{surf_node1 = node1;}
                
                if (surf_node1 < 0) {
                    printf("ERROR: partition_form.c :: PE: %d surf_node1 < 0 :: node1: %d \t surf_node1: %d \n", myid, node1, surf_node1);
                    exit(EXIT_FAILURE);
                }
                iedge1 = xadj[surf_node1] + nedges[surf_node1];
                
                if(g->ndim ==3){adjncy[iedge1] = g->node[node2].global_surf_id;}
                else{adjncy[iedge1] = g->node[node2].gid;}
                adjwgt[iedge1] = (int) (dist_max / dist);
                nedges[surf_node1]++;
            }
            
            if (g->node[node2].resident_pe == myid) {
                if(g->ndim ==3){surf_node2 = g->nodeID_3d_to_2d_sur[node2];}
                else{surf_node2 = node2;}
                
                if (surf_node2 < 0) {
                    printf("ERROR: partition_form.c :: PE: %d surf_node2 < 0 :: node2: %d \t surf_node2: %d \n", myid, node2, surf_node2);
                    exit(EXIT_FAILURE);
                }
                iedge2 = xadj[surf_node2] + nedges[surf_node2];
                
                if(g->ndim ==3){adjncy[iedge2] = g->node[node1].global_surf_id;}
                else{adjncy[iedge2] = g->node[node1].gid;}
                adjwgt[iedge2] = (int) (dist_max / dist);
                nedges[surf_node2]++;
            }
            
            ep = ep->next;
        }
    }
    
    /* partitioning */
    nparts = npes;
    metis_part = (idx_t *) tl_alloc(sizeof(idx_t), local_nnode);
    MPI_Comm_dup(g->smpi->ADH_COMM, &comm);
    
    /* CJT :: debug
     char fn[30+1];
     snprintf(fn, 30, "partition_info_pe%d.dat", g->smpi->myid);
     FILE *fp = fopen(fn, "w");
     fprintf(fp,"local_nnode: %d \t nedge_total: %d\n",local_nnode,nedge_total);
     for (i=0; i<nedge_total; i++) {
     fprintf(fp,"edge: %d \t adjncy: %d \t adjwgt: %d\n",i,adjncy[i],adjwgt[i]);
     }
     for (i=0; i<local_nnode; i++) {
     fprintf(fp,"local_nnode: %d \t vwgt: %d \t xadj: %d\n",i,vwgt[i],xadj[i]);
     }
     for (i=0; i<npes+1; i++) {
     fprintf(fp,"pe: %d \t vtxdist: %d \n",i,vtxdist[i]);
     }
     for (i=0; i<npes; i++) {
     fprintf(fp,"pe: %d \t nnode_pe: %d\n",i,nnode_pe[i]);
     }
     tl_error("for now");
     */
    
    ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, options, &edgecut, metis_part, &comm);
    
    MPI_Comm_free(&comm);
    
    if (verbose == 1 && myid == 0)
        printf("partition_form: nparts = %d, edgecut = %d\n", nparts, edgecut);
    
    /* cast the metis partition into part */
    k=0;
    for (i = 0; i < local_nnode; i++) {
        if(g->ndim ==3){
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
    ierr=0;
    
    metis_part = (idx_t *) tl_free(sizeof(idx_t), local_nnode, metis_part);
    nnode_pe = (int *) tl_free(sizeof(int), npes, nnode_pe);
    vtxdist = (idx_t *) tl_free(sizeof(idx_t), npes + 1, vtxdist);
    xadj = (idx_t *) tl_free(sizeof(idx_t), local_nnode + 1, xadj);
    nedges = (idx_t *) tl_free(sizeof(idx_t), local_nnode, nedges);
    adjncy = (idx_t *) tl_free(sizeof(idx_t), nedge_total, adjncy);
    adjwgt = (idx_t *) tl_free(sizeof(idx_t), nedge_total, adjwgt);
    if(g->ndim==3) vwgt = (idx_t *) tl_free(sizeof(idx_t), local_nnode, vwgt);
    tpwgts = (real_t *) tl_free(sizeof(real_t), npes, tpwgts);
    ubvec = (real_t *) tl_free(sizeof(real_t), ncon, ubvec);
    
    tl_list_free_all(EDGE_LIST);
    
}

