/* calls */

#include "global_header.h"
#include "parmetis.h"
/* so far workds for 2D only */
void partition_adpt(SGRID *g) {
    
    idx_t *vtxdist;		/* the number of nodes belonging to each pe */
    idx_t *adjncy;		/* the nodal connection table */
    idx_t *xadj;		/* the beginning of each nodes list of edges */
    idx_t *vwgt = NULL;		/* node weights for the graph */
    idx_t *adjwgt = NULL;	/* edge weights for the graph */
    idx_t *vsize = NULL;	/* size of vertices wrt distribution costs */
    idx_t wgtflag = 0;		/* flag to indicate graph weighting */
    idx_t numflag = 0;		/* flag to indicate numbering scheme */
    idx_t ncon = 1;			/* number of balance constraints (>= 1) */
    
    idx_t edgecut;	/* number of edge cuts */
    idx_t nparts;			/* number of partitions to be made */
    real_t *tpwgts;		/* ordinarily = 1/npes */
    real_t *ubvec;			/* not used */
    real_t itr = 1000.0;		/* ratio of interproc communication to distribution time */
    int options[4] = { 0, 1, 15, 0 };	/* the options for metis */
    idx_t *metis_part;		/* the partition returned from metis */
    MPI_Comm comm;
    EDGE_LIST_ITEM *ep;		/* points to an edge list item */
    EDGE_LIST_ITEM *edge_hashtab[HASHSIZE];	/* the hash table for the edges */
    int *nedges;			/* number of nodal connections for each node */
    int *nnode_pe;		/* the number of nodes owned by each pe */
    int ie, i, icol, icount, k, ierr;			/* loop counters */
    int ierr_code;		/* the error code from an mpi call */
    int nedge_total;		/* the total number of edge entries */
    int nd1, nd2, nd_tmp;		/* the nodes on an edge */
    int node1, node2;		/* the nodes on an edge */
    int iedge1, iedge2;		/* the current edge (nd1, nd2) and (nd2, nd1) */
    double delx, dely, delz, dist, dist_max;
    int nnode = g->nnodes;
    int my_nnode = g->my_nnodes;
    
    FILE *partfile;
    char filename[MAXLINE];
    int verbose = 0;
    ID_LIST_ITEM *ptr;
    int nd, ncolumns, surf_node1, surf_node2 ;
    FILE *fp;
    
    SMPI *smpi; /* alias */
    if(g->part_smpi != NULL) {
        smpi=g->part_smpi;
    }else{
        smpi=g->smpi;
    }
    int npes = smpi->npes;
    int myid = smpi->myid;
    
    
    //fp = io_fopen(build_filename2(filename, MAXLINE, "PARMETIS_STUFF", "", myid, ".txt", UNSET_INT), "w", TRUE);
    if(g->ndim==3){
        wgtflag=3;
        options[1]=2;
        options[2]=150;
        nnode = g->nnodes_sur;
        my_nnode = g->my_nnodes_sur;
        ncolumns = g->ncolumns;
    }else{
        ncolumns = g->nelems2d;
    }
    
    if(g->interface==0){
        /* allocate memory */
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
        
        /* gather number of actual nodes on each processor to all processors */
        
        vwgt = (idx_t *) tl_alloc(sizeof(idx_t), my_nnode);
        
        
        ierr_code = MPI_Allgather(&my_nnode, 1, MPI_INT, nnode_pe, 1, MPI_INT, smpi->ADH_COMM);
        if(ierr_code != MPI_SUCCESS)
            messg_err(ierr_code);
        
        /* calculate global node numbers */
        vtxdist[0] = 0;
        for(i = 1, vtxdist[0] = 0; i <= smpi->npes; i++)
            vtxdist[i] = vtxdist[i - 1] + nnode_pe[i - 1];
    }
    
    icount=0;
    int *global_nd;
    int *local_map;
    global_nd = (int *) tl_alloc(sizeof(int), g->nnodes);
    local_map = (int *) tl_alloc(sizeof(int), g->nnodes);
    sarray_init_value_int(global_nd, g->nnodes, UNSET_INT);
    sarray_init_value_int(local_map, g->nnodes, UNSET_INT);
    if(g->interface==0){
        if(g->ndim==3){
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
    
    /* load edges in hash table */
    if(g->interface==0){
        for(i = 0; i < HASHSIZE; i++)
            edge_hashtab[i] = NULL;
        
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
                //edge_hash_add_entry(nd1, nd2, edge_hashtab);
                if((global_nd[nd1]!=UNSET_INT) && (global_nd[nd2]!=UNSET_INT)){
                    edge_hash_add_entry(nd1, nd2, edge_hashtab);
                }
            }
        }
        
        if(g->ndim ==2){
            for (ie = 0;ie < g->nelems1d; ie++){
                nd1 = g->elem1d[ie].nodes[0];
                nd2 = g->elem1d[ie].nodes[1];
                
                if (nd1 > nd2) {
                    nd_tmp = nd1;
                    nd1 = nd2;
                    nd2 = nd_tmp;
                }
                //edge_hash_add_entry(nd1, nd2, edge_hashtab);
                if((global_nd[nd1]!=UNSET_INT) && (global_nd[nd2]!=UNSET_INT)){
                    edge_hash_add_entry(nd1, nd2, edge_hashtab);
                }
            }
        }
        
        /* count edges and find maximum edge length */
        
        xadj = (idx_t *) tl_alloc(sizeof(idx_t), my_nnode + 1);
        nedges = (int *)tl_alloc(sizeof(int), my_nnode);
        
        for(i = 0; i < my_nnode; i++)
            nedges[i] = 0;
        
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
                    if (dist > dist_max)
                        dist_max = dist;
                    
                    if (g->node[node1].resident_pe == myid) {
                        if(g->ndim ==3){surf_node1 = local_map[node1];}
                        else{surf_node1 = node1;}
                        
                        if (surf_node1 < 0 || surf_node1 > my_nnode - 1) {
                            printf("ERROR: partition_form.c :: PE: %d surf_node1 < 0 :: node1: %d \t surf_node1: %d \n", myid, node1, surf_node1);
                            exit(EXIT_FAILURE);
                        }
                        nedges[surf_node1]++;
                    }
                    
                    if (g->node[node2].resident_pe == myid) {
                        if(g->ndim ==3){surf_node2 = local_map[node2];}
                        else{surf_node2 = node2;}
                        
                        if (surf_node2 < 0 || surf_node2 > my_nnode - 1) {
                            printf("ERROR: partition_form.c :: PE: %d surf_node2 < 0 :: node2: %d \t surf_node2: %d \n", myid, node2, surf_node2);
                            exit(EXIT_FAILURE);
                        }
                        nedges[surf_node2]++;
                    }
                }else{
                    node1 = ep->nd1;
                    node2 = ep->nd2;
                    if(node1 < my_nnode)
                        nedges[node1]++;
                    if(node2 < my_nnode)
                        nedges[node2]++;
                }
                ep = ep->next;
            }
        }
        if(g->ndim==3) dist_max = messg_dmax(dist_max, smpi->ADH_COMM);
        //fprintf(fp,"xadj\n");
        for(i = 0, nedge_total = 0, xadj[0] = 0; i < my_nnode; i++)
        {
            nedge_total += nedges[i];
            xadj[i + 1] = nedge_total;
            //fprintf(fp,"%d\n",xadj[i+1]);
        }
        //  fprintf(fp,"nedge_total %d\n", nedge_total);
        
        /* load edges into the adjacency list array */
        
        adjncy = (idx_t *) tl_alloc(sizeof(idx_t), nedge_total);
        adjwgt = (idx_t *) tl_alloc(sizeof(idx_t), nedge_total);
        
        
        for(i = 0; i < my_nnode; i++)
            nedges[i] = 0;
        
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
                    if(g->ndim ==3){surf_node1 = local_map[node1];}
                    else{surf_node1 = node1;}
                    
                    if (surf_node1 < 0) {
                        printf("ERROR: partition_form.c :: PE: %d surf_node1 < 0 :: node1: %d \t surf_node1: %d \n", myid, node1, surf_node1);
                        exit(EXIT_FAILURE);
                    }
                    iedge1 = xadj[surf_node1] + nedges[surf_node1];
                    adjncy[iedge1] = global_nd[node2];
                    if(g->ndim ==3){
                        adjwgt[iedge1] = (idx_t) (dist_max / dist);
                    }else{
                        adjwgt[iedge1] = 1;
                    }
                    
                    nedges[surf_node1]++;
                }
                
                if (g->node[node2].resident_pe == myid) {
                    if(g->ndim ==3){surf_node2 = local_map[node2];}
                    else{surf_node2 = node2;}
                    
                    if (surf_node2 < 0) {
                        printf("ERROR: partition_form.c :: PE: %d surf_node2 < 0 :: node2: %d \t surf_node2: %d \n", myid, node2, surf_node2);
                        exit(EXIT_FAILURE);
                    }
                    iedge2 = xadj[surf_node2] + nedges[surf_node2];
                    adjncy[iedge2] = global_nd[node1];
                    if(g->ndim ==3){
                        adjwgt[iedge2] = (idx_t) (dist_max / dist);
                    }else{
                        adjwgt[iedge2] = 1;
                    }
                    
                    nedges[surf_node2]++;
                }
                
                ep = ep->next;
            }
        }
        //  fprintf(fp,"adjncy\n");
        //for(i = 0; i < nedge_total; i++)fprintf(fp,"%d\n",adjncy[i]);
        //fclose(fp);
        
        
        /* deallocate memory consumed by hash table */
        tl_list_free_all(EDGE_LIST);
        
        /* partitioning */
        
        for(i = 0; i < my_nnode; i++)
        {
            metis_part[i] = smpi->myid;
            if(g->ndim==2) vwgt[i] = 1;
        }
        
        
        MPI_Comm_dup(smpi->ADH_COMM, &comm);
        if(verbose >= 1)
        {
            options[0] = 1;		/* assigning options[0] == 1 means user must define options [1,2] */
            options[1] = 1;		/* more debug output can be had by increasing this value -- see ParMetis source. */
            options[2] = 15;		/* default value of random seed */
        }
        global_nd = (int *) tl_free(sizeof(int), g->nnodes, global_nd);
        local_map = (int *) tl_free(sizeof(int), g->nnodes, local_map);
        ParMETIS_V3_AdaptiveRepart(vtxdist, xadj, adjncy, vwgt, vsize, adjwgt, &wgtflag, &numflag,
                                   &ncon, &nparts, tpwgts, ubvec, &itr, options, &edgecut,
                                   metis_part, &comm);
        
        MPI_Comm_free(&comm);
        
        k=0;
        for (i = 0; i < my_nnode; i++) {
            if(g->ndim ==3){
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
        ierr=0;
        
    } else {
        for (i = 0; i < my_nnode; i++) {
            if(g->ndim ==3){
                ptr = g->vertical_list[i];
                smpi->surface_partition_info[i] = smpi->myid;
                
                while (ptr->next != NULL) {
                    nd=ptr->id;
                    
                    smpi->partition_info[nd] = smpi->myid;
                    ptr = ptr->next;
                }
                if( smpi->surface_partition_info[i] != smpi->myid) k++;
            }else{
                smpi->partition_info[i] = smpi->myid;
                if( smpi->partition_info[i] != smpi->myid) {k++;}
            }
        }
    }
    
    
    
    if(verbose == 2) { /* write partition vectors to files */
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

