// CJT :: I hate this file, can we just use read_geo??
// ONLY FOR TETS AND TRIS, even though it looks like both!!!!

/* this routine reads the geometry file */
/* and allocates interface nodes for parallel*/ 
/* strong coupled models*/
#include "global_header.h"
#define DIFF_TOL .0000001

int num_ghosts;

// returns the dimensionality of the grid
int read_interface_geo(SIO *io, SGRID *grid) {
    
    char line[MAXLINE];           /* the input line */
    char *data = NULL;            /* the data */
    char *subdata = NULL;            /* the data */
    NODE_LIST_ITEM *node_list, *node_list_tmp, *ptr_node, *ghost_nodes;
    int my_nnodes,nnodes, macro_nnodes, macro_nnodes_sur;
    double x;
    double y;
    double z;
    double surf_x, surf_y;
    int id, ninterfaces, add_nodes, ielem, ifc, bflag, istring, add_elem;
    int ind1,ind2, check, index, inode;
    int i,j,k,nd1,nd2,nd3,nd4;
    my_nnodes=0;
    nnodes=0;
    macro_nnodes = 0;
    macro_nnodes_sur = 0;
    int *interface_string, *interface_face, *interface_elements, my_node[4], local_id[4];
    int iel_mat, nelems2d, nelems3d, surf_id, local_index, my_nnode_surf, global_surf_id=0;
    ELEM2D_LIST_ITEM *elem2d_list, *ptr_elem2d;
    ELEM3D_LIST_ITEM *elem3d_list, *ptr_elem3d;
    NODE_LIST_ITEM *ptr;
    SVECT coords;
    CARD card;
    int nd_on_fc[NFCPRELM][NDPRFC] = { {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1} };
    
    int myid, ierr;
#ifdef _MESSG
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif
    add_elem=0;
    
    FILE *fp = NULL;
    char *filename = NULL;
    if (io->geo2d.fp != NULL) {
        grid->ndim = 2;
        fp = io->geo2d.fp;
        assert(io->geo2d.fp != NULL);
        filename = io->geo2d.filename;
    } else if (io->geo3d.fp != NULL) {
        grid->ndim = 3;
        fp = io->geo3d.fp;
        assert(io->geo3d.fp != NULL);
        filename = io->geo3d.filename;
    } else {
        tl_error("FATAL ERROR :: 2d or 3d geo file not found!");
    }
    
    assert(fp);
    assert(io->bc.fp);
    if (grid->ndim==3) assert(io->face.fp);

    rewind(fp);
    rewind(io->bc.fp);
    if (grid->ndim==3) rewind(io->face.fp);
    
    /* determine the grid dimension */
    while (fgets(line, MAXLINE, fp) != NULL) {
        io_save_line(io, fp, filename, line);
        if (strip_comments(line) <= 1) {    /* ignore empty ('\n') or comment lines */
            continue;
        }
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
            case CARD_E4T:
            case CARD_TET:
            case CARD_PRISM:
                grid->ndim = 3;
                id = read_int_field(*io, &data);
                grid->macro_nelems3d=id;
                break;
            case CARD_E3T:
            case CARD_TRI:
            case CARD_QUAD:
                grid->ndim = 2;
                id = read_int_field(*io, &data);
                grid->macro_nelems2d=id;
                break;
            case CARD_ND:
                macro_nnodes++;
                
                id = read_int_field(*io, &data);
                x = read_dbl_field(*io, &data);
                y = read_dbl_field(*io, &data);
                z = read_dbl_field(*io, &data);
                if (((fabs(fabs(x) - fabs(surf_x)) > DIFF_TOL) || (fabs(fabs(y) - fabs(surf_y)) > DIFF_TOL)) || grid->macro_nnodes_sur == 0) {
                    surf_x=x;
                    surf_y=y;
                    macro_nnodes_sur++;
                }
                break;
            default:
                break;
        }
    }

    grid->macro_nnodes = macro_nnodes;
    grid->macro_nnodes_sur = macro_nnodes_sur;

    node_list_tmp = tl_alloc(sizeof(NODE_LIST_ITEM), 1);
    node_list_tmp->next = NULL;
    
    // read 1d elements
    io_save_line(io, NULL, "", "");
    rewind(io->bc.fp);
    while (fgets(line, MAXLINE, io->bc.fp) != NULL) {
        io_save_line(io, io->bc.fp, io->bc.filename, line);
        if (strip_comments(line) <= 1) {    /* ignore empty ('\n') or comment lines */
            continue;
        }
        switch (parse_card(line, &data)) {
            case CARD_EGS:
                (grid->macro_nelems1d)++;
                break;
            case CARD_MDS:
                (grid->macro_nelems1d)++;
                break;
            case CARD_FCS:
                printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
                tl_error(">> please list 2d faces for 3d grid in a separate face file.\n");
                break;
            default:
                break;
        }
    }
    
    // read FCS if found in BC found (legacy)
    io_save_line(io, NULL, "", "");
    rewind(io->bc.fp);
    grid->nelems1d=grid->macro_nelems1d;
    if (grid->macro_nelems3d > 0) {
        assert(io->face.fp);
        rewind(io->face.fp);
        (grid->macro_nelems2d) = 0;
        while (fgets(line, MAXLINE, io->face.fp) != NULL) {
            io_save_line(io, io->face.fp, io->face.filename, line);
            if (strip_comments(line) <= 1) {    /* ignore empty ('\n') or comment lines */
                continue;
            }
            switch (parse_card(line, &data)) {
                case CARD_FCS:
                    (grid->macro_nelems2d)++;
                    break;
                default:
                    break;
            }
        }
        io_save_line(io, NULL, "", "");
        rewind(io->face.fp);
    }

    /* count number of interfaces */
    ninterfaces=0;
    while (fgets(line, MAXLINE, io->bc.fp) != NULL) {
        io_save_line(io, io->bc.fp, io->bc.filename, line);
        if (strip_comments(line) <= 1) {
            continue;
        }
        card = parse_card(line, &data);
        switch (card) {
            case CARD_HY:
                switch (parse_card(data, &subdata)){
                    case CARD_INT:
                        ninterfaces++;
                        break;
                }
                break;
            default:
                //tl_error("ERROR: No BC END card");
                break;
        }
    }
    interface_string = (int *)tl_alloc(sizeof(int), ninterfaces);
    ninterfaces=0;

    /* read interface string numbers */
    rewind(io->bc.fp);
    while (fgets(line, MAXLINE, io->bc.fp) != NULL) {
        io_save_line(io, io->bc.fp, io->bc.filename, line);
        if (strip_comments(line) <= 1) {
            continue;
        }
        card = parse_card(line, &data);
        switch (card) {
            case CARD_HY:
                switch (parse_card(data, &subdata)){
                    case CARD_INT:
                        interface_string[ninterfaces++] = read_int_field(*io, &subdata);
                        break;
                }
                break;
            default:
                //tl_error("ERROR: No BC END card");
                break;
        }
    }
    rewind(io->bc.fp);

    /* first deal with 2D grid interface nodes by looking at 1D elements in the BC file */
    if(grid->ndim==2){
        /* count number of nodes on interfaces */
        add_nodes=0;
        while (fgets(line, MAXLINE, io->bc.fp) != NULL) {
            io_save_line(io, io->bc.fp, io->bc.filename, line);
            if (strip_comments(line) <= 1) {
                continue;
            }
            card = parse_card(line, &data);
            switch (card) {
                case CARD_EGS:
                    ind1 = read_int_field(*io, &data);
                    ind2 = read_int_field(*io, &data);
                    istring = read_int_field(*io, &data);
                    ind1--;
                    ind2--;
                    for (i=0;i<ninterfaces;i++){
                        if(istring==interface_string[i]){
                            check = search_ghost_node(node_list_tmp,ind1);
                            if(check<0){
                                add_node_tmp(&node_list_tmp,my_nnodes,ind1);
                                my_nnodes++;
                            }
                            check = search_ghost_node(node_list_tmp,ind2);
                            if(check<0){
                                add_node_tmp(&node_list_tmp,my_nnodes,ind2);
                                my_nnodes++;
                            }
                        }
                    }
                    break;
                default:
                    break;
            }
        }
        rewind(io->bc.fp);
    }else{
        
        while (fgets(line, MAXLINE, io->face.fp) != NULL) {
            io_save_line(io, io->face.fp, io->face.filename, line);
            if (strip_comments(line) <= 1) {    /* ignore empty ('\n') or comment lines */
                continue;
            }
            switch (parse_card(line, &data)) {
                case CARD_FCS:
                    ielem   = read_int_field(*io, &data);
                    ifc     = read_int_field(*io, &data);
                    bflag   = read_int_field(*io, &data);
                    istring = read_int_field(*io, &data);
                    for (i=0;i<ninterfaces;i++){
                        if (istring == interface_string[i]) {add_elem++;}
                    }
                    break;
                case CARD_TRI:
                    ielem   = read_int_field(*io, &data);
                    nd1     = read_int_field(*io, &data); // node 1
                    nd2     = read_int_field(*io, &data); // node 2
                    nd3     = read_int_field(*io, &data); // node 3
                    bflag   = read_int_field(*io, &data);
                    istring = read_int_field(*io, &data);
                    for (i=0;i<ninterfaces;i++){
                        if (istring == interface_string[i]) {add_elem++;}
                    }
                    break;
                case CARD_QUAD:
                    ielem   = read_int_field(*io, &data);
                    nd1     = read_int_field(*io, &data); // node 1
                    nd2     = read_int_field(*io, &data); // node 2
                    nd3     = read_int_field(*io, &data); // node 3
                    nd4     = read_int_field(*io, &data); // node 4
                    bflag   = read_int_field(*io, &data);
                    istring = read_int_field(*io, &data);
                    for (i=0;i<ninterfaces;i++){
                        if (istring == interface_string[i]) {add_elem++;}
                    }
                    break;
                default:
                    break;
            }
        }
        rewind(io->face.fp);

        interface_elements = (int *)tl_alloc(sizeof(int), add_elem);
        interface_face = (int *)tl_alloc(sizeof(int), add_elem);

        add_elem=0;
        while (fgets(line, MAXLINE, io->face.fp) != NULL) {
            io_save_line(io, io->face.fp, io->face.filename, line);
            if (strip_comments(line) <= 1) {    /* ignore empty ('\n') or comment lines */
                continue;
            }
            switch (parse_card(line, &data)) {
                case CARD_TRI:
                    ielem = read_int_field(*io, &data);     ielem--;
                    nd1 = read_int_field(*io, &data);       nd1--;
                    nd2 = read_int_field(*io, &data);       nd2--;
                    nd3 = read_int_field(*io, &data);       nd3--;
                    bflag = read_int_field(*io, &data);
                    istring = read_int_field(*io, &data);
                    for (i=0;i<ninterfaces;i++){
                        if(istring==interface_string[i]){
                            for (i=0; i<4; i++) { // CJT :: ASSUMES TET!!!! BUT WHAT IF A PRISM SURFACE/BOTTOM
                                if      (nd1 == grid->elem3d[ielem].nodes[nd_on_fc[i][0]] &&
                                         nd2 == grid->elem3d[ielem].nodes[nd_on_fc[i][1]] &&
                                         nd3 == grid->elem3d[ielem].nodes[nd_on_fc[i][2]]) {ifc = i;}
                            }
                            interface_elements[add_elem]=ielem;
                            interface_face[add_elem++]=ifc;
                        }
                    }
                    
                case CARD_FCS:
                    ielem = read_int_field(*io, &data);     ielem--;
                    ifc = read_int_field(*io, &data);       ifc--;
                    bflag = read_int_field(*io, &data);
                    istring = read_int_field(*io, &data);
                    for (i=0;i<ninterfaces;i++){
                        if(istring==interface_string[i]){
                            interface_elements[add_elem]=ielem;
                            interface_face[add_elem++]=ifc;
                        }
                    }
                    break;
                default:
                    break;
            }
        }

        /* rewind file pointer; go back and read the data in */
        rewind(fp);
        while (fgets(line, MAXLINE, fp) != NULL){
            io_save_line(io, fp, filename, line);
            if (strip_comments(line) <= 1) continue;
            switch (parse_card(line, &data)) {
                case CARD_PRISM:
                    tl_error(">> prisms not supported for model coupling at the moment");
                    break;
                case CARD_E4T:
                case CARD_TET:
                    /* read element data */
                    index = read_int_field(*io, &data);
                    my_node[0] = read_int_field(*io, &data);
                    my_node[1]= read_int_field(*io, &data);
                    my_node[2] = read_int_field(*io, &data);
                    my_node[3] = read_int_field(*io, &data);
                    /* convert to zero base */
                    index--;
                    my_node[0]--;
                    my_node[1]--;
                    my_node[2]--;
                    my_node[3]--;
                    
                    for (i=0;i<add_elem;i++){
                        if(index==interface_elements[i]){
                            for(j=0;j<3;j++){
                                check = search_ghost_node(node_list_tmp,my_node[nd_on_fc[interface_face[i]][j]]);
                                if(check<0){
                                    add_node_tmp(&node_list_tmp,my_nnodes,my_node[nd_on_fc[interface_face[i]][j]]);
                                    my_nnodes++;
                                }
                            }
                        }
                    }
                    break;
                case CARD_END:
                    break;
                default:
                    break;
            }
        }
    }

    if(grid->ndim==3){
        interface_elements = (int *)tl_free(sizeof(int), add_elem, interface_elements);
        interface_face = (int *)tl_free(sizeof(int), add_elem, interface_face);
    }

    /* Now we know the interface nodes and can read the geo file.
     * All interface nodes and elements will be alocated and the create_interface
     * routine will assign ownership */
    elem2d_list = tl_alloc(sizeof(ELEM2D_LIST_ITEM), 1);
    elem2d_list->next = NULL;
    if (grid->ndim==3) {
        elem3d_list = tl_alloc(sizeof(ELEM3D_LIST_ITEM), 1);
        elem3d_list->next = NULL;
    }
    
    node_list = tl_alloc(sizeof(NODE_LIST_ITEM), 1);
    node_list->next = NULL;
    
    ghost_nodes = tl_alloc(sizeof(NODE_LIST_ITEM), 1);
    ghost_nodes->local=UNSET_INT;
#ifdef _MESSG
    ghost_nodes->global=UNSET_INT;
#endif
    ghost_nodes->next = NULL;
    
    num_ghosts = 0;

    /* rewind file pointer; go back and read the data in */
    rewind(fp);
    grid->nelems2d = 0;
    grid->nelems3d = 0;
    
    while (fgets(line, MAXLINE, fp) != NULL)
    {
        io_save_line(io, fp, filename, line);
        if (strip_comments(line) <= 1) /* ignore empty ('\n') or comment lines */
        {
            continue;
        }
        switch (parse_card(line, &data))
        {
            case CARD_E3T:
            case CARD_TRI:
                /* read element data */
                index = read_int_field(*io, &data);
                my_node[0] = read_int_field(*io, &data);
                my_node[1] = read_int_field(*io, &data);
                my_node[2] = read_int_field(*io, &data);
                iel_mat = read_int_field(*io, &data);
                
                /* convert to zero base */
                index--;
                my_node[0]--;
                my_node[1]--;
                my_node[2]--;
                iel_mat--;
                for (i=0;i<3;i++){
                    local_id[i] = search_ghost_node(node_list_tmp,my_node[i]);
                }
                if((local_id[0] < 0) && (local_id[1] < 0) && (local_id[2] < 0))
                    break;
                for(i=0;i<3;i++){
                    if(local_id[i]<0){
                        local_id[i]=search_ghost_node(ghost_nodes, my_node[i]);
                        if(local_id[i]<0){
                            local_id[i] = add_ghost_node(&ghost_nodes,my_node[i],my_nnodes,&num_ghosts);
                        }
                    }
                }
                add_elem2d_list(&elem2d_list, index, local_id, iel_mat);
                grid->nelems2d++;
                break;
            case CARD_E4T:
            case CARD_TET:
                /* read element data */
                index = read_int_field(*io, &data);
                my_node[0] = read_int_field(*io, &data);
                my_node[1]= read_int_field(*io, &data);
                my_node[2] = read_int_field(*io, &data);
                my_node[3] = read_int_field(*io, &data);
                iel_mat = read_int_field(*io, &data);
                iel_mat--;
                
                /* convert to zero base */
                index--;
                my_node[0]--;
                my_node[1]--;
                my_node[2]--;
                my_node[3]--;
                for (i=0;i<4;i++){
                    local_id[i] = search_ghost_node(node_list_tmp,my_node[i]);
                }
                if((local_id[0] < 0) && (local_id[1] < 0) && (local_id[2] < 0) && (local_id[3] < 0))
                    break;
                for(i=0;i<4;i++){
                    if(local_id[i]<0){
                        local_id[i]=search_ghost_node(ghost_nodes, my_node[i]);
                        if(local_id[i]<0){
                            local_id[i] = add_ghost_node(&ghost_nodes,my_node[i],my_nnodes,&num_ghosts);
                        }
                    }
                }
                add_elem3d_list(&elem3d_list, index, local_id, iel_mat);
                grid->nelems3d++;
                break;
            case CARD_END:
                break;
            default:
                break;
        }
    }
    io_save_line(io, NULL, "", "");
    rewind(fp);

    while (fgets(line, MAXLINE, fp) != NULL)
    {
        io_save_line(io, fp, filename, line);
        if (strip_comments(line) <= 1) /* ignore empty ('\n') or comment lines */
        {
            continue;
        }
        switch (parse_card(line, &data))
        {
            case CARD_ND:
                /*read node data */
                index = read_int_field(*io, &data);
                coords.x = read_dbl_field(*io, &data);
                coords.y = read_dbl_field(*io, &data);
                coords.z = read_dbl_field(*io, &data);
                if ( (grid->macro_nelems2d>0) && (index>(grid->macro_nelems2d * NDONTRI)) ) {
                    tl_error("Something is wrong with the node IDs.");
                }
                if ( (grid->macro_nelems3d>0) && (index>(grid->macro_nelems3d * NDONTET)) ) {
                    tl_error("Something is wrong with the node IDs.");
                }
                /* convert to zero-based ID */
                --index;
                
                /* create surface/bed node map */
                surf_id = UNSET_INT;
                local_index=search_ghost_node(node_list_tmp,index);
                if (((fabs(fabs(coords.x) - fabs(surf_x)) > DIFF_TOL) || (fabs(fabs(coords.y) - fabs(surf_y)) > DIFF_TOL)) || global_surf_id == 0) {
                    surf_x=coords.x;
                    surf_y=coords.y;
                    if(local_index>=0){
                        my_nnode_surf++;
                    }
                    surf_id = global_surf_id;
                    global_surf_id++;
                }
                /* check to see if I own this node */
                if(local_index>=0){
                    
                    add_node(&node_list, local_index, index, coords, -1, local_index, surf_id);
                }
                else if ((local_index = search_ghost_node(ghost_nodes, index)) >= 0)
                {
                    add_node(&node_list, local_index, index, coords, UNSET_INT, UNSET_INT, surf_id);
                }
                
                break;
            case CARD_END:
                break;
            default:
                break;
        }
    }

    grid->my_nnodes = my_nnodes;
    grid->nnodes = grid->my_nnodes + num_ghosts;
    grid->max_nnodes = grid->nnodes;
    
    /* allocate & initialize grid nodes */
    assert (grid->nnodes > 0);     /* must have > 0 nodes */
    snode_init_alloc_array(&(grid->node), grid->nnodes);

    /* allocate & initialize grid 1d elements */
    if (grid->nelems1d > 0) selem1d_init_alloc_array(&(grid->elem1d), grid->nelems1d, NDONSEG);
    
    /* allocate & intialize grid 2d elements */
    if (grid->nelems2d > 0 && grid->ndim == 2) selem2d_init_alloc_array(&(grid->elem2d), grid->nelems2d);
    
    /* allocate & intialize grid 3d elements */
    if (grid->nelems3d > 0) selem3d_init_alloc_array(&(grid->elem3d), grid->nelems3d);
    
    /* allocate & intialize grid 3d elements */
    if (grid->nelems3d > 0) {
        selem2d_init_alloc_array(&(grid->elem2d), grid->macro_nelems2d); // for now
    }

    ptr_node = node_list;
    while (ptr_node->next != NULL) {
        id = ptr_node->local;
        grid->node[id].id = id;
//#ifdef _MESSG
        grid->node[id].original_id = ptr_node->global;
        grid->node[id].x = ptr_node->coords.x;
        grid->node[id].y = ptr_node->coords.y;
        grid->node[id].z = ptr_node->coords.z;
        grid->node[id].gid = ptr_node->global;
        grid->node[id].global_surf_id = ptr_node->global_surf;
//#endif
        grid->node[id].string = NORMAL;
        grid->node[id].myid = -1;
        grid->node[id].resident_pe = ptr_node->sd;
        grid->node[id].resident_id = ptr_node->rnode;
        
        ptr_node = ptr_node->next;
    }

    if (grid->nelems2d > 0) {

        ptr_elem2d = elem2d_list;
        id = grid->nelems2d;
        while (ptr_elem2d->next != NULL) {
            id--;
            
            grid->elem2d[id].nnodes = NDONTRI;
            selem2d_alloc(&(grid->elem2d[id]), NDONTRI);
            grid->elem2d[id].nedges = 3;
            grid->elem2d[id].nnodes_quad = 6;
            grid->elem2d[id].edges = grid->nd_on_TriEdge;
            grid->haveTris = TRUE;
            
            grid->elem2d[id].id = id;
            grid->elem2d[id].gid = ptr_elem2d->ielem;
            grid->elem2d[id].mat = ptr_elem2d->mat;
            grid->elem2d[id].nodes[0] = ptr_elem2d->nd1;
            grid->elem2d[id].nodes[1] = ptr_elem2d->nd2;
            grid->elem2d[id].nodes[2] = ptr_elem2d->nd3;
            ptr_elem2d = ptr_elem2d->next;
        }
        
        if (id != 0) {
            printf("myid %d id %d nelem2sd %d \n", grid->smpi->myid, id, nelems2d);
            tl_error(">> Unexpected number of 2D elements.");
        }
    } else if (grid->nelems3d > 0) {
        
        /* fill in the elem3d array */
        ptr_elem3d = elem3d_list;
        id = grid->nelems3d;
        
        while (ptr_elem3d->next != NULL) {
            id--;
            
            grid->elem3d[id].nnodes = NDONTET;
            selem3d_alloc(&(grid->elem3d[id]), NDONTET);
            grid->elem3d[id].nedges = 6;
            grid->elem3d[id].nnodes_quad = 10;
            grid->elem3d[id].edges = grid->nd_on_TetEdge;
            grid->haveTets = TRUE;
            
            grid->elem3d[id].id = id;
            grid->elem3d[id].gid = ptr_elem3d->ielem;
            grid->elem3d[id].mat = ptr_elem3d->mat;
            grid->elem3d[id].nodes[0] = ptr_elem3d->nd1;
            grid->elem3d[id].nodes[1] = ptr_elem3d->nd2;
            grid->elem3d[id].nodes[2] = ptr_elem3d->nd3;
            grid->elem3d[id].nodes[3] = ptr_elem3d->nd4;
            ptr_elem3d = ptr_elem3d->next;
            
        }
    }
    
    io_save_line(io, NULL, "", "");
    //fclose(fp);

    /* read 2d faces if grid is 3d */
    if (grid->nelems3d > 0) {
        rewind(io->face.fp);
        read_faces(io, grid, NO);
        grid->elem2d = (SELEM_2D *) tl_realloc(sizeof(SELEM_2D), grid->nelems2d, grid->macro_nelems2d, grid->elem2d);
    }
    
    
    
    /* check for gaps in nodes*/
    for (i = 0; i < grid->nnodes; i++) {
        if (grid->node[i].string == UNSET_INT) {
            printf("MYID %d \n",grid->smpi->myid);
            printf("nnodes: %d \n",grid->nnodes);
            printf("i: %d\n",i);
            snode_printScreen(grid->node[i]);
            tl_error("Missing node in input geometry file.");
        }
    }
    /* check for gaps in elements */
    for (i = 0; i < grid->nelems2d; i++) {
        if (grid->elem2d[i].mat == UNSET_INT ) {
            tl_error("Missing 2d element in input geometry file.");
        }
    }
    for (i = 0; i < grid->nelems3d; i++) {
        if (grid->elem3d[i].mat == UNSET_INT) {
            tl_error("Missing 3d element in input geometry file.");
        }
    }
    
    // close the geo and face file
//    if (io->geo2d.fp != NULL) {
//        fclose(io->geo2d.fp);
//    } else if (io->geo3d.fp != NULL) {
//        fclose(io->geo3d.fp);
//    }
//    if (grid->nelems3d > 0) fclose(io->face.fp);
    
    
    /* Tidy up */
    if(grid->nelems3d > 0)  free_elem3d_list(elem3d_list);
    free_elem2d_list(elem2d_list);
    free_node_list(node_list);
    free_node_list(node_list_tmp);
    free_node_list(ghost_nodes);
    interface_string = (int *)tl_free(sizeof(int), ninterfaces, interface_string);
    return grid->ndim;
    
    
#ifdef _DEBUG
    fflush(stdout);
    //messg_barrier(MPI_COMM_WORLD);
    printf("-- GLOBAL PE: %d FINISHED Reading INTERFACE geo file %s\n",myid,filename);
#endif
}


