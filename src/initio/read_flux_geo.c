/* this routine reads the geometry file */
/* and allocates interface nodes for parallel*/ 
/* strong coupled models*/
#include "global_header.h"
#include "unistd.h"
#define NEWIO_HASHSIZE 10069
#define DIFF_TOL .0000001


int num_ghosts;

// returns the dimensionality of the grid
void read_flux_geo(int super, int sub, SIO *io, SGRID *grid) {
    
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
    
    SGRID *grid2;
    grid2=tl_alloc(sizeof(SGRID),1);
#ifdef _MESSG
    MPI_Comm dummy;
    sgrid_init(grid,0,dummy);
#else
    sgrid_init(grid,0);
#endif
    int myid, ierr;
#ifdef _MESSG
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif
    add_elem=0;
    
    /* determine the grid dimension */
    if( access( io->face.filename, F_OK ) != -1 ) {
        grid2->ndim = 3;
    }else{
        grid2->ndim=2;
    }
    
    SFILE_IN geo; // not sure if shallow copy will work with FILE *
    if (grid->ndim == 2) {
        strcpy(geo.filename, io->geo2d.filename);
        geo.fp = io->geo2d.fp;
    } else if (grid->ndim ==3) {
        strcpy(geo.filename, io->geo3d.filename);
        geo.fp = io->geo3d.fp;
    }
    
    node_list_tmp = tl_alloc(sizeof(NODE_LIST_ITEM), 1);
    node_list_tmp->next = NULL;
    
    rewind(geo.fp);
    rewind(io->bc.fp);
    assert(io->bc.fp);
    
    io_save_line(io, NULL, "", "");
    rewind(io->bc.fp);
    /* count number of flux interfaces */
    ninterfaces=0;
    while (fgets(line, MAXLINE, io->bc.fp) != NULL) {
        io_save_line(io, io->bc.fp, io->bc.filename, line);
        if (strip_comments(line) <= 1) {
            continue;
        }
        card = parse_card(line, &data);
        switch (card) {
            case CARD_NB:
                switch (parse_card(data, &subdata)){
                    case CARD_CPL:
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
            case CARD_NB:
                switch (parse_card(data, &subdata)){
                    case CARD_CPL:
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
    
    /* first deal with 2D grid interface nodes by looking
     * at 1D elements in the BC file */
    if(grid2->ndim==2){
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
        rewind(geo.fp);
        while (fgets(line, MAXLINE, geo.fp) != NULL){
            io_save_line(io, geo.fp, geo.filename, line);
            if (strip_comments(line) <= 1) /* ignore empty ('\n') or comment lines */
            {
                continue;
            }
            switch (parse_card(line, &data))
            {
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
    rewind(geo.fp);
    grid->nelems2d = 0;
    grid->nelems3d = 0;
    
    while (fgets(line, MAXLINE, geo.fp) != NULL)
    {
        io_save_line(io, geo.fp, geo.filename, line);
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
    rewind(geo.fp);
    
    while (fgets(line, MAXLINE, geo.fp) != NULL)
    {
        io_save_line(io, geo.fp, geo.filename, line);
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
                
                
                break;
            case CARD_END:
                break;
            default:
                break;
        }
    }
    grid2->my_nnodes = my_nnodes;
    grid2->nnodes = grid2->my_nnodes + num_ghosts;
    grid2->max_nnodes = grid2->nnodes;
    
    /* allocate & initialize grid nodes */
    assert (grid2->nnodes > 0);     /* must have > 0 nodes */
    grid2->node = (SNODE *) tl_alloc(sizeof(SNODE), grid2->nnodes);
    for (inode=0; inode < grid2->nnodes; inode++) {
        snode_init(&(grid2->node[inode]));
    }
    
    /* allocate & initialize grid 1d elements */
    if (grid2->nelems1d > 0) {
        grid2->elem1d = (SELEM_1D *) tl_alloc(sizeof(SELEM_1D), grid2->nelems1d);
        for (ielem=0; ielem < grid2->nelems1d; ielem++) {
            selem1d_init(&(grid2->elem1d[ielem]));
        }
    }
    
    /* allocate & intialize grid 2d elements */
    if (grid2->nelems2d > 0) {
        grid2->elem2d = (SELEM_2D *) tl_alloc(sizeof(SELEM_2D), grid2->nelems2d);
        for (ielem=0; ielem < grid2->nelems2d; ielem++) {
            selem2d_init(&(grid2->elem2d[ielem]));
        }
    }
    
    /* allocate & intialize grid 3d elements */
    if (grid2->nelems3d > 0) {
        grid2->elem3d = (SELEM_3D *) tl_alloc(sizeof(SELEM_3D), grid2->nelems3d);
        for (ielem=0; ielem<grid2->nelems3d; ielem++) {
            selem3d_init(&(grid2->elem3d[ielem]));
        }
        grid2->elem2d = (SELEM_2D *) tl_alloc(sizeof(SELEM_2D), grid2->macro_nelems2d); //all 2d elems for now, will reaclloc after faces file
        for (ielem=0; ielem < grid2->macro_nelems2d; ielem++) {
            selem2d_init(&(grid2->elem2d[ielem]));
        }
    }
    
    ptr_node = node_list;
    while (ptr_node->next != NULL)
    {
        id = ptr_node->local;
        grid2->node[id].id = id;
#ifdef _MESSG
        grid2->node[id].original_id = ptr_node->global; // cjt :: not sure if we need this line serially
        grid2->node[id].x = ptr_node->coords.x;
        grid2->node[id].y = ptr_node->coords.y;
        grid2->node[id].z = ptr_node->coords.z;
        grid2->node[id].gid = ptr_node->global;
        grid2->node[id].global_surf_id = ptr_node->global_surf;
#endif
        grid2->node[id].string = NORMAL;
        grid2->node[id].myid = -1;
        grid2->node[id].resident_pe = ptr_node->sd;
        grid2->node[id].resident_id = ptr_node->rnode;
        
        
        ptr_node = ptr_node->next;
    }
    
    if (grid2->nelems2d > 0)
    {
        
        ptr_elem2d = elem2d_list;
        id = grid2->nelems2d;
        
        while (ptr_elem2d->next != NULL)
        {
            id--;
            grid2->elem2d[id].id = id;
            grid2->elem2d[id].gid = ptr_elem2d->ielem;
            grid2->elem2d[id].mat = ptr_elem2d->mat;
            grid2->elem2d[id].nodes[0] = ptr_elem2d->nd1;
            grid2->elem2d[id].nodes[1] = ptr_elem2d->nd2;
            grid2->elem2d[id].nodes[2] = ptr_elem2d->nd3;
            ptr_elem2d = ptr_elem2d->next;
        }
        
        if (id != 0)
        {
            printf("myid %d id %d nelem2sd %d \n", grid2->smpi->myid, id, nelems2d);
            tl_error("Unexpected number of 2D elements.");
        }
    }
    else if (grid2->nelems3d > 0)
    {
        
        /* fill in the elem3d array */
        ptr_elem3d = elem3d_list;
        id = grid2->nelems3d;
        
        while (ptr_elem3d->next != NULL)
        {
            id--;
            grid2->elem3d[id].id = id;
            grid2->elem3d[id].gid = ptr_elem3d->ielem;
            grid2->elem3d[id].mat = ptr_elem3d->mat;
            grid2->elem3d[id].nodes[0] = ptr_elem3d->nd1;
            grid2->elem3d[id].nodes[1] = ptr_elem3d->nd2;
            grid2->elem3d[id].nodes[2] = ptr_elem3d->nd3;
            grid2->elem3d[id].nodes[3] = ptr_elem3d->nd4;
            ptr_elem3d = ptr_elem3d->next;
            
        }
    }
    
    io_save_line(io, NULL, "", "");
    fclose(geo.fp);
    
    /* read 2d faces if grid is 3d */
    if (grid2->nelems3d > 0) {
        rewind(io->face.fp);

        read_faces(io, grid, NO);
        grid2->elem2d = (SELEM_2D *) tl_realloc(sizeof(SELEM_2D), grid2->nelems2d, grid2->macro_nelems2d, grid2->elem2d);
    }
    /* read 1D elements for 2D grid flux */
    int local_index1,local_index2;
    if(grid2->ndim==2){
        grid2->nelems1d=0;
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
                            local_index1 = is_my_node(ind1, grid2);
                            local_index2 = is_my_node(ind2,grid2);
                            ind1 = local_index1;
                            ind2 = local_index2;
                            
                            if(grid2->nelems1d>=grid2->max_nelems1d) grid2->nelems1d = elem1d_new(grid, 40);
                            grid2->elem1d[grid2->nelems1d].id = grid2->nelems1d;
                            grid2->elem1d[grid2->nelems1d].string = istring;
                            grid2->elem1d[grid2->nelems1d].mat = UNSET_INT;
                            grid2->elem1d[grid2->nelems1d].nodes[0] = ind1;
                            grid2->elem1d[grid2->nelems1d].nodes[1] = ind2;
                        }
                    }
                    grid2->nelems1d++;
                    break;
                case CARD_MDS:
                    grid2->nelems1d++;
                    break;
                default:
                    break;
            }
        }
        rewind(io->bc.fp);
    }
    int overlay = 0, count, *cpl_elems, ie1,ie2;
    double x3d[NDONTRI], y3d[NDONTRI], z3d[NDONTRI],dist1, dist2, dist3;
    
    /* allocate temporary storage for coupled elements */
    if(grid2->ndim==3){
        cpl_elems = (int *) tl_alloc(sizeof(int),grid2->nelems2d);
    }else{
        cpl_elems = (int *) tl_alloc(sizeof(int),grid2->nelems1d);
    }
    /* Deremine flux coupled element numbers */
    if(grid->ndim==2){
        /* My grid is 2d */
        for(ie1=0;ie1<grid->nelems1d;ie1++){
            count = 0;
            if(grid2->ndim==3){
                /* for 2D-3D coupled grid is 3D, look at 2D faces */
                for (ie2 = 0; ie2 < grid2->nelems2d; ie2++) {
                    for (ifc=0;ifc<ninterfaces;ifc++){
                        if(grid2->elem2d[ie2].string==interface_string[ifc]){
                            
                            overlay = 0;
                            
                            for (i=0; i<NDONTRI; i++) {
                                x3d[i] = grid2->node[grid2->elem2d[ie2].nodes[i]].x;
                                y3d[i] = grid2->node[grid2->elem2d[ie2].nodes[i]].y;
                            }
                            
                            dist1 = fabs(grid->node[grid->elem1d[ie1].nodes[0]].x - x3d[0]) + fabs(grid->node[grid->elem1d[ie1].nodes[0]].y - y3d[0]);
                            dist2 = fabs(grid->node[grid->elem1d[ie1].nodes[0]].x - x3d[1]) + fabs(grid->node[grid->elem1d[ie1].nodes[0]].y - y3d[1]);
                            dist3 = fabs(grid->node[grid->elem1d[ie1].nodes[0]].x - x3d[2]) + fabs(grid->node[grid->elem1d[ie1].nodes[0]].y - y3d[2]);
                            if (dist1 < 1e-6 || dist2 < 1e-6 || dist3 < 1e-6) overlay++;
                            
                            dist1 = fabs(grid->node[grid->elem1d[ie1].nodes[1]].x - x3d[0]) + fabs(grid->node[grid->elem1d[ie1].nodes[1]].y - y3d[0]);
                            dist2 = fabs(grid->node[grid->elem1d[ie1].nodes[1]].x - x3d[1]) + fabs(grid->node[grid->elem1d[ie1].nodes[1]].y - y3d[1]);
                            dist3 = fabs(grid->node[grid->elem1d[ie1].nodes[1]].x - x3d[2]) + fabs(grid->node[grid->elem1d[ie1].nodes[1]].y - y3d[2]);
                            if (dist1 < 1e-6 || dist2 < 1e-6 || dist3 < 1e-6) overlay++;
                            
                            if (overlay == 2) { // two of the 3D faces nodes match different edge nodes, add
                                cpl_elems[count] = ie2;
                                count++;
                            }
                        }
                    }
                }
            }else{
                /* For 2D-2D coupled grid 2D, look at 1D elements */
                for (ie2 = 0; ie2 < grid2->nelems1d; ie2++) {
                    for (ifc=0;ifc<ninterfaces;ifc++){
                        if(grid2->elem1d[ie2].string==interface_string[ifc]){
                            
                            overlay = 0;
                            
                            for (i=0; i<2; i++) {
                                x3d[i] = grid2->node[grid2->elem1d[ie2].nodes[i]].x;
                                if(grid2->ndim==3){
                                    cpl_elems = (int *) tl_alloc(sizeof(int),grid2->nelems2d);
                                }else{
                                    cpl_elems = (int *) tl_alloc(sizeof(int),grid2->nelems1d);
                                }              y3d[i] = grid2->node[grid2->elem1d[ie2].nodes[i]].y;
                            }
                            
                            dist1 = fabs(grid->node[grid->elem1d[ie1].nodes[0]].x - x3d[0]) + fabs(grid->node[grid->elem1d[ie1].nodes[0]].y - y3d[0]);
                            dist2 = fabs(grid->node[grid->elem1d[ie1].nodes[0]].x - x3d[1]) + fabs(grid->node[grid->elem1d[ie1].nodes[0]].y - y3d[1]);
                            if (dist1 < 1e-6 || dist2 < 1e-6 ) overlay++;
                            
                            dist1 = fabs(grid->node[grid->elem1d[ie1].nodes[1]].x - x3d[0]) + fabs(grid->node[grid->elem1d[ie1].nodes[1]].y - y3d[0]);
                            dist2 = fabs(grid->node[grid->elem1d[ie1].nodes[1]].x - x3d[1]) + fabs(grid->node[grid->elem1d[ie1].nodes[1]].y - y3d[1]);
                            if (dist1 < 1e-6 || dist2 < 1e-6 ) overlay++;
                            
                            if (overlay == 2) { // two of the 3D faces nodes match different edge nodes, add
                                cpl_elems[count] = ie2;
                                count++;
                            }
                        }
                    }
                }
            }
            /* shouldn't need re_alloc but maybe more than 1 supermodel is coupled? */
            if(grid->elem1d[ie1].flux_elem !=NULL ){
                grid->elem1d[ie1].flux_elem = (FLUX_ELEM *) tl_realloc(sizeof(FLUX_ELEM), grid->elem1d[ie1].flux_elem_tot, (grid->elem1d[ie1].flux_elem_tot+count), grid->elem1d[ie1].flux_elem);
            }else{
                grid->elem1d[ie1].flux_elem = (FLUX_ELEM *) tl_alloc(sizeof(FLUX_ELEM), (grid->elem1d[ie1].flux_elem_tot+count));
            }
            for(i=0;i<count;i++){
                grid->elem1d[ie1].flux_elem[grid->elem1d[ie1].flux_elem_tot+count].elem  = cpl_elems[i];
                grid->elem1d[ie1].flux_elem[grid->elem1d[ie1].flux_elem_tot+count].sub  = sub;
                grid->elem1d[ie1].flux_elem[grid->elem1d[ie1].flux_elem_tot+count].super= super;
                printf("1d elem %d sup %d sub %d elem %d \n", ie1,super,sub,cpl_elems[i]);
            }
            grid->elem1d[ie1].flux_elem_tot = count;
        }
    }else{
        /* My grid is 3D look at 2D faces */
        for(ie1=0;ie1<grid->nelems2d;ie1++){
            count=0;
            if(grid2->ndim==3){
                /* For 3D-3D, coupled grid is 3D, match faces */
                for (ie2 = 0; ie2 < grid2->nelems2d; ie2++) {
                    for (ifc=0;ifc<ninterfaces;ifc++){
                        if(grid2->elem2d[ie2].string==interface_string[ifc]){
                            
                            overlay = 0;
                            
                            for (i=0; i<NDONTRI; i++) {
                                x3d[i] = grid2->node[grid2->elem2d[ie2].nodes[i]].x;
                                y3d[i] = grid2->node[grid2->elem2d[ie2].nodes[i]].y;
                                z3d[i] = grid2->node[grid2->elem2d[ie2].nodes[i]].z;
                            }
                            
                            for (i=0; i<NDONTRI; i++) {
                                /* compare againts each node to find exact face match*/
                                
                                dist1 = fabs(grid->node[grid->elem2d[ie1].nodes[i]].x - x3d[0])
                                + fabs(grid->node[grid->elem2d[ie1].nodes[i]].y - y3d[0])
                                + fabs(grid->node[grid->elem2d[ie1].nodes[i]].z - z3d[0]);
                                
                                dist2 = fabs(grid->node[grid->elem2d[ie1].nodes[i]].x - x3d[1])
                                + fabs(grid->node[grid->elem2d[ie1].nodes[i]].y - y3d[1])
                                + fabs(grid->node[grid->elem2d[ie1].nodes[i]].z - z3d[1]);
                                
                                dist3 = fabs(grid->node[grid->elem2d[ie1].nodes[i]].x - x3d[2])
                                + fabs(grid->node[grid->elem2d[ie1].nodes[i]].y - y3d[2])
                                + fabs(grid->node[grid->elem2d[ie1].nodes[i]].z - z3d[2]);
                                if (dist1 < 1e-6 || dist2 < 1e-6 || dist3 < 1e-6) overlay++;
                            }
                            
                            
                            
                            if (overlay == 3) { /* found match for all 3 nodes, couple faces */
                                cpl_elems[count] = ie2;
                                count++;
                            }
                        }
                    }
                }
            }else{
                /* For 3D-2D, coupled grid is 2D, look at 1D elements */
                for (ie2 = 0; ie2 < grid2->nelems1d; ie2++) {
                    for (ifc=0;ifc<ninterfaces;ifc++){
                        if(grid2->elem1d[ie2].string==interface_string[ifc]){
                            
                            overlay = 0;
                            
                            for (i=0; i<2; i++) {
                                x3d[i] = grid2->node[grid2->elem1d[ie2].nodes[i]].x;
                                y3d[i] = grid2->node[grid2->elem1d[ie2].nodes[i]].y;
                            }
                            
                            dist1 = fabs(grid->node[grid->elem2d[ie1].nodes[0]].x - x3d[0]) + fabs(grid->node[grid->elem2d[ie1].nodes[0]].y - y3d[0]);
                            dist2 = fabs(grid->node[grid->elem2d[ie1].nodes[0]].x - x3d[1]) + fabs(grid->node[grid->elem2d[ie1].nodes[0]].y - y3d[1]);
                            if (dist1 < 1e-6 || dist2 < 1e-6 ) overlay++;
                            
                            dist1 = fabs(grid->node[grid->elem2d[ie1].nodes[1]].x - x3d[0]) + fabs(grid->node[grid->elem1d[ie1].nodes[1]].y - y3d[0]);
                            dist2 = fabs(grid->node[grid->elem2d[ie1].nodes[1]].x - x3d[1]) + fabs(grid->node[grid->elem1d[ie1].nodes[1]].y - y3d[1]);
                            if (dist1 < 1e-6 || dist2 < 1e-6 ) overlay++;
                            
                            dist1 = fabs(grid->node[grid->elem2d[ie1].nodes[2]].x - x3d[0]) + fabs(grid->node[grid->elem1d[ie1].nodes[2]].y - y3d[0]);
                            dist2 = fabs(grid->node[grid->elem2d[ie1].nodes[2]].x - x3d[1]) + fabs(grid->node[grid->elem1d[ie1].nodes[2]].y - y3d[1]);
                            if (dist1 < 1e-6 || dist2 < 1e-6 ) overlay++;
                            
                            if (overlay == 2) { // two of the 3D faces nodes match different edge nodes, add
                                cpl_elems[count] = ie2;
                                count++;
                            }
                        }
                    }
                }
            }
            /* 3D faces should only couple to 1 1D element or 2D face, but maybe more than 1 supermodel is coupled? */
            if(grid->elem2d[ie1].flux_elem !=NULL ){
                grid->elem2d[ie1].flux_elem = (FLUX_ELEM *) tl_realloc(sizeof(FLUX_ELEM), grid->elem2d[ie1].flux_elem_tot, (grid->elem2d[ie1].flux_elem_tot+count), grid->elem2d[ie1].flux_elem);
            }else{
                grid->elem2d[ie1].flux_elem = (FLUX_ELEM *) tl_alloc(sizeof(FLUX_ELEM), (grid->elem2d[ie1].flux_elem_tot+count));
            }
            for(i=0;i<count;i++){
                grid->elem2d[ie1].flux_elem[grid->elem2d[ie1].flux_elem_tot+count].elem = cpl_elems[i];
                grid->elem2d[ie1].flux_elem[grid->elem2d[ie1].flux_elem_tot+count].sub  = sub;
                grid->elem2d[ie1].flux_elem[grid->elem2d[ie1].flux_elem_tot+count].super= super;
                printf("2d elem %d sup %d sub %d elem %d \n", ie1,super,sub,cpl_elems[i]);
            }
            grid->elem2d[ie1].flux_elem_tot = count;
        }
    }
    
    
    
    
    
    /* Tidy up */
    if(grid->nelems3d > 0)  free_elem3d_list(elem3d_list);
    free_elem2d_list(elem2d_list);
    free_node_list(node_list);
    free_node_list(node_list_tmp);
    free_node_list(ghost_nodes);
    interface_string = (int *)tl_free(sizeof(int), ninterfaces, interface_string);
    if(grid2->ndim==3){
        cpl_elems = (int *) tl_free(sizeof(int),grid2->nelems2d, cpl_elems);
    }else{
        cpl_elems = (int *) tl_free(sizeof(int),grid2->nelems1d, cpl_elems);
    }
    sgrid_free(grid2);
    return;
}




