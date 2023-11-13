#include "global_header.h"

static double pi = 3.14159265359;

//****************************************************************************//
//****************************************************************************//
// cjt :: find all nodes on the slice x,y = constant, sort and return them
int tl_find_nodes_on_sliceX(SGRID *grid, int nd0, int nd1, int **GnodeIDs_on_slice) {
    
    int i,j,inode = UNSET_INT;
    double x, x0, x1, y, y0, y1, residual = 0.;
    
    int node_count = 0; // total nodes found on slice

#ifdef _MESSG
    int root0=-1,root1=-1, ierr;
    for(i=0;i<grid->my_nnodes;i++){
      if((grid->node[i].gid == nd0) && (grid->node[i].resident_pe == grid->smpi->myid)) {
          root0 = grid->smpi->myid;
          x0 = grid->node[i].x;
          y0 = grid->node[i].y;
        }

      if((grid->node[i].gid == nd1) && (grid->node[i].resident_pe == grid->smpi->myid)) {
          root1 = grid->smpi->myid;
          x1 = grid->node[i].x;
          y1 = grid->node[i].y;
       }

    }
    root0 =  messg_imax(root0, grid->smpi->ADH_COMM);
    root1 =  messg_imax(root1, grid->smpi->ADH_COMM);
    ierr = MPI_Bcast(&x0, 1, MPI_DOUBLE, root0, grid->smpi->ADH_COMM);
    ierr = MPI_Bcast(&y0, 1, MPI_DOUBLE, root0, grid->smpi->ADH_COMM);
    ierr = MPI_Bcast(&x1, 1, MPI_DOUBLE, root1, grid->smpi->ADH_COMM);
    ierr = MPI_Bcast(&y1, 1, MPI_DOUBLE, root1, grid->smpi->ADH_COMM);

    i=UNSET_INT;
#else
    // lagrange linear interpolate {y = y0 * ((x-x1)/(x0-x1)) + y1 * ((x-x0)/(x1-x0))}
    x0 = grid->node[nd0].x; x1 = grid->node[nd1].x;
    y0 = grid->node[nd0].y; y1 = grid->node[nd1].y;
#endif   

    for (inode=0; inode<grid->my_nnodes; inode++) {
           
        // is line equation satisfied?
        x = grid->node[inode].x; y = grid->node[inode].y;
        if (fabs(x0 - x1)>1e-6) {
            residual = y - (y0 * ((x-x1)/(x0-x1)) + y1 * ((x-x0)/(x1-x0)));
            if (fabs(residual) < 1e-6) {
                node_count++;
            }
        } else {
            if (fabs(x - x0)<1e-6) {
                node_count++;
            }
        }
    }

  if ((node_count == 0) && (grid->smpi->npes==1)) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> There are no nodes on this 2d slice.");
    }
  if(node_count>0){
    (*GnodeIDs_on_slice) = (int *) tl_alloc(sizeof(int), node_count);
    int *nodeIDs_on_slice = (*GnodeIDs_on_slice);
    
    int count = 0;
    for (inode=0; inode<grid->my_nnodes; inode++) {

        // is line equation satisfied?
        x = grid->node[inode].x; y = grid->node[inode].y;
        if (fabs(x0 - x1)>1e-6) {
            residual = y - (y0 * ((x-x1)/(x0-x1)) + y1 * ((x-x0)/(x1-x0)));
            if (fabs(residual) < 1e-6) {
                nodeIDs_on_slice[count] = inode;
                count++;
            }
        } else {
            if (fabs(x - x0)<1e-6) {
                nodeIDs_on_slice[count] = inode;
                count++;
            }
        }
    }
    
    // Now order the nodes from the varying coordinate size using bubble sort
    int temp = 0;
    for(i=0; i<node_count; i++) {
        for(j=0; j<node_count-1; j++) {
            if (grid->node[nodeIDs_on_slice[j]].y > grid->node[nodeIDs_on_slice[j+1]].y) {
                temp = nodeIDs_on_slice[j+1];
                nodeIDs_on_slice[j+1] = nodeIDs_on_slice[j];
                nodeIDs_on_slice[j] = temp;
            }
        }
    }
  }

   return node_count;
}

//****************************************************************************//
//****************************************************************************//
// cjt :: find all nodes on the slice x,y = constant, sort and return them
int tl_find_nodes_on_sliceY(SGRID *grid, int nd0, int nd1, int **GnodeIDs_on_slice) {
    
    int i,j,inode = UNSET_INT;
    double x, x0, x1, y, y0, y1, residual = 0.;

#ifdef _MESSG
    int root0=-1,root1=-1, ierr;
    for(i=0;i<grid->my_nnodes;i++){
      if((grid->node[i].gid == nd0) && (grid->node[i].resident_pe == grid->smpi->myid)) {
          root0 = grid->smpi->myid;
          x0 = grid->node[i].x;
          y0 = grid->node[i].y;
        }

      if((grid->node[i].gid == nd1) && (grid->node[i].resident_pe == grid->smpi->myid)) {
          root1 = grid->smpi->myid;
          x1 = grid->node[i].x;
          y1 = grid->node[i].y;
       }

    }
    root0 =  messg_imax(root0, grid->smpi->ADH_COMM);
    root1 =  messg_imax(root1, grid->smpi->ADH_COMM);
    ierr = MPI_Bcast(&x0, 1, MPI_DOUBLE, root0, grid->smpi->ADH_COMM);
    ierr = MPI_Bcast(&y0, 1, MPI_DOUBLE, root0, grid->smpi->ADH_COMM);
    ierr = MPI_Bcast(&x1, 1, MPI_DOUBLE, root1, grid->smpi->ADH_COMM);
    ierr = MPI_Bcast(&y1, 1, MPI_DOUBLE, root1, grid->smpi->ADH_COMM);
    i=UNSET_INT;
#else
    // lagrange linear interpolate {y = y0 * ((x-x1)/(x0-x1)) + y1 * ((x-x0)/(x1-x0))}
    x0 = grid->node[nd0].x; x1 = grid->node[nd1].x;
    y0 = grid->node[nd0].y; y1 = grid->node[nd1].y;
#endif

    int node_count = 0; // total nodes found on slice
    for (inode=0; inode<grid->my_nnodes; inode++) {
        
       
        // is line equation satisfied?
        x = grid->node[inode].x; y = grid->node[inode].y;
        if (fabs(x0 - x1)>1e-6) {
            residual = y - (y0 * ((x-x1)/(x0-x1)) + y1 * ((x-x0)/(x1-x0)));
            if (fabs(residual) < 1e-6) {
                node_count++;
            }
        } else {
            if (fabs(x - x0)<1e-6) {
                node_count++;
            }
        }
    }
    
    if ((node_count == 0) &&(grid->smpi->npes==1)) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> There are no nodes on this 2d slice.");
    }
    if (node_count>0){ 
    (*GnodeIDs_on_slice) = (int *) tl_alloc(sizeof(int), node_count);
    int *nodeIDs_on_slice = (*GnodeIDs_on_slice);
    
    int count = 0;
    for (inode=0; inode<grid->my_nnodes; inode++) {
         
        // is line equation satisfied?
        x = grid->node[inode].x; y = grid->node[inode].y;
        if (fabs(x0 - x1)>1e-6) {
            residual = y - (y0 * ((x-x1)/(x0-x1)) + y1 * ((x-x0)/(x1-x0)));
            if (fabs(residual) < 1e-6) {
                nodeIDs_on_slice[count] = inode;
                count++;
            }
        } else {
            if (fabs(x - x0)<1e-6) {
                nodeIDs_on_slice[count] = inode;
                count++;
            }
        }
    }
    
    // Now order the nodes from the varying coordinate size using bubble sort
    int temp = 0;
    for(i=0; i<node_count; i++) {
        for(j=0; j<node_count-1; j++) {
            if (grid->node[nodeIDs_on_slice[j]].x > grid->node[nodeIDs_on_slice[j+1]].x) {
                temp = nodeIDs_on_slice[j+1];
                nodeIDs_on_slice[j+1] = nodeIDs_on_slice[j];
                nodeIDs_on_slice[j] = temp;
            }
        }
    }
    }
    else{printf("MYID %d has zero nodes on the slice \n", grid->smpi->myid);}
    
    return node_count;
}

//****************************************************************************//
//****************************************************************************//
// cjt :: find all nodes on the slice x,y = constant, sort and return them
int tl_find_surface_nodes_on_sliceX(SGRID *grid, int nd0, int nd1, int **GnodeIDs_on_slice) {
    
    int i,j,id,inode = UNSET_INT;
    double x, x0, x1, y, y0, y1, residual = 0.;
    ID_LIST_ITEM *ptr;

#ifdef _MESSG
    int root0=-1,root1=-1, ierr;

    for(i=0;i<grid->nnodes_sur;i++){
      if((grid->node[grid->nodeID_2d_to_3d_sur[i]].gid == nd0) && (grid->node[grid->nodeID_2d_to_3d_sur[i]].resident_pe == grid->smpi->myid)) {
          root0 = grid->smpi->myid;
          x0 = grid->node[grid->nodeID_2d_to_3d_sur[i]].x;
          y0 = grid->node[grid->nodeID_2d_to_3d_sur[i]].y;
        }
       
      if((grid->node[grid->nodeID_2d_to_3d_sur[i]].gid == nd1) && (grid->node[grid->nodeID_2d_to_3d_sur[i]].resident_pe == grid->smpi->myid)) {
          root1 = grid->smpi->myid;
          x1 = grid->node[grid->nodeID_2d_to_3d_sur[i]].x;
          y1 = grid->node[grid->nodeID_2d_to_3d_sur[i]].y;
       }
      
    }
    root0 =  messg_imax(root0, grid->smpi->ADH_COMM);
    root1 =  messg_imax(root1, grid->smpi->ADH_COMM);
    ierr = MPI_Bcast(&x0, 1, MPI_DOUBLE, root0, grid->smpi->ADH_COMM);
    ierr = MPI_Bcast(&y0, 1, MPI_DOUBLE, root0, grid->smpi->ADH_COMM);
    ierr = MPI_Bcast(&x1, 1, MPI_DOUBLE, root1, grid->smpi->ADH_COMM);
    ierr = MPI_Bcast(&y1, 1, MPI_DOUBLE, root1, grid->smpi->ADH_COMM);
    i=UNSET_INT;
#else
    // lagrange linear interpolate {y = y0 * ((x-x1)/(x0-x1)) + y1 * ((x-x0)/(x1-x0))}
    x0 = grid->node[nd0].x; x1 = grid->node[nd1].x;
    y0 = grid->node[nd0].y; y1 = grid->node[nd1].y;
#endif 

    int node_count = 0; // total nodes found on slice
    for (inode=0; inode<grid->nnodes_sur; inode++) {
        ptr = grid->vertical_list[inode];
        id = ptr->id;
        
        // lagrange linear interpolate {y = y0 * ((x-x1)/(x0-x1)) + y1 * ((x-x0)/(x1-x0))}       
        // is line equation satisfied?
#ifdef _MESSG
      if(grid->node[id].resident_pe == grid->smpi->myid){
#endif
        x = grid->node[id].x; y = grid->node[id].y;
        if (fabs(x0 - x1)>1e-6) {
            residual = y - (y0 * ((x-x1)/(x0-x1)) + y1 * ((x-x0)/(x1-x0)));
            if (fabs(residual) < 1e-6) {
                node_count++;
            }
        } else {
            if (fabs(x - x0)<1e-6) {
                node_count++;
            }
        }
#ifdef _MESSG
      }
#endif
    }
    
    if ((node_count == 0) &&(grid->smpi->npes==1)) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> There are no nodes on this 2d slice.");
    }

    if(node_count != 0)  (*GnodeIDs_on_slice) = (int *) tl_alloc(sizeof(int), node_count);
    int *nodeIDs_on_slice = (*GnodeIDs_on_slice);

    if(node_count !=0 ){
    int count = 0;
    for (inode=0; inode<grid->nnodes_sur; inode++) {
        ptr = grid->vertical_list[inode];
        id = ptr->id;
        
        // lagrange linear interpolate {y = y0 * ((x-x1)/(x0-x1)) + y1 * ((x-x0)/(x1-x0))}
        // is line equation satisfied?
#ifdef _MESSG
      if(grid->node[id].resident_pe == grid->smpi->myid){
#endif
        x = grid->node[id].x; y = grid->node[id].y;
        if (fabs(x0 - x1)>1e-6) {
            residual = y - (y0 * ((x-x1)/(x0-x1)) + y1 * ((x-x0)/(x1-x0)));
            if (fabs(residual) < 1e-6) {
                nodeIDs_on_slice[count] = id;
                count++;
            }
        } else {
            if (fabs(x - x0)<1e-6) {
                nodeIDs_on_slice[count] = id;
                count++;
            }
        }
#ifdef _MESSG
      }
#endif
    }
    
    
    // Now order the nodes from the varying coordinate size using bubble sort
    int temp = UNSET_INT;
    for(i=0; i<node_count; i++) {
        for(j=0; j<node_count-1; j++) {
            if (grid->node[nodeIDs_on_slice[j]].y > grid->node[nodeIDs_on_slice[j+1]].y) {
                temp = nodeIDs_on_slice[j+1];
                nodeIDs_on_slice[j+1] = nodeIDs_on_slice[j];
                nodeIDs_on_slice[j] = temp;
            }
        }
    }
    }    
    /*
    printf("node_count: %d\n",node_count);
    for(i=0; i<node_count; i++) {
        printf("id: %d nodeIDs_on_slice: %d\n",i+1,nodeIDs_on_slice[i]+1);
    }
    */
    
    
    return node_count;
}

//****************************************************************************//
//****************************************************************************//
// cjt :: find all nodes on the slice x,y = constant, sort and return them
int tl_find_surface_nodes_on_sliceY(SGRID *grid, int nd0, int nd1, int **GnodeIDs_on_slice) {
    
    int i,j,id,inode = UNSET_INT;
    double x, x0, x1, y, y0, y1, residual = 0.;
    ID_LIST_ITEM *ptr;
#ifdef _MESSG
    int root0 = -1,root1 = -1, ierr;
    for(i=0;i<grid->nnodes_sur;i++){
      if((grid->node[grid->nodeID_2d_to_3d_sur[i]].gid == nd0) && (grid->node[grid->nodeID_2d_to_3d_sur[i]].resident_pe == grid->smpi->myid)) {
          root0 = grid->smpi->myid;
          x0 = grid->node[grid->nodeID_2d_to_3d_sur[i]].x;
          y0 = grid->node[grid->nodeID_2d_to_3d_sur[i]].y;
        }

      if((grid->node[grid->nodeID_2d_to_3d_sur[i]].gid == nd1) && (grid->node[grid->nodeID_2d_to_3d_sur[i]].resident_pe == grid->smpi->myid)) {
          root1 = grid->smpi->myid;
          x1 = grid->node[grid->nodeID_2d_to_3d_sur[i]].x;
          y1 = grid->node[grid->nodeID_2d_to_3d_sur[i]].y;
       }

    }
    root0 =  messg_imax(root0, grid->smpi->ADH_COMM);
    root1 =  messg_imax(root1, grid->smpi->ADH_COMM);
    ierr = MPI_Bcast(&x0, 1, MPI_DOUBLE, root0, grid->smpi->ADH_COMM);
    ierr = MPI_Bcast(&y0, 1, MPI_DOUBLE, root0, grid->smpi->ADH_COMM);
    ierr = MPI_Bcast(&x1, 1, MPI_DOUBLE, root1, grid->smpi->ADH_COMM);
    ierr = MPI_Bcast(&y1, 1, MPI_DOUBLE, root1, grid->smpi->ADH_COMM);
    i=UNSET_INT;
#else
    // lagrange linear interpolate {y = y0 * ((x-x1)/(x0-x1)) + y1 * ((x-x0)/(x1-x0))}
    x0 = grid->node[nd0].x; x1 = grid->node[nd1].x;
    y0 = grid->node[nd0].y; y1 = grid->node[nd1].y;
#endif
    int node_count = 0; // total nodes found on slice
    for (inode=0; inode<grid->nnodes_sur; inode++) {
        ptr = grid->vertical_list[inode];
        id = ptr->id;
        
        // lagrange linear interpolate {y = y0 * ((x-x1)/(x0-x1)) + y1 * ((x-x0)/(x1-x0))}
        // is line equation satisfied?
#ifdef _MESSG
      if(grid->node[id].resident_pe == grid->smpi->myid){
#endif
        x = grid->node[id].x; y = grid->node[id].y;
        if (fabs(x0 - x1)>1e-6) {
            residual = y - (y0 * ((x-x1)/(x0-x1)) + y1 * ((x-x0)/(x1-x0)));
            if (fabs(residual) < 1e-6) {
                node_count++;
            }
        } else {
            if (fabs(x - x0)<1e-6) {
                node_count++;
            }
        }
#ifdef _MESSG
      }
#endif
    }
    
    if (node_count == 0&&(grid->smpi->npes==1)) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> There are no nodes on this 2d slice.");
    }
    
    if(node_count != 0 ) (*GnodeIDs_on_slice) = (int *) tl_alloc(sizeof(int), node_count);
    int *nodeIDs_on_slice = (*GnodeIDs_on_slice);
    
    if(node_count != 0 ){
    int count = 0;
    for (inode=0; inode<grid->nnodes_sur; inode++) {
        ptr = grid->vertical_list[inode];
        id = ptr->id;
        
        // lagrange linear interpolate {y = y0 * ((x-x1)/(x0-x1)) + y1 * ((x-x0)/(x1-x0))}
        // is line equation satisfied?
#ifdef _MESSG
      if(grid->node[id].resident_pe == grid->smpi->myid){
#endif
        x = grid->node[id].x; y = grid->node[id].y;
        if (fabs(x0 - x1)>1e-6) {
            residual = y - (y0 * ((x-x1)/(x0-x1)) + y1 * ((x-x0)/(x1-x0)));
            if (fabs(residual) < 1e-6) {
                nodeIDs_on_slice[count] = id;
                count++;
            }
        } else {
            if (fabs(x - x0)<1e-6) {
                nodeIDs_on_slice[count] = id;
                count++;
            }
        }
#ifdef _MESSG
      }
#endif
    }
    
    // Now order the nodes from the varying coordinate size using bubble sort
    int temp = 0;
    for(i=0; i<node_count; i++) {
        for(j=0; j<node_count-1; j++) {
            if (grid->node[nodeIDs_on_slice[j]].x > grid->node[nodeIDs_on_slice[j+1]].x) {
                temp = nodeIDs_on_slice[j+1];
                nodeIDs_on_slice[j+1] = nodeIDs_on_slice[j];
                nodeIDs_on_slice[j] = temp;
            }
        }
    }
    }
    return node_count;
}

                                                                
                                                                
                                                            
                                                                
                                        
                                                                

