/* writes the Geometry information to HDF5 file */
// CJT :: NOT SURE HOW THIS IS GOING OT HANDLE MIXED GRIDS!!!


#include "global_header.h"
#ifdef _ADH_HDF5
#include "adh_hdf5.h"
#endif

void ps_print_geo_xdmf(
#ifdef _ADH_HDF5
                       SGRID *grid
#endif
)
{
#ifdef _ADH_HDF5
    
    int ii,jj,ielem, owner;
    
    hid_t dataspace, dataset, aid, attr;
    herr_t status;
    hsize_t dims[2];
    char string[MAXLINE];
    
    HDF5 *hdf5 = &(grid->hdf5);
    int nnode = grid->nnodes;
    int nelem3d = grid->nelems3d;
    int nelem2d = grid->nelems2d;
    int npes = 1;
    int myid = 0;
    
    int compression = 0;
    hdf5->current_mesh++;
    
    int nelem = UNSET_INT, nodes_per_elem = UNSET_INT;
    if (grid->ndim == 3) {
        nelem = nelem3d;
        nodes_per_elem = 4;
        
    } else {
        nelem = nelem2d;
        nodes_per_elem = 3;
    }
    assert(nelem>0);
    assert(nodes_per_elem>0);
    
    // fill in connetivity
    int total_nelems = 0;
    int *elem_conn = (int *) malloc(sizeof(int) * nelem * nodes_per_elem);
    if (grid->ndim == 3) {
        for (ielem = 0; ielem < nelem; ielem++) {
            owner = npes;
            for (jj = 0; jj < 4; jj++) {
                elem_conn[ielem * nodes_per_elem + jj] = grid->elem3d[ielem].nodes[jj];
                owner = myid;
            }
            if (owner == myid){++total_nelems;}
        }
    } else {
        for (ielem = 0; ielem < nelem; ielem++) {
            owner = npes;
            for (jj = 0; jj < 3; jj++) {
                elem_conn[ielem * nodes_per_elem + jj] = grid->elem2d[ielem].nodes[jj];
                owner = myid;
            }
            if (owner == myid){++total_nelems;}
        }
    }
    assert(total_nelems>0);
    
    //printf("nelem: %d nodes_per_elem: %d total_nelems: %d\n",nelem,nodes_per_elem,total_nelems);
    
    /* Fill in nodal coordinate */
    double *node_coords = (double *) malloc(sizeof(double) * 3 * nnode);
    if (node_coords==NULL) {
        fprintf(stderr, "Could not malloc memory for nodes; %s:%d\n", __FILE__,__LINE__);
        exit(-1);
    }
    for (ii = 0; ii < nnode; ii++) {
        node_coords[3 * ii]     = grid->node[ii].x;
        node_coords[3 * ii + 1] = grid->node[ii].y;
        node_coords[3 * ii + 2] = grid->node[ii].z;
    }
    
    /* Write out the element connectivity */
    dims[0] = nelem;
    dims[1] = nodes_per_elem;
    snprintf(string, MAXLINE, "%d", hdf5->current_mesh);
    xdmf_write_dataset(hdf5->group_elements, string, H5T_NATIVE_INT, compression, 2, dims, elem_conn);
    
    /* Add the global number of elements as an attribute to element connectivity data set */
    dims[0] = 1;
    dataset = H5Dopen(hdf5->group_elements, string, H5P_DEFAULT);
    aid = H5Screate_simple(1, dims, NULL);
    attr = H5Acreate(dataset, "global_nelem", H5T_NATIVE_INT, aid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &total_nelems);
    status = H5Aclose(attr);
    status = H5Sclose(aid);
    status = H5Dclose(dataset);
    
    /* Write out the nodal coordinates */
    dims[0] = nnode;
    dims[1] = 3;
    snprintf(string, MAXLINE, "%d", hdf5->current_mesh);
    xdmf_write_dataset(hdf5->group_nodes, string, H5T_NATIVE_DOUBLE, compression, 2, dims, node_coords);
    

    free(elem_conn);
    free(node_coords);
#endif
    return;
}

