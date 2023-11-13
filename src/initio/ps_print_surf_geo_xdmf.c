/* writes the surface Geometry information to HDF5 file */
/* NOTE: actually we are printing out all 2D elements */

#include "global_header.h"
#ifdef _ADH_HDF5
#include "adh_hdf5.h"
#endif
void ps_print_surf_geo_xdmf(
#ifdef _ADH_HDF5
    SMODEL *mod /* The container group for mesh */
#endif
    )
{
#ifdef _ADH_HDF5
    hid_t dataspace, dataset, aid, attr;
    herr_t status;
    hsize_t dims[2];
    char string[MAXLINE];
    int i;
    int ie, ie2d, ii, ielem, jj;
    int nelem;
    int num_ghosts;
    int nodes_per_elem;
    int *elem_conn = NULL;        /* element connectivity table; has dimensions (nodes_per_elem) x (nelem) */
    int *node_map_data = NULL;    /* used to pack node_map structure data */
    int *elem3d_material;         /* temporary vector to store material ids for printing */
    double *node_coords = NULL;
    int *ghost_nodes = NULL;
    int compression;              /* turn on/off compression */
    int *surf_flag = NULL;
    int icount;
    int index;
    int my_nelems=0, total_nelems=0;
    int npes, myid, owner;
    int my_nnode, nnode, total_nnodes=0;

    assert(mod->flag.SW3_FLOW
#ifdef _ADH_GROUNDWATER
            || mod->flag.GW_FLOW
#endif
            );
    SGRID *grid = mod->grid;
    HDF5 *hdf5 = &(mod->hdf5);
    npes = grid->smpi->npes;
    myid = grid->smpi->myid;
    /* Note that current_mesh must be set elsewhere (namely in ps_print_geo_xdmf).
     * So, that routine MUST be called before this one
     */

    /* compression = print_control.hdf5_compression; */
    /* no compression yet LP*/
    compression=0;
    /* Adding 2D grid in case of SW3D models */
    /**********************************************************************************************************************************************************************************************/
    /**********************************************************************************************************************************************************************************************/
    /* Writing 2D Grid */
    //hdf5->current_mesh++;
    my_nnode = grid->my_nnodes_sur;
    nnode = grid->nnodes_sur;
    nelem = grid->nelems2d_sur;
    nodes_per_elem=3;

#ifdef _MESSG
    total_nnodes = messg_isum(my_nnode, grid->smpi->ADH_COMM);
#else
    total_nnodes = grid->nnodes_sur;
#endif

    // now write element and conectivities
    int nodeID_3d;
    elem_conn = (int *) tl_alloc(sizeof(int), nelem * nodes_per_elem);
    my_nelems=0;
    for (ie = 0; ie < grid->nelems2d_sur; ie++) {
        owner = npes;
        ie2d = grid->elem2d_sur[ie]; // 2d element id in the complete list of 2d elements
        for (jj=0; jj<nodes_per_elem; jj++){
            nodeID_3d = grid->elem2d[ie2d].nodes[jj];
            elem_conn[ie * nodes_per_elem + jj] = grid->nodeID_3d_to_2d_sur[nodeID_3d];
#ifdef _MESSG
            owner = MIN(grid->node[nodeID_3d].resident_pe, owner);
#else
            owner = myid;
#endif
        }
        if (owner == myid) my_nelems++;
    }
    total_nelems = 0;
#ifdef _MESSG
    total_nelems = messg_isum(my_nelems, grid->smpi->ADH_COMM);
#else
    total_nelems = my_nelems;
#endif
    /* Fill in nodal coordinate */
    node_coords = (double *) tl_alloc(sizeof(double), 3 * nnode);
#ifdef _ADH_GROUNDWATER
    if (mod->flag.SW3_FLOW){
#endif
        for (ii = 0; ii < nnode; ii++){
            nodeID_3d = grid->nodeID_2d_to_3d_bed[ii];
            node_coords[3 * ii]     = grid->node[nodeID_3d].x;
            node_coords[3 * ii + 1] = grid->node[nodeID_3d].y;
            node_coords[3 * ii + 2] = grid->node[nodeID_3d].z;
        }
#ifdef _ADH_GROUNDWATER
    }
    else if (mod->flag.GW_FLOW){
        for (ii = 0; ii < nnode; ii++){
            nodeID_3d = grid->nodeID_2d_to_3d_sur[ii]; /* Note: GW needs sur! */
            node_coords[3 * ii]     = grid->node[nodeID_3d].x;
            node_coords[3 * ii + 1] = grid->node[nodeID_3d].y;
            node_coords[3 * ii + 2] = grid->node[nodeID_3d].z;
        }
    }
#endif

    /* Write out the element connectivity */
    dims[0] = nelem;
    dims[1] = nodes_per_elem;
    snprintf(string, MAXLINE, "%d", hdf5->current_mesh);
    xdmf_write_dataset(hdf5->group_surf_elements, string, H5T_NATIVE_INT, compression, 2, dims, elem_conn);
    /* Add the global number of elements as an attribute to element connectivity data set */
    dims[0] = 1;
    dataset = H5Dopen(hdf5->group_surf_elements, string, H5P_DEFAULT);
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
    xdmf_write_dataset(hdf5->group_surf_nodes, string, H5T_NATIVE_DOUBLE, compression, 2, dims, node_coords);

    /* Write out the ghost nodes */
    num_ghosts = nnode - my_nnode;
    if (num_ghosts > 0)
      {
        ghost_nodes = (int *) tl_alloc(sizeof(int), num_ghosts);
        for (ii = 0; ii < num_ghosts; ii++)
          {
            ghost_nodes[ii] = my_nnode + ii;
          }
        dims[0] = num_ghosts;
        snprintf(string, MAXLINE, "%d", hdf5->current_mesh);
        xdmf_write_dataset(hdf5->group_surf_ghost_nodes, string, H5T_NATIVE_INT, compression, 1, dims, ghost_nodes);
        ghost_nodes = (int *) tl_free(sizeof(int), num_ghosts, ghost_nodes);
      }

    /* Add the global and processor's number of nodes as attributes to the nodes data set */
    dims[0] = 1;
    dataset = H5Dopen(hdf5->group_surf_nodes, string, H5P_DEFAULT);
    aid = H5Screate_simple(1, dims, NULL);
    attr = H5Acreate(dataset, "global_nnode", H5T_NATIVE_INT, aid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &total_nnodes);
    attr = H5Acreate(dataset, "my_nnode", H5T_NATIVE_INT, aid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &my_nnode);
    status = H5Sclose(aid);
    status = H5Aclose(attr);
    status = H5Dclose(dataset);

    /* Write out the original node numbers */
    //dims[0] = nnode;
    //int *orig_nd_number = (int *) tl_alloc(sizeof(int), nnode);
    //for (ii=0; ii<nnode; ii++){
    //    nodeID_3d = grid->nodeID_2d_to_3d_sur[ii];
    //    orig_nd_number[ii] = grid->node[nodeID_3d].global_surf_id;
    //    //orig_nd_number[ii] = grid->node[ii].original_id;
    //}
    //xdmf_write_dataset(hdf5->group_surf_orig_nd, string, H5T_NATIVE_INT, compression, 1, dims, orig_nd_number);
    //orig_nd_number = (int*) tl_free(sizeof(int), nnode, orig_nd_number);

    /* Write out original material ID */
    //dims[0] = nelem;
    //snprintf(string, MAXLINE, "%d", hdf5->current_mesh);
    //int *value_array = (int *) tl_alloc(sizeof(int), nelem);
    //for (ie = 0; ie < nelem; ie++) {
    //    ie2d = grid->elem2d_sur[ie]; // 2d element id in the complete list of 2d elements
    //    value_array[ie] = grid->elem2d[ie2d].mat + 1;
    //}
    //xdmf_write_dataset(hdf5->group_surf_mat_id, string, H5T_NATIVE_INT, compression, 1, dims, value_array);
    //value_array = (int *) tl_free(sizeof(int), nelem, value_array);

    /* Pack node_map data and write it out */
    node_map_data = (int *) tl_alloc(sizeof(int), 2 * nnode);
    for (ii = 0; ii < nnode; ii++) {
        nodeID_3d = grid->nodeID_2d_to_3d_sur[ii];
#ifdef _MESSG
        node_map_data[2 * ii] = grid->node[nodeID_3d].resident_pe; 
        node_map_data[2 * ii + 1] = grid->node[nodeID_3d].global_surf_id;
#else
        node_map_data[2 * ii] = myid;
        node_map_data[2 * ii + 1] = grid->node[nodeID_3d].global_surf_id;
#endif
    }

    dims[0] = nnode;
    dims[1] = 2;
    snprintf(string, MAXLINE, "%d", hdf5->current_mesh);
    xdmf_write_dataset(hdf5->group_surf_node_map, string, H5T_NATIVE_INT, compression, 2, dims, node_map_data);
    node_map_data = (int *) tl_free(sizeof(int), 2 * nnode, node_map_data);
    elem_conn = (int *) tl_free(sizeof(int), nelem * nodes_per_elem, elem_conn);
    node_coords = (double *) tl_free(sizeof(double), 3 * nnode, node_coords);


#endif
  return;
}

