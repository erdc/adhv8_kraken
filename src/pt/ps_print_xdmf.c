/* write the solution to a XDMF file */

#include "global_header.h"

void ps_print_xdmf(SGRID *grid, SVECT *vel, double *dpl, double time) {
    
#ifdef _ADH_HDF5
    double *index;                /* the index for scaling the error */
    double c_time = time;           /* converted time */
    int i, ie, ibl;               /* loop counter and end of loop */
    double conv1, convr, conv3, conv4;
    char string[MAXLINE];
    char buff1[MAXLINE] = "";
    char file_backup[MAXLINE];
    herr_t status;
    xmlNodePtr node1 = NULL, node2 = NULL;
    FILE *fp;
    
    HDF5 *hdf5 = &(grid->hdf5);
    int nnode = grid->nnodes;
    
    /* Open the HDF file and groups */
    hdf5->h5_file_ptr    = H5Fopen(hdf5->h5_filename, H5F_ACC_RDWR, H5P_DEFAULT);
    hdf5->group_mesh     = H5Gopen(hdf5->h5_file_ptr, "Geometry", H5P_DEFAULT);
    hdf5->group_data     = H5Gopen(hdf5->h5_file_ptr, "Data", H5P_DEFAULT);
    hdf5->group_elements = H5Gopen(hdf5->group_mesh, "elements", H5P_DEFAULT);
    hdf5->group_nodes    = H5Gopen(hdf5->group_mesh, "nodes", H5P_DEFAULT);
    
    //printf("Current time = %f",time);
    
    /* Set XDMF Element Type */
    char elem_type[20];
    int nelem = UNSET_INT;
    if (grid->ndim == 3){
        nelem = grid->nelems3d;
        sprintf(elem_type, "Tetrahedron");
    } else {
        nelem = grid->nelems2d;
        sprintf(elem_type, "Triangle");
    }
    
    /* update XML tree <Grid/> */
    hdf5->current_grid_node = xmlNewChild(hdf5->collection_node, NULL, BAD_CAST "Grid", NULL);
    xmlNewProp(hdf5->current_grid_node, BAD_CAST "GridType", BAD_CAST "Uniform");
    
    /* add time stamp <Time/> */
    node1 = xmlNewNode(NULL, BAD_CAST "Time");
    sprintf(buff1, "%f", time);
    xmlNewProp(node1, BAD_CAST "Value", BAD_CAST buff1);
    xmlAddChild(hdf5->current_grid_node, node1);
    
    /* TOPOLOGY <Topology/> data */
    node1 = xmlNewNode(NULL, BAD_CAST "Topology");
    sprintf(buff1, "%d", nelem);
    xmlNewProp(node1, BAD_CAST "NumberOfElements", BAD_CAST buff1);
    xmlNewProp(node1, BAD_CAST "Type", BAD_CAST elem_type);
    sprintf(buff1, "%s:/Geometry/elements/%d", hdf5->h5_filename, hdf5->current_mesh);
    node2 = xmlNewChild(node1, NULL, BAD_CAST "DataItem", BAD_CAST buff1);
    xmlNewProp(node2, BAD_CAST "DataType", BAD_CAST "Int");
    if (grid->ndim == 3) {
        sprintf(buff1, "%d 4", grid->nelems3d);
    } else {
        sprintf(buff1, "%d 3", grid->nelems2d);
    }
    xmlNewProp(node2, BAD_CAST "Dimensions", BAD_CAST buff1);
    xmlNewProp(node2, BAD_CAST "Format", BAD_CAST "HDF");
    xmlAddChild(hdf5->current_grid_node, node1);
    
    
    /* GEOMETRY <Geomtery/> data */
    node1 = xmlNewNode(NULL, BAD_CAST "Geometry");
    xmlNewProp(node1, BAD_CAST "Type", BAD_CAST "XYZ");
    sprintf(buff1, "%s:/Geometry/nodes/%d", hdf5->h5_filename, hdf5->current_mesh);
    node2 = xmlNewChild(node1, NULL, BAD_CAST "DataItem", BAD_CAST buff1);
    xmlNewProp(node2, BAD_CAST "DataType", BAD_CAST "Float");
    sprintf(buff1, "%d 3", nnode);
    xmlNewProp(node2, BAD_CAST "Dimensions", BAD_CAST buff1);
    xmlNewProp(node2, BAD_CAST "Format", BAD_CAST "HDF");
    xmlNewProp(node2, BAD_CAST "Precision", BAD_CAST "8");
    xmlAddChild(hdf5->current_grid_node, node1);
    
    /* add Information element to supply grid->my_nnodes (number of nodes owned by this process) */
    node2 = xmlNewChild(node1, NULL, BAD_CAST "Information", NULL);
    xmlNewProp(node2, BAD_CAST "Name", BAD_CAST "my_nnode");
    sprintf(buff1, "%d", nnode);
    xmlNewProp(node2, BAD_CAST "Value", BAD_CAST buff1);
    
    xmlAddChild(hdf5->current_grid_node, node1);
    
    /*************************************************************************************************************************************************/
    
    if (vel==NULL) fprintf(stderr, "Memory Error %s:%d\n", __FILE__,__LINE__);
    ps_print_ts_xdmf(hdf5,"Velocity", vel, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, VECTOR3D_DATA);
    if (dpl != NULL) ps_print_ts_xdmf(hdf5,"Displacement", dpl, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
    
    /*************************************************************************************************************************************************/
    /* close HDF groups and file */
    status = H5Gclose(hdf5->group_nodes);
    status = H5Gclose(hdf5->group_elements);
    status = H5Gclose(hdf5->group_data);
    status = H5Gclose(hdf5->group_mesh);
    status = H5Fclose(hdf5->h5_file_ptr);
    
    /* If the xmf file exists, then move rename it before writing the current file.
     In this way, if the program terminates prematurely, we still have a copy of the
     state of of the output up until this point */
    
    fp = fopen(hdf5->xmf_filename, "r");
    
    if (fp != NULL) {
        fclose(fp);
        snprintf(string, MAXLINE - 1, "%s_back", hdf5->xmf_filename);
        rename(hdf5->xmf_filename, string);
    }
    /* Write the xmf file */
    
    xmlSaveFormatFileEnc(hdf5->xmf_filename, hdf5->doc, "UTF-8", 1);
    
    /*************************************************************************************************************************************************/
    
    /* increments the printstep */
    hdf5->iprint_step++;
    
    /* fflush(NULL); */
#endif
    return;
}



