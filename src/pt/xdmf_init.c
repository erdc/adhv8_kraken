#define HDF5_INIT
#include "global_header.h"
/* ADH Version 2.0.0 6-04 */
/* Initialize use of XMDF/HDF5 */

#ifdef _ADH_HDF5
void xdmf_init(SGRID *grid, int npes, int myid, const char *proj_name)
{
    HDF5 *hdf5= &(grid->hdf5);
    hid_t aid, attr;
    hsize_t dims[1];
    FILE *fp;
    //xmlDtdPtr dtd = NULL;
    xmlNodePtr node1 = NULL;
    char file_backup[HDF5_NAME_LEN+6];
        
    hdf5->iprint_step = 0;
    
    /* Initialize time step counter and mesh counter */
    hdf5->current_mesh = -1;
    hdf5->current_ts = 0;
    
    /* Form filenames for the HDF5 file and the XMF file
     If we are running in parallel, add "_myid" to filenames */
    if (npes > 1) {
        snprintf(hdf5->xmf_filename, HDF5_NAME_LEN, "%s_grid_p%d.xmf", proj_name, myid);
        snprintf(hdf5->h5_filename, HDF5_NAME_LEN, "%s_grid_p%d.h5", proj_name, myid);
    } else {
        snprintf(hdf5->xmf_filename, HDF5_NAME_LEN, "%s_grid_pall1.xmf", proj_name);
        snprintf(hdf5->h5_filename, HDF5_NAME_LEN, "%s_grid_p0.h5", proj_name);
    }
    
    /* test that we can create the xmf file */
    fp = fopen(hdf5->xmf_filename, "w");
    if (fp == NULL){fprintf(stderr,"Unable to create xmf file; %s:%d", __FILE__, __LINE__);}
    fclose(fp);
    
    /* If the *h5 file exists, then move it to backup copy */
    fp = fopen(hdf5->h5_filename, "r");
    if (fp != NULL) {
        fclose(fp);
        snprintf(file_backup, HDF5_NAME_LEN + 5, "%s_back", hdf5->h5_filename);
        rename(hdf5->h5_filename, file_backup);
    }
    
    /* Open the HDF5 output file */
    hdf5->h5_file_ptr = H5Fcreate(hdf5->h5_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_ID(hdf5->h5_file_ptr);
    
    /* Add the number of processors as an attributes */
    dims[0] = 1;
    aid = H5Screate_simple(1, dims, NULL);
    HDF5_ID(aid);
    attr = H5Acreate(hdf5->h5_file_ptr, "npes", H5T_NATIVE_INT, aid, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_ID(attr);
    HDF5_ERR(H5Awrite(attr, H5T_NATIVE_INT, &npes));
    HDF5_ERR(H5Aclose(attr));
    HDF5_ERR(H5Sclose(aid));
    
    /* Create the container group for the meshes */
    hdf5->group_mesh = H5Gcreate(hdf5->h5_file_ptr, "/Geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_ID(hdf5->group_mesh);
    
    hdf5->group_nodes = H5Gcreate(hdf5->group_mesh, "nodes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_ID(hdf5->group_nodes);
    
    hdf5->group_elements = H5Gcreate(hdf5->group_mesh, "elements", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_ID(hdf5->group_elements);
    
    /* Create the container group for the data sets */
    hdf5->group_data = H5Gcreate(hdf5->h5_file_ptr, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_ID(hdf5->group_data);
    
    /* Print initial mesh */
    ps_print_geo_xdmf(grid);
    
    /* close everything (reopen when needed --- avoids having a truncated file) */
    HDF5_ERR(H5Gclose(hdf5->group_data));
    HDF5_ERR(H5Gclose(hdf5->group_elements));
    HDF5_ERR(H5Gclose(hdf5->group_nodes));
    HDF5_ERR(H5Gclose(hdf5->group_mesh));
    HDF5_ERR(H5Fclose(hdf5->h5_file_ptr));
    
    /* Create the XML tree */
    hdf5->doc = xmlNewDoc(BAD_CAST "1.0");
    hdf5->root_node = xmlNewNode(NULL, BAD_CAST "Xdmf");
    xmlNewProp(hdf5->root_node, BAD_CAST "Version", BAD_CAST "2.0");
    xmlNewProp(hdf5->root_node, BAD_CAST "xmlns:xi", BAD_CAST "http://www.w3.org/2001/XInclude");
    xmlDocSetRootElement(hdf5->doc, hdf5->root_node);
    
    /* create a DTD declaration */
    /*dtd =*/ xmlCreateIntSubset(hdf5->doc, BAD_CAST "Xdmf", NULL, BAD_CAST "Xdmf.dtd");
    
    /* Create DOMAIN */
    node1 = xmlNewChild(hdf5->root_node, NULL, BAD_CAST "Domain", NULL);
    
    /* GRID element for TEMPORAL collection */
    hdf5->collection_node = xmlNewChild(node1, NULL, BAD_CAST "Grid", NULL);
    xmlNewProp(hdf5->collection_node, BAD_CAST "CollectionType", BAD_CAST "Temporal");
    xmlNewProp(hdf5->collection_node, BAD_CAST "GridType", BAD_CAST "Collection");
    xmlNewProp(hdf5->collection_node, BAD_CAST "Name", BAD_CAST "Mesh");
    
    return;
}
#endif
