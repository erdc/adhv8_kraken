#include "global_header.h"

#ifdef _ADH_HDF5
#include "adh_hdf5.h"
#include <libxml/xpath.h>

// -------------------------------------------------
// -------------------------------------------------

void xdmf_init_particles(HDF5 *hdf5, int np, SPARTICLE *p, int npes, int myid, const char *proj_name) {
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
        snprintf(hdf5->xmf_filename, HDF5_NAME_LEN, "%s_ptout_p%d.xmf", proj_name, myid);
        snprintf(hdf5->h5_filename, HDF5_NAME_LEN, "%s_ptout_p%d.h5", proj_name, myid);
    } else {
        snprintf(hdf5->xmf_filename, HDF5_NAME_LEN, "%s_ptout_pall1.xmf", proj_name);
        snprintf(hdf5->h5_filename, HDF5_NAME_LEN, "%s_ptout_p0.h5", proj_name);
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
    
    /* Create the container group for the data sets */
    hdf5->group_data = H5Gcreate(hdf5->h5_file_ptr, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_ID(hdf5->group_data);
    
    /* Print initial mesh */
    ps_print_particle_t0_xdmf(hdf5,np,p);
    
    /* close everything (reopen when needed --- avoids having a truncated file) */
    HDF5_ERR(H5Gclose(hdf5->group_data));
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

// -------------------------------------------------
// -------------------------------------------------
// writes initial particle locations
void ps_print_particle_t0_xdmf(HDF5 *hdf5, int np, SPARTICLE *p) {
    
    int ii,jj,ielem, owner;
    
    hid_t dataspace, dataset, aid, attr;
    herr_t status;
    hsize_t dims[2];
    char string[MAXLINE];
    
    int npes = 1;
    int myid = 0;
    
    int compression = 0;
    hdf5->current_mesh++;
    
    // pack particle data
    double *node_coords = (double *) malloc(sizeof(double) * 3 * np);
    if (node_coords==NULL) {
        fprintf(stderr, "Could not malloc memory for nodes; %s:%d\n", __FILE__,__LINE__);
        exit(-1);
    }
    for (ii = 0; ii < np; ii++) {
        node_coords[3 * ii]     = p[ii].r.x;
        node_coords[3 * ii + 1] = p[ii].r.y;
        node_coords[3 * ii + 2] = p[ii].r.z;
    }
    
    /* Write out the nodal coordinates */
    dims[0] = np;
    dims[1] = 3;
    snprintf(string, MAXLINE, "%d", hdf5->current_mesh);
    xdmf_write_dataset(hdf5->group_nodes, string, H5T_NATIVE_DOUBLE, compression, 2, dims, node_coords);

    free(node_coords);
    return;
}

// -------------------------------------------------
// -------------------------------------------------

void ps_print_particle_xdmf(HDF5 *hdf5, SVECT *dpl, int np, double time) {
    
    double *index;                /* the index for scaling the error */
    int i, ie, ibl;               /* loop counter and end of loop */
    double conv1, convr, conv3, conv4;
    char string[MAXLINE];
    char buff1[MAXLINE] = "";
    char file_backup[MAXLINE];
    herr_t status;
    xmlNodePtr node1 = NULL, node2 = NULL;
    FILE *fp;
    
    /* Open the HDF file and groups */
    hdf5->h5_file_ptr = H5Fopen(hdf5->h5_filename, H5F_ACC_RDWR, H5P_DEFAULT);
    hdf5->group_mesh = H5Gopen(hdf5->h5_file_ptr, "Geometry", H5P_DEFAULT);
    hdf5->group_data = H5Gopen(hdf5->h5_file_ptr, "Data", H5P_DEFAULT);
    hdf5->group_nodes = H5Gopen(hdf5->group_mesh, "nodes", H5P_DEFAULT);
    
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
    sprintf(buff1, "%d", np);
    xmlNewProp(node1, BAD_CAST "NumberOfElements", BAD_CAST buff1);
    char elem_type[20];
    sprintf(elem_type, "Polyvertex");
    xmlNewProp(node1, BAD_CAST "Type", BAD_CAST elem_type);
    xmlAddChild(hdf5->current_grid_node, node1);
    
    
    /* GEOMETRY <Geomtery/> data */
    node1 = xmlNewNode(NULL, BAD_CAST "Geometry");
    xmlNewProp(node1, BAD_CAST "Type", BAD_CAST "XYZ");
    sprintf(buff1, "%s:/Geometry/nodes/%d", hdf5->h5_filename, hdf5->current_mesh);
    node2 = xmlNewChild(node1, NULL, BAD_CAST "DataItem", BAD_CAST buff1);
    xmlNewProp(node2, BAD_CAST "DataType", BAD_CAST "Float");
    sprintf(buff1, "%d 3", np);
    xmlNewProp(node2, BAD_CAST "Dimensions", BAD_CAST buff1);
    xmlNewProp(node2, BAD_CAST "Format", BAD_CAST "HDF");
    xmlNewProp(node2, BAD_CAST "Precision", BAD_CAST "8");
    xmlAddChild(hdf5->current_grid_node, node1);
    
    /* add Information element to supply np (number of nodes owned by this process) */
    node2 = xmlNewChild(node1, NULL, BAD_CAST "Information", NULL);
    xmlNewProp(node2, BAD_CAST "Name", BAD_CAST "np");
    sprintf(buff1, "%d", np);
    xmlNewProp(node2, BAD_CAST "Value", BAD_CAST buff1);
    
    xmlAddChild(hdf5->current_grid_node, node1);
    
    xdmf_particles_print_ts(hdf5,"Particle Displacement", dpl, np, time, hdf5->iprint_step, 1);
    
    /* close HDF groups and file */
    status = H5Gclose(hdf5->group_nodes);
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
    return;
}

// -------------------------------------------------
// -------------------------------------------------
void xdmf_particles_print_ts(
                             HDF5 *hdf5,
                             char *data_name,  /* String to be written as the data set identifier */
                             SVECT *dpl,       /* trajectory displacements */
                             int np,           /* gkc adding for elem centered data */
                             double time,      /* the time of the print */
                             int ntime,        /* index of time step */
                             int data_type     /* 0=scalar, 1=vector */ // should be 1 for now
) {
    
    int i;
    
    char type_string[MAXLINE] = "";
    sprintf(type_string, "Node");
    
    int *write_node_index = (int *)malloc(sizeof(int) * np);
    if (write_node_index==NULL){fprintf(stderr,"Memory Error: %s:%d",__FILE__,__LINE__);}
    for (i = 0; i < np; i++){
        write_node_index[i]=i;
    }
    
    // for now, no compression, but can be added as a feature later
    int compression = 0;
    
    // If this is the first time step, we need to create the HDF group
    hid_t group_id = UNSET_INT, dataspace, dataset, aid, attr, dataspace_times, filespace, cparms;
    hsize_t ndims, dims[2], maxdims[2], chunkdims[2], size[2], offset[2], slabdims[2];
    herr_t status;
    if (ntime == 0) {
        group_id = H5Gcreate(hdf5->group_data, data_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        // create "Times" data set which contains the time stamps
        dims[0] = 1;
        maxdims[0] = H5S_UNLIMITED;
        dataspace_times = H5Screate_simple(1, dims, maxdims);
        cparms = H5Pcreate(H5P_DATASET_CREATE);
        chunkdims[0] = 20;
        status = H5Pset_chunk(cparms, 1, chunkdims);
        dataset = H5Dcreate(group_id, "Times", H5T_NATIVE_DOUBLE, dataspace_times, H5P_DEFAULT, cparms, H5P_DEFAULT);
        status  = H5Pclose(cparms);
        status  = H5Sclose(dataspace_times);
        status  = H5Dclose(dataset);
    } else {
        group_id = H5Gopen(hdf5->group_data, data_name, 0);
    }
    
    // pack the initial positions
    double *node_coords = (double *) malloc(sizeof(double) * 3 * np);
    if (node_coords==NULL) {
        fprintf(stderr, "Could not malloc memory for nodes; %s:%d\n", __FILE__,__LINE__);
        exit(-1);
    }
    for (i = 0; i < np; i++) {
        node_coords[3 * i]     = dpl[write_node_index[i]].x;
        node_coords[3 * i + 1] = dpl[write_node_index[i]].y;
        node_coords[3 * i + 2] = dpl[write_node_index[i]].z;
    }
    
    // Use the time step index in the name of the HDF group
    char number[MAXLINE] = "";
    sprintf(number, "%d", ntime);
    
    // Write the HDF file
    ndims = 2;
    dims[0] = np;
    dims[1] = 3;
    xdmf_write_dataset(group_id, number, H5T_NATIVE_DOUBLE, compression, ndims, dims, node_coords);
    
    // Add attribute to data set which contains the corresponding mesh id
    xdmf_write_attribute(group_id, number, "Mesh_ID", H5T_NATIVE_INT, &(hdf5->current_mesh));
    
    // Add attribute to data set which contains the time stamp
    xdmf_write_attribute(group_id, number, "time", H5T_NATIVE_DOUBLE, &time);
    
    // write the time stamp to the Times dataset
    dataset = H5Dopen(group_id, "Times", H5P_DEFAULT);
    if (ntime != 0) {
        size[0] = ntime + 1;
        status = H5Dset_extent(dataset, size);
    }
    filespace = H5Dget_space(dataset);
    offset[0] = ntime;
    slabdims[0] = 1;
    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, slabdims, NULL);
    dims[0] = 1;
    dataspace_times = H5Screate_simple(1, dims, NULL);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace_times, filespace, H5P_DEFAULT, &time);
    status = H5Sclose(filespace);
    status = H5Dclose(dataset);
    status = H5Sclose(dataspace_times);
    status = H5Gclose(group_id);
    
    /* Add data entries to XML tree */
    char buff1[MAXLINE] = "";
    sprintf(buff1, "Vector");
    xmlNodePtr node1 = NULL, node2 = NULL;
    node1 = xmlNewNode(NULL, BAD_CAST "Attribute");
    xmlNewProp(node1, BAD_CAST "AttributeType", BAD_CAST buff1);
    xmlNewProp(node1, BAD_CAST "Center", BAD_CAST type_string);
    xmlNewProp(node1, BAD_CAST "Name", BAD_CAST data_name);
    sprintf(buff1, "%s:Data/%s/%d", hdf5->h5_filename, data_name, ntime);
    node2 = xmlNewChild(node1, NULL, BAD_CAST "DataItem", BAD_CAST buff1);
    xmlNewProp(node2, BAD_CAST "DataType", BAD_CAST "Float");
    sprintf(buff1, "%d 3", np);
    xmlNewProp(node2, BAD_CAST "Dimensions", BAD_CAST buff1);
    xmlNewProp(node2, BAD_CAST "Format", BAD_CAST "HDF");
    xmlNewProp(node2, BAD_CAST "Precision", BAD_CAST "8");
    xmlAddChild(hdf5->current_grid_node, node1);
    
    // free memory
    free(node_coords);
    free(write_node_index);
}

// -------------------------------------------------
// -------------------------------------------------

void xdmf_particle_finalize(HDF5 *hdf5, const char *proj_name, int npes, int myid, int flag) {
    
    herr_t status;
    //xmlDtdPtr dtd = NULL;
    xmlNodePtr node1 = NULL;
    xmlNodePtr node2 = NULL;
    xmlNodePtr temporal_node = NULL;    /* the temporal COLLECTION node */
    xmlNodePtr spatial_node = NULL; /* the spatial COLLECTION node */
    xmlNodePtr root_node = NULL, root_node_surf=NULL;
    xmlDocPtr doc_out = NULL;     /* the concatenated document */
    xmlDocPtr doc_in = NULL;      /* pointer to individual processor document */
    xmlNodePtr current_collection = NULL; /* pointer in current doc to collection grid node */
    xmlNodePtr current_temporal = NULL;   /* pointer into current doc at level to copy */
    xmlXPathContextPtr xpathCtx;
    xmlChar *xpathExpr;
    xmlXPathObjectPtr xpathObj;
    xmlAttr *props = NULL;
    int i,j;
    int igrid = 0;
    int icount = 0;
    int num_time_steps = 0;
    char *basename;
    char filename[128];
    FILE *fp;
    
    if (myid > 0) {
        return;
    }
    /* If we are serial, then we don't need to create the master XMF file */
    if (npes <= 1) {
        return;
    }
    
    xmlInitParser();
    //basename = malloc(128);
    //basename = proj_name;
    
    /*
     * Creates a new document, a node and set it as a root node
     */
    doc_out = xmlNewDoc(BAD_CAST "1.0");
    root_node = xmlNewNode(NULL, BAD_CAST "Xdmf");
    xmlNewProp(root_node, BAD_CAST "Version", BAD_CAST "2.0");
    xmlNewProp(root_node, BAD_CAST "xmlns:xi", BAD_CAST "http://www.w3.org/2001/XInclude");
    xmlDocSetRootElement(doc_out, root_node);
    
    /* loop over the xmf files */
    for (igrid = 0; igrid < npes; igrid++)  {
        /* Advance to next processor's file */
        sprintf(filename, "%s_ptout_p%d.xmf", proj_name, igrid);
        doc_in = xmlReadFile(filename, NULL, 0);
        
        /* Create xpath evaluation context */
        xpathCtx = xmlXPathNewContext(doc_in);
        if (xpathCtx == NULL)
        {
            //tl_error("XDMF Finalize Error: unable to create new XPath context\n");
            fprintf(stderr,"XDMF Particle Finalize Error: unable to create new XPath context; %s:%d\n", __FILE__, __LINE__);
            xmlFree(doc_out);
            exit(0);
        }
        xpathExpr = xmlCharStrdup("/Xdmf/Domain/Grid[@CollectionType=\"Temporal\"]");
        xpathObj = xmlXPathEvalExpression(xpathExpr, xpathCtx);
        if (xpathObj == NULL)
        {
            fprintf(stderr,"XDMF Finalize Error: unable to evaluate xpath expression; %s:%d\n", __FILE__, __LINE__);
            xmlXPathFreeContext(xpathCtx);
            xmlFreeDoc(doc_out);
            exit(0);
        }
        
        /* Gather all timestep nodes */
        current_collection = xpathObj->nodesetval->nodeTab[0];
        
        free(xpathExpr);
        free(xpathObj);
        
        /* the number of children at this point = no. of time steps */
        num_time_steps = xmlChildElementCount(current_collection);
        printf("num_time_Setps: %d\n",num_time_steps);
        
        /* set first time step collection in tree */
        current_temporal = current_collection->children;
        while (current_temporal->type != XML_ELEMENT_NODE)
        {
            current_temporal = current_temporal->next;
        }
        
        /* Create file template */
        if (igrid == 0)
        {
            /* Creates a DTD declaration. Isn't mandatory. */
            /*dtd =*/ xmlCreateIntSubset(hdf5->doc, BAD_CAST "Xdmf", NULL, BAD_CAST "Xdmf.dtd");
            
            /* xmlNewChild() creates a new node, which is "attached" as child node
             * of root_node node. */
            node1 = xmlNewChild(root_node, NULL, BAD_CAST "Domain", NULL);
            
            /* GRID element for TEMPORAL collection */
            temporal_node = xmlNewChild(node1, NULL, BAD_CAST "Grid", NULL);
            xmlNewProp(temporal_node, BAD_CAST "CollectionType", BAD_CAST "Temporal");
            xmlNewProp(temporal_node, BAD_CAST "GridType", BAD_CAST "Collection");
            xmlNewProp(temporal_node, BAD_CAST "Name", BAD_CAST "Mesh");
            
            /* START OF SPATIAL COLLECTION */
            for (i = 0; i < num_time_steps; i++)
            {
                spatial_node = xmlNewChild(temporal_node, NULL, BAD_CAST "Grid", NULL);
                xmlNewProp(spatial_node, BAD_CAST "CollectionType", BAD_CAST "Spatial");
                xmlNewProp(spatial_node, BAD_CAST "GridType", BAD_CAST "Collection");
            }
        }
        
        /* Start with first spatial node */
        spatial_node = temporal_node->children;
        while (spatial_node->type != XML_ELEMENT_NODE)
        {
            spatial_node = spatial_node->next;
        }
        
        for (i = 0; i < num_time_steps; i++)
        {
            /* Copy Time stamp above collection for XDMF animation in ParaView */
            if (igrid == 0)
            {
                node2 = xmlFirstElementChild(current_temporal);
                node1 = xmlCopyNode(node2, 1);
                xmlAddChild(spatial_node, node1);
            }
            
            node1 = xmlCopyNode(current_temporal, 1);
            xmlAddChild(spatial_node, node1);
            spatial_node = xmlNextElementSibling(spatial_node);
            current_temporal = xmlNextElementSibling(current_temporal);
        }
        
        xmlFreeDoc(doc_in);
    }
    
    /* Create Master file */
    sprintf(filename,"%s_ptout_pall_%d.xmf", proj_name, igrid);
    xmlSaveFormatFileEnc(filename, doc_out, "UTF-8", 1);
    
    /*free the document */
    xmlFreeDoc(doc_out);
    
    return;
}
#endif
