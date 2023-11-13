#define HDF5_INIT
#include "global_header.h"
#ifdef _ADH_HDF5
#include "adh_hdf5.h"
#endif
/* ADH Version 2.0.0 6-04 */
/* Initialize use of XMDF/HDF5 */

#ifdef _ADH_HDF5
void xdmf_init(SMODEL *mod, int npes, int myid)
{
    HDF5 *hdf5= &(mod->hdf5);
    hid_t aid, attr;
    hsize_t dims[1];
    FILE *fp;
    xmlDtdPtr dtd = NULL;
    xmlNodePtr node1 = NULL;
    char file_backup[HDF5_NAME_LEN];
    SIO *io = mod->io;
    
    /* If adaption is on, then XDMF must print out new grids */
    
    hdf5->iprint_step = 0;
    
    /* Initialize time step counter and mesh counter */
    hdf5->current_mesh = -1;
    hdf5->current_ts = 0;
    
    /* Form filenames for the HDF5 file and the XMF file
     If we are running in parallel, add "_myid" to filenames */
    if (npes > 1) {
        snprintf(hdf5->xmf_filename, HDF5_NAME_LEN, "%s_p%d.xmf", io->proj_name, myid);
        snprintf(hdf5->h5_filename, HDF5_NAME_LEN, "%s_p%d.h5", io->proj_name, myid);
    } else {
        snprintf(hdf5->xmf_filename, HDF5_NAME_LEN, "%s_pall1.xmf", io->proj_name);
        snprintf(hdf5->h5_filename, HDF5_NAME_LEN, "%s_p0.h5", io->proj_name);
    }
    
    /* test that we can create the xmf file */
    fp = io_fopen(hdf5->xmf_filename, "w", TRUE);
    if (fp == NULL){tl_error("Unable to create xmf file");}
    fclose(fp);
    
    /* If the *h5 file exists, then move it to backup copy */
    fp = fopen(hdf5->h5_filename, "r");
    if (fp != NULL)
    {
        fclose(fp);
        snprintf(file_backup, HDF5_NAME_LEN - 1, "%s_back", hdf5->h5_filename);
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
    hdf5->group_ghost_nodes = H5Gcreate(hdf5->group_mesh, "ghost_nodes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_ID(hdf5->group_ghost_nodes);
    hdf5->group_elements = H5Gcreate(hdf5->group_mesh, "elements", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_ID(hdf5->group_elements);
    hdf5->group_orig_nd = H5Gcreate(hdf5->group_mesh, "orig_node_number", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_ID(hdf5->group_orig_nd);
    hdf5->group_mat_id = H5Gcreate(hdf5->group_mesh, "mat_id", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_ID(hdf5->group_mat_id);
    hdf5->group_node_map = H5Gcreate(hdf5->group_mesh, "node_map", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_ID(hdf5->group_node_map);
    
    /* Create the container group for the data sets */
    hdf5->group_data = H5Gcreate(hdf5->h5_file_ptr, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_ID(hdf5->group_data);
    
    /* Print initial mesh */
    ps_print_geo_xdmf(mod);
    
    /* close everything (reopen when needed --- avoids having a truncated file) */
    HDF5_ERR(H5Gclose(hdf5->group_data));
    HDF5_ERR(H5Gclose(hdf5->group_node_map));
    HDF5_ERR(H5Gclose(hdf5->group_mat_id));
    HDF5_ERR(H5Gclose(hdf5->group_orig_nd));
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
    dtd = xmlCreateIntSubset(hdf5->doc, BAD_CAST "Xdmf", NULL, BAD_CAST "Xdmf.dtd");
    
    /* Create DOMAIN */
    node1 = xmlNewChild(hdf5->root_node, NULL, BAD_CAST "Domain", NULL);
    
    /* GRID element for TEMPORAL collection */
    hdf5->collection_node = xmlNewChild(node1, NULL, BAD_CAST "Grid", NULL);
    xmlNewProp(hdf5->collection_node, BAD_CAST "CollectionType", BAD_CAST "Temporal");
    xmlNewProp(hdf5->collection_node, BAD_CAST "GridType", BAD_CAST "Collection");
    xmlNewProp(hdf5->collection_node, BAD_CAST "Name", BAD_CAST "Mesh");
    
    
    /***************************************************************************************************************/
    /* Writing 2D mesh data */
    /* Handle the xdmf file for the surface mesh independently */
    if (mod->flag.SW3_FLOW
#ifdef _ADH_GROUNDWATER
        || mod->flag.GW_FLOW
#endif
        )
    {
        if (npes > 1)
        {
            snprintf(hdf5->xmf_surf_filename, HDF5_NAME_LEN, "%s_surf_p%d.xmf", io->proj_name, myid);
            snprintf(hdf5->h5_surf_filename, HDF5_NAME_LEN, "%s_surf_p%d.h5", io->proj_name, myid);
        }
        else
        {
            snprintf(hdf5->xmf_surf_filename, HDF5_NAME_LEN, "%s_surf_pall1.xmf", io->proj_name);
            snprintf(hdf5->h5_surf_filename, HDF5_NAME_LEN, "%s_surf_p0.h5", io->proj_name);
        }
        
        /* test that we can create the xmf file */
        fp = io_fopen(hdf5->xmf_surf_filename, "w", TRUE);
        if (fp == NULL)
        {
            tl_error("Unable to create xmf file");
        }
        fclose(fp);
        
        /* If the *h5 file exists, then move it to backup copy */
        fp = fopen(hdf5->h5_surf_filename, "r");
        if (fp != NULL)
        {
            fclose(fp);
            snprintf(file_backup, HDF5_NAME_LEN - 1, "%s_back", hdf5->h5_surf_filename);
            rename(hdf5->h5_surf_filename, file_backup);
        }
        hdf5->h5_surf_file_ptr = H5Fcreate(hdf5->h5_surf_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        HDF5_ID(hdf5->h5_surf_file_ptr);
        
        /* Create the container group for the meshes */
        hdf5->group_surf_mesh = H5Gcreate(hdf5->h5_surf_file_ptr, "/Geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        HDF5_ID(hdf5->group_surf_mesh);
        hdf5->group_surf_nodes = H5Gcreate(hdf5->group_surf_mesh, "nodes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        HDF5_ID(hdf5->group_surf_nodes);
        hdf5->group_surf_ghost_nodes = H5Gcreate(hdf5->group_surf_mesh, "ghost_nodes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        HDF5_ID(hdf5->group_surf_ghost_nodes);
        hdf5->group_surf_elements = H5Gcreate(hdf5->group_surf_mesh, "elements", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        HDF5_ID(hdf5->group_surf_elements);
        //hdf5->group_surf_orig_nd = H5Gcreate(hdf5->group_surf_mesh, "orig_node_number", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //HDF5_ID(hdf5->group_surf_orig_nd);
        //hdf5->group_surf_mat_id = H5Gcreate(hdf5->group_surf_mesh, "mat_id", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //HDF5_ID(hdf5->group_surf_mat_id);
        hdf5->group_surf_node_map = H5Gcreate(hdf5->group_surf_mesh, "node_map", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        HDF5_ID(hdf5->group_surf_node_map);
        
        //hdf5->group_surf_nodes = H5Gcreate(hdf5->group_surf_mesh, "nodes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //HDF5_ID(hdf5->group_surf_nodes);
        //hdf5->group_surf_elements = H5Gcreate(hdf5->group_surf_mesh, "elements", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        //HDF5_ID(hdf5->group_surf_elements);
        
        /* Create the container group for the data sets */
        hdf5->group_surf_data = H5Gcreate(hdf5->h5_surf_file_ptr, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        HDF5_ID(hdf5->group_surf_data);
        
        /* Print initial surface mesh */
        ps_print_surf_geo_xdmf(mod);
        
        /* close everything (reopen when needed --- avoids having a truncated file) */
        HDF5_ERR(H5Gclose(hdf5->group_surf_data));
        HDF5_ERR(H5Gclose(hdf5->group_surf_node_map));
        //HDF5_ERR(H5Gclose(hdf5->group_surf_mat_id));
        //HDF5_ERR(H5Gclose(hdf5->group_surf_orig_nd));
        HDF5_ERR(H5Gclose(hdf5->group_surf_elements));
        HDF5_ERR(H5Gclose(hdf5->group_surf_ghost_nodes));
        HDF5_ERR(H5Gclose(hdf5->group_surf_nodes));
        HDF5_ERR(H5Gclose(hdf5->group_surf_mesh));
        HDF5_ERR(H5Fclose(hdf5->h5_surf_file_ptr));
        
        /* Create the XML tree */
        hdf5->doc_surf = xmlNewDoc(BAD_CAST "1.0");
        hdf5->root_node_surf = xmlNewNode(NULL, BAD_CAST "Xdmf");
        xmlNewProp(hdf5->root_node_surf, BAD_CAST "Version", BAD_CAST "2.0");
        xmlNewProp(hdf5->root_node_surf, BAD_CAST "xmlns:xi", BAD_CAST "http://www.w3.org/2001/XInclude");
        xmlDocSetRootElement(hdf5->doc_surf, hdf5->root_node_surf);
        
        /* create a DTD declaration */
        dtd = xmlCreateIntSubset(hdf5->doc_surf, BAD_CAST "Xdmf", NULL, BAD_CAST "Xdmf.dtd");
        
        /* Create DOMAIN */
        node1 = xmlNewChild(hdf5->root_node_surf, NULL, BAD_CAST "Domain", NULL);
        
        /* GRID element for TEMPORAL collection */
        hdf5->collection_node_surf = xmlNewChild(node1, NULL, BAD_CAST "Grid", NULL);
        xmlNewProp(hdf5->collection_node_surf, BAD_CAST "CollectionType", BAD_CAST "Temporal");
        xmlNewProp(hdf5->collection_node_surf, BAD_CAST "GridType", BAD_CAST "Collection");
        xmlNewProp(hdf5->collection_node_surf, BAD_CAST "Name", BAD_CAST "Mesh");
    }
    
    
    /* Gajanan gkc adding this in order to ensure the first density output (t=0) is correct. */
    /* This, however, is probably not the appropriate location for adding all these lines. */
    if (mod->flag.SW2_FLOW){
        switch(mod->flag.BAROCLINIC) {
            case 1:
                tl_density_calculator_metric(mod->density, NULL, 1., mod->con[mod->salinity_id].concentration, mod->con[mod->salinity_id].property[0],  mod->grid->nnodes, mod->sw->d2->density, mod->flag.EOS);
                break;
            case 10:
                tl_density_calculator_metric(mod->density, mod->con[mod->temperature_id].concentration, 1., NULL, mod->con[mod->temperature_id].property[0], mod->grid->nnodes, mod->sw->d2->density, mod->flag.EOS);
                break;
            case 11:
                tl_density_calculator_metric(mod->density, mod->con[mod->temperature_id].concentration, 1., mod->con[mod->salinity_id].concentration, mod->con[mod->salinity_id].property[0], mod->grid->nnodes, mod->sw->d2->density, mod->flag.EOS);
                break;
        }
        
        int inode=0;
        for (inode = 0; inode < mod->grid->nnodes; inode++) {
            mod->sw->d2->density[inode] /= mod->density;
            mod->sw->d2->density[inode] -= 1.;
        }
    }
    else if (mod->flag.SW3_FLOW
#ifdef _ADH_GROUNDWATER
             || mod->flag.GW_FLOW
#endif
             ){
        switch(mod->flag.BAROCLINIC) {
            case 1:
                tl_density_calculator_metric(mod->density, NULL, 1., mod->con[mod->salinity_id].concentration, mod->con[mod->salinity_id].property[0], mod->grid->nnodes, mod->sw->d3->density, mod->flag.EOS);
                break;
            case 10:
                tl_density_calculator_metric(mod->density, mod->con[mod->temperature_id].concentration, 1., NULL, mod->con[mod->temperature_id].property[0], mod->grid->nnodes, mod->sw->d3->density, mod->flag.EOS);
                break;
            case 11:
                tl_density_calculator_metric(mod->density, mod->con[mod->temperature_id].concentration, 1., mod->con[mod->salinity_id].concentration, mod->con[mod->salinity_id].property[0], mod->grid->nnodes, mod->sw->d3->density, mod->flag.EOS);
                break;
        }
    }
    
    return;
}
#endif
