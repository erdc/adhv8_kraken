#include "global_header.h"
#ifdef _ADH_HDF5
#include "adh_hdf5.h"
#endif

void geo_xdmf_write(char *adh_root) {

#ifdef _ADH_HDF5

#ifdef _MESSG
    debug_initialize(MPI_COMM_WORLD);
#else
    debug_initialize();
#endif

    xmlNodePtr node1 = NULL, node2 = NULL;
    FILE *fp;

    char buff1[MAXLINE] = "";
    char elem_type[20];
    char string[MAXLINE];

    SMODEL *mod;
    HDF5 *hdf5 = &(mod->hdf5);
    mod = (SMODEL *) tl_alloc(sizeof(SMODEL), 1);
    /* standard model set up incase results need to get converted*/
    smodel_defaults(mod);
    mod->proc_flag = 1;
    mod->io = (SIO *) tl_alloc(sizeof(SIO), 1);
    sio_init(mod->io, adh_root);

    // check which geo file exists
    if (doesFileExist(mod->io->geo2d.filename)) {
        open_input_file(&(mod->io->geo2d), "2d geometry file", 0);
        mod->flag.SW2_FLOW = TRUE;
    } else if (doesFileExist(mod->io->geo3d.filename)) {
        open_input_file(&(mod->io->geo3d), "3d geometry file", 0);
        open_input_file( &(mod->io->face), "3d boundary face file", 0);
    } else {
        tl_error("2d or 3d geo file needed\n");
    }

    // open face file is 3D
    if (doesFileExist(mod->io->face.filename)) {
        ssw_3d_open_input(mod);
    }

    #ifdef _MESSG
        sgrid_alloc_init(&(mod->grid), mod->io, mod->proc_flag, mod->file_output, mod->model_comm);
    #else
        sgrid_alloc_init(&(mod->grid), mod->io, mod->proc_flag, mod->file_output);
    #endif
    
    
    if(mod->grid->ndim==3) mod->flag.SW3_FLOW=TRUE;
    SGRID *grid = mod->grid;
    SIO *io = mod->io;

    int nnode = grid->nnodes;
    int my_nnode = grid->my_nnodes;
    int nelems3d = grid->nelems3d;
    int nelems2d = grid->nelems2d;
    int nelem;
    xdmf_init(mod, 1, 0);
    printf("Initial grid converted\n");
    double time=0.;
    double c_time=0.;
    double outfact=1.;
    hdf5 = &(mod->hdf5); 
    if (doesFileExist(mod->io->fout_sw2_head.filename)||doesFileExist(mod->io->fout_sw3_depth.filename)) {
        outfact=mod->series_out->outfact;
    }else{
        /* update XML tree <Grid/> */
        mod->hdf5.current_grid_node = xmlNewChild(mod->hdf5.collection_node, NULL, BAD_CAST "Grid", NULL);
        xmlNewProp(hdf5->current_grid_node, BAD_CAST "GridType", BAD_CAST "Uniform");
        
        if (grid->haveTets == TRUE) {
            nelem = grid->nelems3d;
            sprintf(elem_type, "Tetrahedron");
        } else if (grid->havePrisms == TRUE) {
            nelem = grid->nelems3d;
            sprintf(elem_type, "Wedge");
        } else if (grid->haveTris == TRUE) {
            nelem = grid->nelems2d;
            sprintf(elem_type, "Triangle");
        } else if (grid->haveQuads == TRUE) {
            nelem = grid->nelems2d;
            sprintf(elem_type, "Quadrilaterals");
        }
        
        /* add time stamp <Time/> */
        node1 = xmlNewNode(NULL, BAD_CAST "Time");
        sprintf(buff1, "%f", c_time);
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
        if (grid->haveTets == TRUE) {
            sprintf(buff1, "%d 4", grid->nelems3d);
        } else if (grid->havePrisms == TRUE) {
            sprintf(buff1, "%d 6", grid->nelems3d);
        } else if (grid->haveTris == TRUE) {
            sprintf(buff1, "%d 3", grid->nelems2d);
        } else if (grid->haveQuads == TRUE) {
            sprintf(buff1, "%d 4", grid->nelems2d);
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
        sprintf(buff1, "%d 3", grid->nnodes);
        xmlNewProp(node2, BAD_CAST "Dimensions", BAD_CAST buff1);
        xmlNewProp(node2, BAD_CAST "Format", BAD_CAST "HDF");
        xmlNewProp(node2, BAD_CAST "Precision", BAD_CAST "8");
        xmlAddChild(hdf5->current_grid_node, node1);
        
        /* add Information element to supply grid->my_nnodes (number of nodes owned by this process) */
        node2 = xmlNewChild(node1, NULL, BAD_CAST "Information", NULL);
        xmlNewProp(node2, BAD_CAST "Name", BAD_CAST "my_nnode");
        sprintf(buff1, "%d", grid->my_nnodes);
        xmlNewProp(node2, BAD_CAST "Value", BAD_CAST buff1);
        xmlAddChild(hdf5->current_grid_node, node1);
        fp = io_fopen(hdf5->xmf_filename, "r", TRUE);
        
        if (fp != NULL) {
            fclose(fp);
            snprintf(string, MAXLINE - 1, "%s_back", hdf5->xmf_filename);
            rename(hdf5->xmf_filename, string);
        }
        /* Write the xmf file */
        
        xmlSaveFormatFileEnc(hdf5->xmf_filename, hdf5->doc, "UTF-8", 1);
        if (mod->flag.SW3_FLOW) {
            /* update XML tree <Grid/> */
            hdf5->current_grid_node_surf = xmlNewChild(hdf5->collection_node_surf, NULL, BAD_CAST "Grid", NULL);
            xmlNewProp(hdf5->current_grid_node_surf, BAD_CAST "GridType", BAD_CAST "Uniform");
            
            /* add time stamp <Time/> */
            node1 = xmlNewNode(NULL, BAD_CAST "Time");
            sprintf(buff1, "%f", c_time);
            xmlNewProp(node1, BAD_CAST "Value", BAD_CAST buff1);
            xmlAddChild(hdf5->current_grid_node_surf, node1);
            
            /* TOPOLOGY <Topology/> data */
            node1 = xmlNewNode(NULL, BAD_CAST "Topology");
            sprintf(buff1, "%d", nelem);
            xmlNewProp(node1, BAD_CAST "NumberOfElements", BAD_CAST buff1);
            xmlNewProp(node1, BAD_CAST "Type", BAD_CAST elem_type);
            sprintf(buff1, "%s:/Geometry/elements/%d", hdf5->h5_surf_filename, hdf5->current_mesh);
            node2 = xmlNewChild(node1, NULL, BAD_CAST "DataItem", BAD_CAST buff1);
            xmlNewProp(node2, BAD_CAST "DataType", BAD_CAST "Int");
            sprintf(buff1, "%d 3", nelem); /* Compulsorily for Triangles */
            xmlNewProp(node2, BAD_CAST "Dimensions", BAD_CAST buff1);
            xmlNewProp(node2, BAD_CAST "Format", BAD_CAST "HDF");
            xmlAddChild(hdf5->current_grid_node_surf, node1);
            
            /* GEOMETRY <Geomtery/> data */
            node1 = xmlNewNode(NULL, BAD_CAST "Geometry");
            xmlNewProp(node1, BAD_CAST "Type", BAD_CAST "XYZ");
            sprintf(buff1, "%s:/Geometry/nodes/%d", hdf5->h5_surf_filename, hdf5->current_mesh);
            node2 = xmlNewChild(node1, NULL, BAD_CAST "DataItem", BAD_CAST buff1);
            xmlNewProp(node2, BAD_CAST "DataType", BAD_CAST "Float");
            sprintf(buff1, "%d 3", grid->nnodes_sur);
            xmlNewProp(node2, BAD_CAST "Dimensions", BAD_CAST buff1);
            xmlNewProp(node2, BAD_CAST "Format", BAD_CAST "HDF");
            xmlNewProp(node2, BAD_CAST "Precision", BAD_CAST "8");
            xmlAddChild(hdf5->current_grid_node_surf, node1);
            
            /* add Information element to supply grid->my_nnodes (number of nodes owned by this process) */
            node2 = xmlNewChild(node1, NULL, BAD_CAST "Information", NULL);
            xmlNewProp(node2, BAD_CAST "Name", BAD_CAST "my_nnode");
            sprintf(buff1, "%d", grid->my_nnodes_sur);
            xmlNewProp(node2, BAD_CAST "Value", BAD_CAST buff1);
            xmlAddChild(hdf5->current_grid_node_surf, node1);
            fp = io_fopen(hdf5->xmf_surf_filename, "r", TRUE);
            if (fp != NULL) {
                fclose(fp);
                snprintf(string, MAXLINE - 1, "%s_back", hdf5->xmf_surf_filename);
                rename(hdf5->xmf_surf_filename, string);
            }
            xmlSaveFormatFileEnc(hdf5->xmf_surf_filename, hdf5->doc_surf, "UTF-8", 1);
        }
        xdmf_finalize(&(mod->hdf5), mod->io, 1, 0, mod->flag.SW3_FLOW);
        printf("XDMF files finalized\n");
        //   smodel_free(mod);printf("Here %d\n",,__LINE__);
        sgrid_free(mod->grid);
        
    }
#endif
}



