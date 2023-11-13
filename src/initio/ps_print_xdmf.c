/* write the solution to a XDMF file */

#include "global_header.h"

void ps_print_xdmf(
                   SMODEL *mod,
                   double time,       /* the time of the print */
                   double outfact
                   )
{
#ifdef _ADH_HDF5
    SVECT *vel_dim;                /* the dimensional velocity (NS) */
    SVECT *vel_grad_dim;           /* the dimensional gradients of each vel. comp. */
    double *prs_dim;              /* the dimensional pressure (NS) */
    double *dpl_dim;              /* the dimensional displacement (NS), only used if moving grid */
    double *spd_dim;              /* the dimensional grid speed (NS), only used if moving grid */
    double *conc_dim;             /* the dimensional concentration */
    double *bl_conc_dim;          /* the dimensional concentration */
    double *pconc_dim;            /* the dimensional concentration */
    double *pbl_conc_dim;         /* the dimensional concentration */
    double *index;                /* the index for scaling the error */
    double *node_error;
    double *nodal_sat;            /* nodal saturation values */
    SVECT *nodal_gw_vel;           /* nodal velocities for printing */
    double c_time;                /* converted time */
    int itrns;                    /* loop counter over the number of transported items */
    int sed_number;               /* the clay transport index */
    int i, ie, ibl;               /* loop counter and end of loop */
    double conv1, convr, conv3, conv4;
    char string[MAXLINE];
    char buff1[MAXLINE] = "";
    char file_backup[MAXLINE];
    char elem_type[20];
    herr_t status;
    xmlNodePtr node1 = NULL, node2 = NULL;
    FILE *fp;
    int nelem;
    
    /* Gajanan gkc adding missing declarations */
    SGRID *grid = mod->grid;
    SIO *io = mod->io;
    HDF5 *hdf5 = &(mod->hdf5);
    int nnode = grid->nnodes;
    int my_nnode = grid->my_nnodes;
    int nelems3d = grid->nelems3d;
    int nelems2d = grid->nelems2d;
    
    int REF_FLAG_RECV = NO;
    int UNREF_FLAG_RECV = NO;
#ifdef _MESSG
    int REF_FLAG = mod->flag.GRID_REFINED;
    int UNREF_FLAG = mod->flag.GRID_UNREFINED;
    MPI_Allreduce(&REF_FLAG, &REF_FLAG_RECV, 1, MPI_INT, MPI_MAX, mod->grid->smpi->ADH_COMM);
    MPI_Allreduce(&UNREF_FLAG, &UNREF_FLAG_RECV, 1, MPI_INT, MPI_MAX, mod->grid->smpi->ADH_COMM);
#else
    REF_FLAG_RECV = mod->flag.GRID_REFINED;
    UNREF_FLAG_RECV = mod->flag.GRID_UNREFINED;
#endif
    
    /* Open the HDF file and groups */
    hdf5->h5_file_ptr = H5Fopen(hdf5->h5_filename, H5F_ACC_RDWR, H5P_DEFAULT);
    hdf5->group_mesh = H5Gopen(hdf5->h5_file_ptr, "Geometry", H5P_DEFAULT);
    hdf5->group_data = H5Gopen(hdf5->h5_file_ptr, "Data", H5P_DEFAULT);
    hdf5->group_elements = H5Gopen(hdf5->group_mesh, "elements", H5P_DEFAULT);
    hdf5->group_nodes = H5Gopen(hdf5->group_mesh, "nodes", H5P_DEFAULT);
    hdf5->group_ghost_nodes = H5Gopen(hdf5->group_mesh, "ghost_nodes", H5P_DEFAULT);
    hdf5->group_orig_nd = H5Gopen(hdf5->group_mesh, "orig_node_number", H5P_DEFAULT);
    hdf5->group_mat_id = H5Gopen(hdf5->group_mesh, "mat_id", H5P_DEFAULT);
    hdf5->group_node_map = H5Gopen(hdf5->group_mesh, "node_map", H5P_DEFAULT);
    
    /* convert time to specified output units - dss */
    c_time = time * outfact;
    //printf("Current time = %f",c_time);
    
    /* Set XDMF Element Type */
    if (grid->haveTets == TRUE)
    {
        nelem = grid->nelems3d;
        sprintf(elem_type, "Tetrahedron");
    }
    else if (grid->havePrisms == TRUE)
    {
        nelem = grid->nelems3d;
        sprintf(elem_type, "Prism");
    }
    else if (grid->haveTris == TRUE)
    {
        nelem = grid->nelems2d;
        sprintf(elem_type, "Triangle");
    }
    
    /* if grid was adapted, things are more complicated */
    if (mod->flag.ADAPTED_THE_GRID == YES || REF_FLAG_RECV == YES || UNREF_FLAG_RECV == YES) {
        //tag();
        ps_print_geo_xdmf(mod);
    }
    
    
    /* update XML tree <Grid/> */
    hdf5->current_grid_node = xmlNewChild(hdf5->collection_node, NULL, BAD_CAST "Grid", NULL);
    xmlNewProp(hdf5->current_grid_node, BAD_CAST "GridType", BAD_CAST "Uniform");
    
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
    if (grid->haveTets == TRUE)
    {
        sprintf(buff1, "%d 4", grid->nelems3d);
    }
    else if (grid->havePrisms == TRUE)
    {
        sprintf(buff1, "%d 6", grid->nelems3d);
    }
    else if (grid->haveTris == TRUE)
    {
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
    
    if (grid->nnodes > grid->my_nnodes)
    {
        /* GHOST NODE <Set/> information */
        node1 = xmlNewNode(NULL, BAD_CAST "Set");
        xmlNewProp(node1, BAD_CAST "Ghost", BAD_CAST "1");
        xmlNewProp(node1, BAD_CAST "SetType", BAD_CAST "Node");
        sprintf(buff1, "%s:/Geometry/ghost_nodes/%d", hdf5->h5_filename, hdf5->current_mesh);
        node2 = xmlNewChild(node1, NULL, BAD_CAST "DataItem", BAD_CAST buff1);
        xmlNewProp(node2, BAD_CAST "DataType", BAD_CAST "Int");
        sprintf(buff1, "%d", grid->nnodes - grid->my_nnodes);
        xmlNewProp(node2, BAD_CAST "Dimensions", BAD_CAST buff1);
        xmlNewProp(node2, BAD_CAST "Format", BAD_CAST "HDF");
    }
    xmlAddChild(hdf5->current_grid_node, node1);
    
    
    /*************************************************************************************************************************************************/
    if (mod->flag.SW2_FLOW && (mod->flag.GW_FLOW==OFF)) {
        SSW_2D *sw = mod->sw->d2;
        double *wse = (double *) tl_alloc(sizeof(double), grid->nnodes);
        for (i=0; i<grid->nnodes; i++){
            wse[i] = sw->head[i] + grid->node[i].z;
        }
        ps_print_ts_xdmf(hdf5,"Surface_Elevation", wse, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        wse = (double *) tl_free(sizeof(double), grid->nnodes, wse);
        /* prints the time step */
        ps_print_ts_xdmf(hdf5,"Depth", sw->head, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        ps_print_ts_xdmf(hdf5,"Velocity", sw->vel, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, VECTOR2D_DATA);
        if (mod->flag.BAROCLINIC){
            double *density = (double *) tl_alloc(sizeof(double), grid->nnodes);
            for (i=0; i<grid->nnodes; i++){
                density[i] = (sw->density[i]+1.0)*mod->density;
            }
            ps_print_ts_xdmf(hdf5,"Density", density, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            density = (double *) tl_free(sizeof(double), grid->nnodes, density);
        }
        if (mod->tau_temporal != 0.0){
            ps_print_ts_xdmf(hdf5,"Old_Depth", sw->old_head, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            ps_print_ts_xdmf(hdf5,"Old_Velocity", sw->old_vel, grid->nnodes, c_time,hdf5->iprint_step, NODE_CENTERED, VECTOR2D_DATA);
        }
        if ((mod->flag.WIND == ON || mod->flag.WIND_LIBRARY == ON) && mod->file_output.wind == ON) {
            ps_print_ts_xdmf(hdf5,"Wind_Stress", sw->winds, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, VECTOR2D_DATA);
        }
    }
    else if (mod->flag.SW3_FLOW) {
        /* prints the time step */
        SSW_3D *sw = mod->sw->d3;
        //ID_LIST_ITEM *ptr;
        ////double *node_dep3d = (double *) tl_alloc(sizeof(double), grid->nnodes);
        //double *node_dep2d = (double *) tl_alloc(sizeof(double), grid->nnodes);
        //double *wse2d = (double *) tl_alloc(sizeof(double), grid->nnodes);
        //for (i=0; i<grid->nnodes_sur; i++){
        //    int inode=-1, botnode=-1;
        //    ptr = grid->vertical_list[i];
        //    int surfnode = ptr->id;
        //    while (ptr->next!=NULL){
        //        botnode = ptr->id;
        //        ptr= ptr->next;
        //    }
        //    ptr = grid->vertical_list[i];
        //    while (ptr->next!=NULL){
        //        inode = ptr->id;
        //        //node_dep3d[inode] = sw->displacement[inode] + grid->node[inode].z -(sw->displacement[botnode] + grid->node[botnode].z);
        //        node_dep2d[inode] = sw->displacement[surfnode] + grid->node[surfnode].z -(sw->displacement[botnode] + grid->node[botnode].z);
        //        wse2d[inode] = sw->displacement[surfnode] + grid->node[surfnode].z;
        //        ptr= ptr->next;
        //    }
        //}
        // gkc - I'm wondering how to choose names here. What I'm doing with node_dep2d and wse2d is just assigning the same depth and wse along the entire column of nodes to print the solution. We may want to just delete it entirely
        //ps_print_ts_xdmf(hdf5,"3D_Depth", node_dep3d, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        //ps_print_ts_xdmf(hdf5,"2D_Surface_Elevation", wse2d, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        //ps_print_ts_xdmf(hdf5,"2D_Depth", node_dep2d, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        ////node_dep3d = (double *) tl_free(sizeof(double), grid->nnodes, node_dep3d);
        //node_dep2d = (double *) tl_free(sizeof(double), grid->nnodes, node_dep2d);
        //wse2d = (double *) tl_free(sizeof(double), grid->nnodes, wse2d);
        ps_print_ts_xdmf(hdf5,"Velocity", sw->vel, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, VECTOR3D_DATA);
        ps_print_ts_xdmf(hdf5,"Displacement", sw->displacement, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        //ps_print_ts_xdmf(hdf5,"Grid_Speed", sw->grid_speed, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        //ps_print_ts_xdmf(hdf5,"Pressure", sw->prs, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        //ps_print_ts_xdmf(hdf5,"Nodal_Flux", sw->vertical_node_flux, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        if (mod->flag.BAROCLINIC){
            ps_print_ts_xdmf(hdf5,"Density", sw->density, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        }
        if (mod->tau_temporal!=0.0){
            ps_print_ts_xdmf(hdf5,"Old_Velocity", sw->old_vel, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, VECTOR3D_DATA);
            ps_print_ts_xdmf(hdf5,"Old_Displacement", sw->old_displacement, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        }
    }
    else if (mod->flag.NS3_FLOW) {
        SNS_3D *ns = mod->ns->d3;
        ps_print_ts_xdmf(hdf5,"Velocity", ns->vel, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, VECTOR3D_DATA);
        ps_print_ts_xdmf(hdf5,"Displacement", ns->displacement, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        //ps_print_ts_xdmf(hdf5,"Grid_Speed", sw->grid_speed, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        //ps_print_ts_xdmf(hdf5,"Pressure", sw->prs, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
    }
#ifdef _ADH_GROUNDWATER
    else if (mod->flag.GW_FLOW) {
        SGW * gw = mod->sgw;
        assert(gw);
        /** pressure head -- assuming already calculated **/
        
        ps_print_ts_xdmf(hdf5,"Pressure Head", gw->gw_phead, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        
        
        /** total head -- assuming already calculated **/
        double * nodal_tmp = (double *) tl_alloc(sizeof(double), grid->nnodes);
        int * nodal_cnt    = (int *) tl_alloc(sizeof(int), grid->nnodes);
        for (i=0; i<grid->nnodes; i++){
            nodal_tmp[i] = gw->gw_phead[i]+grid->node[i].z;
        }
        ps_print_ts_xdmf(hdf5,"Total Head", nodal_tmp, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        
        /** density -- assuming already calculated **/
        ps_print_ts_xdmf(hdf5,"Density", gw->gw_density, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        /** average saturations at nodes **/
        sgw_nodal_sat_avg(gw,grid,grid->nnodes, nodal_tmp, nodal_cnt);
        ps_print_ts_xdmf(hdf5,"Saturation", nodal_tmp, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        /** error -- assuming already calculated **/
        ps_print_ts_xdmf(hdf5,"Estimated GW Error", gw->error, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        ps_print_ts_xdmf(hdf5,"GW Error", grid->elem_error, grid->nelems3d, c_time, hdf5->iprint_step, ELEM_3D_CENTERED, SCALAR_DATA);
        /** Darcy flux on elements */
        ps_print_ts_xdmf(hdf5,"Darcy Flux", gw->elem_gw_flux, grid->nelems3d, c_time, hdf5->iprint_step, ELEM_3D_CENTERED, VECTOR3D_DATA);
        /** Nodally averaged Velocity at nodes calculated from Darcy flux */
        SVECT * nodal_vel = (SVECT *) tl_alloc(sizeof(SVECT), grid->nnodes);
        //tl_elem3d_SVECT_var_to_nodal_avg(grid->nnodes, grid->nelems3d, grid->elem3d, gw->elem_gw_flux, nodal_vel, grid->smpi);
        sgw_elem3d_flux_to_nodal_vel(grid, mod->mat, gw, nodal_vel);
        ps_print_ts_xdmf(hdf5,"Velocity", nodal_vel, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, VECTOR3D_DATA);
        nodal_vel = (SVECT *) tl_free(sizeof(SVECT), grid->nnodes, nodal_vel);
        
        nodal_tmp = (double *) tl_free(sizeof(double), grid->nnodes, nodal_tmp);
        nodal_cnt = (int *) tl_free(sizeof(int), grid->nnodes, nodal_cnt);
    }
#endif
    
    if (mod->ntransport > 0) { /* Mostly same for both 2D and 3D SW models */
        SCON *con=NULL;
        con = mod->con;
        conc_dim = (double *) tl_alloc(sizeof(double), grid->nnodes);
        pconc_dim = (double *) tl_alloc(sizeof(double), grid->nnodes);
        for (itrns = 0; itrns < mod->ntransport; itrns++) {
            sarray_scale_dbl(conc_dim, con[itrns].concentration, con[itrns].property[0], grid->nnodes);
            sarray_scale_dbl(pconc_dim, con[itrns].concentration, con[itrns].property[0], grid->nnodes);
            if (con[itrns].type == SAL) {
                snprintf(string, MAXLINE, "Salt_Concentration");
            }
            else if (con[itrns].type == TMP) {
                snprintf(string, MAXLINE, "Temperature_Distribution");
            }
            else {
                snprintf(string, MAXLINE, "Concentration %d", itrns+1);
            }
            ps_print_ts_xdmf(hdf5,string, conc_dim, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            if(mod->tau_temporal != 0.0){
                if (con[itrns].type == SAL) {
                    snprintf(string, MAXLINE, "Old_Salt_Concentration");
                }
                else if (con[itrns].type == TMP) {
                    snprintf(string, MAXLINE, "Old_Temperature_Distribution");
                }
                else {
                    snprintf(string, MAXLINE, "Old_Concentration %d", itrns+1);
                }
                ps_print_ts_xdmf(hdf5, string, pconc_dim, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            }
        }
        
#ifdef _SEDIMENT
        if (mod->flag.SEDIMENT && mod->flag.SW2_FLOW) {
            tl_error("Gajanan gkc - XDMF Sediment output still needs to be developed!");
            // /* prints the time step */
            // ps_print_ts_xdmf(hdf5, "Displacement", displacement, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            // ps_print_ts_xdmf(hdf5, "Bed_Shear_Stress", bed_shear_stress, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            // ps_print_ts_xdmf(hdf5, "Active_Layer_Thickness", active_layer_thickness, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            // for (ibl = 0; ibl < number_bed_layers; ibl++) {
            //   snprintf(string, MAXLINE, "BLT%d", itrns);
            //   ps_print_ts_xdmf(hdf5, string, bed_layer_thickness[ibl], grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            //   snprintf(string, MAXLINE, "BLD%d", itrns);
            //   ps_print_ts_xdmf(hdf5, string, bed_layer_distribution[ibl], grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            //   snprintf(string, MAXLINE, "BLBD%d", itrns);
            //   ps_print_ts_xdmf(hdf5, string, bed_layer_bulk_density[ibl], grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            //   snprintf(string, MAXLINE, "BLCES%d", itrns);
            //   ps_print_ts_xdmf(hdf5, string, bed_layer_critical_erosion_shear[ibl], grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            //   snprintf(string, MAXLINE, "BLERC%d", itrns);
            //   ps_print_ts_xdmf(hdf5, string, bed_layer_erosion_rate_constant[ibl], grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            //   snprintf(string, MAXLINE, "BLERE%d", itrns);
            //   ps_print_ts_xdmf(hdf5, string, bed_layer_erosion_rate_exponent[ibl], grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            // }
            // ps_print_ts_xdmf(hdf5, "Active_Layer_Dist", active_layer_distribution, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            // ps_print_ts_xdmf(hdf5, "Sed_MR", sed_mr, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            // ps_print_ts_xdmf(hdf5, "Bedload_Vector", bedload_vector, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, VECTOR2D_DATA);
            // ps_print_ts_xdmf(hdf5, "Susload_Vector", susload_vector, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, VECTOR2D_DATA);
        }
        else if (mod->flag.SEDIMENT && mod->flag.SW3_FLOW) {
            tl_error("Gajanan gkc - XDMF Sediment output still needs to be developed!");
            // //cjt-- modified to write 2 d data on 3 d grid(all x - sections are the same)
            // double *buffer = NULL;        /* cjt -- created this to store 2d data on 3d grid */
            // buffer = (double *) tl_alloc(sizeof(double), grid->nnodes);
            // for (i = 0; i < grid->nnodes; i++) {
            //   buffer[i] = 0.0;
            //   if (mesh_info.bottom_nodes[i] == YES) {
            //     buffer[i] = bed_shear_stress[map_3d2bottom[i]];
            //   }
            // }
            // ps_print_ts_xdmf(hdf5,"Bed_Shear_Stress", buffer, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            // buffer = (double *) tl_free(sizeof(double), grid->nnodes, buffer);
        }
#endif
        tl_free(sizeof(double), grid->nnodes, (void *) conc_dim);
        tl_free(sizeof(double), grid->nnodes, (void *) pconc_dim);
    }
    
    /*************************************************************************************************************************************************/
    
    int inode;
    if (mod->flag.SW2_FLOW && (mod->flag.GW_FLOW==OFF)) {
        SSW_2D *sw = mod->sw->d2;
        //        /* allocate space to compute the index */
        //        sarray_init_int(sw->iarray, grid->nnodes);
        //        sarray_init_dbl(sw->darray, grid->nnodes);
        //        for (ie=0; ie<grid->nelems2d; ie++) {
        //            for (inode=0; inode<NDONTRI; inode++) {
        //                sw->darray[grid->elem2d[ie].nodes[inode]] += grid->elem_error[ie];
        //                sw->iarray[grid->elem2d[ie].nodes[inode]] += 1;
        //            }
        //        }
        //        for (inode=0; inode<grid->nnodes; inode++) {
        //            sw->darray[inode] /= (double)sw->iarray[inode];
        //        }
#ifdef _MESSG
        //        comm_update_double(sw->darray, 1, grid->smpi);
        //        comm_update_int(sw->iarray, 1, grid->smpi);
        comm_update_double(sw->error, 1, grid->smpi);
#endif
        /* Print error */
        ps_print_ts_xdmf(hdf5, "Error", grid->elem_error, grid->nelems2d, c_time, hdf5->iprint_step, ELEM_2D_CENTERED, SCALAR_DATA);
        //ps_print_ts_xdmf(hdf5, "Error", sw->darray, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        ps_print_ts_xdmf(hdf5, "Hydro Error", sw->error, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
    }
    else if (mod->flag.SW3_FLOW) {
        SSW_3D *sw = mod->sw->d3;
        //        /* allocate space to compute the index */
        //        sarray_init_int(sw->iarray, grid->nnodes);
        //        sarray_init_dbl(sw->darray, grid->nnodes);
        //        for (ie=0; ie<grid->nelems3d; ie++) {
        //            for (inode=0; inode<NDONTET; inode++) {
        //                sw->darray[grid->elem3d[ie].nodes[inode]] += grid->elem_error[ie];
        //                sw->iarray[grid->elem3d[ie].nodes[inode]] += 1;
        //            }
        //        }
        //        for (inode=0; inode<grid->nnodes; inode++) {
        //            sw->darray[inode] /= (double)sw->iarray[inode];
        //        }
#ifdef _MESSG
        //        comm_update_double(sw->darray, 1, grid->smpi);
        //        comm_update_int(sw->iarray, 1, grid->smpi);
        comm_update_double(sw->error, 1, grid->smpi);
#endif
        /* Print error */
        ps_print_ts_xdmf(hdf5, "Error", grid->elem_error, grid->nelems3d, c_time, hdf5->iprint_step, ELEM_3D_CENTERED, SCALAR_DATA);
        //ps_print_ts_xdmf(hdf5, "Error", sw->darray, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        ps_print_ts_xdmf(hdf5, "Hydro Error", sw->error, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
    }
    else if (mod->flag.NS3_FLOW) {
        SNS_3D *ns = mod->ns->d3;
#ifdef _MESSG
        comm_update_double(ns->error, 1, grid->smpi);
#endif
        /* Print error */
        ps_print_ts_xdmf(hdf5, "Error", grid->elem_error, grid->nelems3d, c_time, hdf5->iprint_step, ELEM_3D_CENTERED, SCALAR_DATA);
        //ps_print_ts_xdmf(hdf5, "Error", sw->darray, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        ps_print_ts_xdmf(hdf5, "Hydro Error", ns->error, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
    }
    if (mod->ntransport > 0) { /* Mostly same for both 2D and 3D SW models */
        SCON *con=NULL;
        con = mod->con;
        for (itrns = 0; itrns < mod->ntransport; itrns++) {
            if (con[itrns].type == SAL) {
                snprintf(string, MAXLINE, "Salt_Concentration_Error");
            }
            else if (con[itrns].type == TMP) {
                snprintf(string, MAXLINE, "Temperature_Distribution_Error");
            }
            else {
                snprintf(string, MAXLINE, "Concentration_Error %d", itrns+1);
            }
#ifdef _MESSG
            comm_update_double(con[itrns].error, 1, grid->smpi);
#endif
            ps_print_ts_xdmf(hdf5, string, con[itrns].error, grid->nnodes, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
        }
    }
    /* close HDF groups and file */
    status = H5Gclose(hdf5->group_node_map);
    status = H5Gclose(hdf5->group_mat_id);
    status = H5Gclose(hdf5->group_orig_nd);
    status = H5Gclose(hdf5->group_nodes);
    status = H5Gclose(hdf5->group_ghost_nodes);
    status = H5Gclose(hdf5->group_elements);
    status = H5Gclose(hdf5->group_data);
    status = H5Gclose(hdf5->group_mesh);
    status = H5Fclose(hdf5->h5_file_ptr);
    
    /* If the xmf file exists, then move rename it before writing the current file.
     In this way, if the program terminates prematurely, we still have a copy of the
     state of of the output up until this point */
    
    fp = io_fopen(hdf5->xmf_filename, "r", TRUE);
    
    if (fp != NULL) {
        fclose(fp);
        snprintf(string, MAXLINE - 1, "%s_back", hdf5->xmf_filename);
        rename(hdf5->xmf_filename, string);
    }
    /* Write the xmf file */
    
    xmlSaveFormatFileEnc(hdf5->xmf_filename, hdf5->doc, "UTF-8", 1);
    
    
    /*************************************************************************************************************************************************/
    /*************************************************************************************************************************************************/
    /*************************************************************************************************************************************************/
    /*************************************************************************************************************************************************/
    /*************************************************************************************************************************************************/
    /*************************************************************************************************************************************************/
    /*************************************************************************************************************************************************/
    /*************************************************************************************************************************************************/
    if (mod->flag.SW3_FLOW || mod->flag.GW_FLOW) {
        /**********************************************************************************************/
        /* Open the HDF file and groups */
        hdf5->h5_surf_file_ptr = H5Fopen(hdf5->h5_surf_filename, H5F_ACC_RDWR, H5P_DEFAULT);
        hdf5->group_surf_mesh = H5Gopen(hdf5->h5_surf_file_ptr, "Geometry", H5P_DEFAULT);
        hdf5->group_surf_data = H5Gopen(hdf5->h5_surf_file_ptr, "Data", H5P_DEFAULT);
        hdf5->group_surf_elements = H5Gopen(hdf5->group_surf_mesh, "elements", H5P_DEFAULT);
        hdf5->group_surf_nodes = H5Gopen(hdf5->group_surf_mesh, "nodes", H5P_DEFAULT);
        hdf5->group_surf_ghost_nodes = H5Gopen(hdf5->group_surf_mesh, "ghost_nodes", H5P_DEFAULT);
        //hdf5->group_surf_orig_nd = H5Gopen(hdf5->group_surf_mesh, "orig_node_number", H5P_DEFAULT);
        //hdf5->group_surf_mat_id = H5Gopen(hdf5->group_surf_mesh, "mat_id", H5P_DEFAULT);
        hdf5->group_surf_node_map = H5Gopen(hdf5->group_surf_mesh, "node_map", H5P_DEFAULT);
        /* Set XDMF Element Type */
        sprintf(elem_type, "Triangle");
        nelem=grid->nelems2d_sur;
        
        /* if grid was adapted, things are more complicated */
        if (mod->flag.ADAPTED_THE_GRID == YES || REF_FLAG_RECV == YES || UNREF_FLAG_RECV == YES) {
            ps_print_surf_geo_xdmf(mod);
        }
        
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
        
        if (grid->nnodes_sur > grid->my_nnodes_sur)
        {
            /* GHOST NODE <Set/> information */
            node1 = xmlNewNode(NULL, BAD_CAST "Set");
            xmlNewProp(node1, BAD_CAST "Ghost", BAD_CAST "1");
            xmlNewProp(node1, BAD_CAST "SetType", BAD_CAST "Node");
            sprintf(buff1, "%s:/Geometry/ghost_nodes/%d", hdf5->h5_surf_filename, hdf5->current_mesh);
            node2 = xmlNewChild(node1, NULL, BAD_CAST "DataItem", BAD_CAST buff1);
            xmlNewProp(node2, BAD_CAST "DataType", BAD_CAST "Int");
            sprintf(buff1, "%d", grid->nnodes_sur - grid->my_nnodes_sur);
            xmlNewProp(node2, BAD_CAST "Dimensions", BAD_CAST buff1);
            xmlNewProp(node2, BAD_CAST "Format", BAD_CAST "HDF");
        }
        xmlAddChild(hdf5->current_grid_node_surf, node1);
        /**********************************************************************************************/
        
        
        
        
        /**********************************************************************************************/
        if (mod->flag.SW3_FLOW) {
            SSW_3D *sw = mod->sw->d3;
            ID_LIST_ITEM *ptr;
            double *wse2d = (double *) tl_alloc(sizeof(double), grid->nnodes_sur);
            for (i=0; i<grid->nnodes_sur; i++){
                ptr = grid->vertical_list[i];
                int surfnode = ptr->id;
                wse2d[i] = sw->displacement[surfnode] + grid->node[surfnode].z;
            }
            ps_print_surf_ts_xdmf(hdf5,"Surface_Elevation", wse2d, grid->nnodes_sur, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            wse2d = (double *) tl_free(sizeof(double), grid->nnodes_sur, wse2d);
            ps_print_surf_ts_xdmf(hdf5, "Depth", sw->depth, grid->nnodes_sur, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
            ps_print_surf_ts_xdmf(hdf5, "Depth_Averaged_Velocity", sw->depth_avg_vel, grid->nnodes_sur, c_time, hdf5->iprint_step, NODE_CENTERED, VECTOR2D_DATA);
        }
        /**********************************************************************************************/
#ifdef _ADH_GROUNDWATER
        else if (mod->flag.GW_FLOW) {
#ifdef _DWGW_COUPLING
            if (mod->flag.SW2_FLOW) {
                int nodeID_3d;
                SSW_2D *sw = mod->sw->d2;
                double *tempscalar = (double *) tl_alloc(sizeof(double), grid->nnodes_sur);
                SVECT2D *tempvec2d = (SVECT2D *) tl_alloc(sizeof(SVECT2D), grid->nnodes_sur);
                
                /* Print surface elevation */
                for (i=0; i<grid->nnodes_sur; i++){
                    nodeID_3d = grid->nodeID_2d_to_3d_sur[i];
                    tempscalar[i] = sw->head[nodeID_3d] + grid->node[nodeID_3d].z;
                }
                ps_print_surf_ts_xdmf(hdf5,"Surface_Elevation", tempscalar, grid->nnodes_sur, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
                
                /* Print Depth */
                for (i=0; i<grid->nnodes_sur; i++){
                    nodeID_3d = grid->nodeID_2d_to_3d_sur[i];
                    tempscalar[i] = sw->head[nodeID_3d];
                    //printf("\nnodeID_3d = %2i, Depth[%2i] = %f",nodeID_3d, i, tempscalar[i]);
                }
                ps_print_surf_ts_xdmf(hdf5,"Depth", tempscalar, grid->nnodes_sur, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
                
                /* Print Velocity */
                for (i=0; i<grid->nnodes_sur; i++){
                    nodeID_3d = grid->nodeID_2d_to_3d_sur[i];
                    tempvec2d[i].x = sw->vel[nodeID_3d].x;
                    tempvec2d[i].y = sw->vel[nodeID_3d].y;
                }
                ps_print_surf_ts_xdmf(hdf5,"Velocity", tempvec2d, grid->nnodes_sur, c_time, hdf5->iprint_step, NODE_CENTERED, VECTOR2D_DATA);
                
                /* Print Density */
                if (mod->flag.BAROCLINIC){
                    for (i=0; i<grid->nnodes_sur; i++){
                        nodeID_3d = grid->nodeID_2d_to_3d_sur[i];
                        tempscalar[i] = (sw->density[nodeID_3d]+1.0)*mod->density;
                    }
                    ps_print_surf_ts_xdmf(hdf5,"Density", tempscalar, grid->nnodes_sur, c_time, hdf5->iprint_step, NODE_CENTERED, SCALAR_DATA);
                }
                tempscalar = (double *) tl_free(sizeof(double), grid->nnodes_sur, tempscalar);
                tempvec2d = (SVECT2D *) tl_free(sizeof(SVECT2D), grid->nnodes_sur, tempvec2d);
            }
#endif
        }
#endif
        
        /**********************************************************************************************/
        status = H5Gclose(hdf5->group_surf_node_map);
        //status = H5Gclose(hdf5->group_surf_mat_id);
        //status = H5Gclose(hdf5->group_surf_orig_nd);
        status = H5Gclose(hdf5->group_surf_ghost_nodes);
        status = H5Gclose(hdf5->group_surf_nodes);
        status = H5Gclose(hdf5->group_surf_elements);
        status = H5Gclose(hdf5->group_surf_data);
        status = H5Gclose(hdf5->group_surf_mesh);
        status = H5Fclose(hdf5->h5_surf_file_ptr);
        
        fp = io_fopen(hdf5->xmf_surf_filename, "r", TRUE);
        if (fp != NULL) {
            fclose(fp);
            snprintf(string, MAXLINE - 1, "%s_back", hdf5->xmf_surf_filename);
            rename(hdf5->xmf_surf_filename, string);
        }
        xmlSaveFormatFileEnc(hdf5->xmf_surf_filename, hdf5->doc_surf, "UTF-8", 1);
    }
    
    /*************************************************************************************************************************************************/
    
    /* increments the printstep */
    hdf5->iprint_step++;
    
    /* fflush(NULL); */
#endif
    return;
}
