/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  partition_main.c This file partitions the grid for HPC usage.   */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "global_header.h"

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Paritions a model grid domain for distributed computing.
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout]  mod        (SMODEL *) a pointer to an AdH model structure
 * @param[in]  flag           (int) a flag to indicate = 1 no division of processors necessary or = 0 processors get divided in superfile read
 *
 * \note CJT\:: This file only currently works for tets and triangles!
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static int DEBUG;

void partition_main(SMODEL *mod, int flag) {
    
    //assert(mod->grid->type == COLUMNAR);
    
    DEBUG = OFF;
    
    int ierr = 0;
    int wrong_node1,wrong_node2;
#ifdef _MESSG
    SGRID *g = mod->grid;
    int ntransport;
    SCON *con;
    int i, ii, inode, ie, j;
    int old_nnode = g->macro_nnodes;		/* the old number of nodes */
    int nnodes_sur_prev = g->nnodes_sur;
    int min_id=g->smpi->npes, min_id_wgt=0, my_count=0;
    char str[25];
    int time;
    double wave_stress;
    int new_column_flag=YES;
    MIDPT_LIST_ITEM *myptr;
    
    
    /* ssupermodel uses cleanup routine below but doesn't need to split the mesh */
    if (flag < 1) {
        /* calculates the partition (gives partition # for each local node) */
        if((g->smpi->partition_flag) == 1){
            partition_adpt_final(g);
        } else {
            //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
            partition_form_final(g);
            //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
        }
        /* transfers the nodal data */
        //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
        partition_transfer(mod);
        //tag(MPI_COMM_WORLD); MPI_Barrier(MPI_COMM_WORLD);
    }
    
    if (g->type == COLUMNAR) {
        for (i=0; i<g->num_midpts; i++)  {
            myptr = g->midpt_list[i];
            if (myptr->node1<0 || myptr->node2<0) {
                printf("partition main1 :: midpoint: %d \t top_node: %d \t bot_node: %d\n",i,myptr->node1,myptr->node2);
                tl_error("");
            }
        }
    }
    
    /* cleanup partition */
    partition_cleanup(g, g->smpi->myid);
    
    /* renumbers the mesh */
    node_renumber(mod, 1);
    
    if (g->nelems3d>0) elem3d_renumber(g);
    elem2d_renumber(g);
    elem1d_renumber(g);
    
    if (g->type == COLUMNAR) {
        for (i=0; i<g->num_midpts; i++)  {
            myptr = g->midpt_list[i];
            if (myptr->node1<0 || myptr->node2<0) {
                printf("partition main1 :: midpoint: %d \t top_node: %d \t bot_node: %d\n",i,myptr->node1,myptr->node2);
                tl_error("");
            }
        }
    }
    
    // cjt :: this must be done before columns are built
    for (ie=0; ie<g->nelems2d; ie++) {
        if (g->elem2d[ie].nnodes == NDONTRI) {
            g->elem2d[ie].nedges = 3;
            g->elem2d[ie].edges = g->nd_on_TriEdge;
            g->elem2d[ie].nnodes_quad = NDONTRIQUAD;
        } else {
            tl_error("quadrilaterals not supported!");
        }
    }
    
    for (ie=0; ie<g->nelems3d; ie++) {
        if (g->elem3d[ie].nnodes == NDONTET) {
            g->elem3d[ie].nedges = 6;
            g->elem3d[ie].edges = g->nd_on_TetEdge;
            g->elem3d[ie].nnodes_quad = NDONTETQUAD;
        } else {
            tl_error("prisms not supported!");
        }
    }
    
    if(g->type == COLUMNAR) {
        build_columns(g,0);
#ifdef _DEBUG
        if (DEBUG) {
            tl_check_all_pickets(__FILE__,__LINE__);
        }
#endif
        for (i=0; i<g->num_midpts; i++)  {
            myptr = g->midpt_list[i];
            if (myptr->node1<0 || myptr->node2<0) {
                printf("partition main1 :: midpoint: %d \t top_node: %d \t bot_node: %d\n",i,myptr->node1,myptr->node2);
                tl_error("");
            }
        }
        node_renumber_surface(mod);
    }  else {
        if (mod->grid->ndim == 3) {
            classify_2d_elements(&mod->grid);  // CJT Added to make sure GW/NS/etc surface arrays are renumbered
        }
    }
    
#ifdef _DEBUG
    if (DEBUG && g->type == COLUMNAR) {
        for (i=0; i<g->num_midpts; i++)  {
            myptr = g->midpt_list[i];
            if (myptr->node1<0 || myptr->node2<0) {
                printf("partition main1 :: midpoint: %d \t top_node: %d \t bot_node: %d\n",i,myptr->node1,myptr->node2);
                tl_error("ERROR :: midpoint calculation bad!");
            }
        }
    }
#endif
    
    /* updates the message keys for the repartitioned mesh */
    comm_set_keys(g);
    comm_update_GN(1, g);
    comm_update_snode(g);
#ifdef _ADH_GROUNDWATER
    if (mod->flag.GW_FLOW)  comm_update_gw(mod->sgw, mod->grid);
#endif
    if (mod->flag.SW2_FLOW) comm_update_sw2(mod->sw->d2, mod->grid);
    if (mod->flag.SW3_FLOW) {
        comm_update_sw3(mod->sw->d3, mod->grid);
        comm_update_sw3_surface(mod->sw->d3, mod->grid);
    }
#ifdef _SEDIMENT
    if(mod->sed != NULL) comm_update_sed(mod->sed, mod->grid, mod->nconti, mod->nsed, mod->nlayers);
#endif
    if(mod->ntransport > 0) comm_update_con(mod->con, mod->grid, mod->ntransport);
    
    /* recalculate norms might not need to do whole mesh but being safe */
    if (g->ndim == 2) {
        int mine;
        for (i = 0; i < g->nelems2d; i++) {
            assert(g->elem2d[i].nnodes == NDONTRI);
            g->elem2d[i].my_pe = 0;
            if ((g->node[g->elem2d[i].nodes[0]].resident_pe != g->node[g->elem2d[i].nodes[1]].resident_pe) &&
                (g->node[g->elem2d[i].nodes[0]].resident_pe != g->node[g->elem2d[i].nodes[2]].resident_pe) &&
                (g->node[g->elem2d[i].nodes[1]].resident_pe != g->node[g->elem2d[i].nodes[2]].resident_pe)){
                
                for (inode = 0; inode < 3; inode++){
                    if (g->node[g->elem2d[i].nodes[inode]].resident_pe > g->smpi->myid) g->elem2d[i].my_pe++;
                }
                
            }
            else{
                for (inode = 0; inode < 3; inode++){
                    if (g->node[g->elem2d[i].nodes[inode]].resident_pe == g->smpi->myid) g->elem2d[i].my_pe++;
                }
            }
            g->wd_flag[i] = 1;
            
            // get 2D element djacs, normals and basis function gradients (elemental constants)
            int inode = UNSET_INT, nnodes2d = g->elem2d[i].nnodes;
            SVECT nds[nnodes2d];
            for (inode=0; inode<nnodes2d; inode++) {
                nds[inode].x = g->node[ g->elem2d[i].nodes[inode] ].x;
                nds[inode].y = g->node[ g->elem2d[i].nodes[inode] ].y;
                nds[inode].z = g->node[ g->elem2d[i].nodes[inode] ].z; // initial displacement should be added here
            }
            
            if (nnodes2d == NDONTRI) {
                get_triangle_linear_djac_nrml_gradPhi(&(g->elem2d[i]), NULL, nds);
            } else {
                g->elem2d[i].nrml = get_elem2d_normals(nds);
            }
            
            if (g->elem2d[i].djac < SMALL6) {
                fprintf(stderr, "WARNING :: MYID %d file %s line %d  :: Improperly numbered triangle=%d || Nodes %d %d %d Jacobian is: %20.10f \n",g->smpi->myid, __FILE__,__LINE__,(i + 1), g->elem2d[i].nodes[0], g->elem2d[i].nodes[1], g->elem2d[i].nodes[2],g->elem2d[i].djac);
                
                snode_printScreen(g->node[g->elem2d[i].nodes[0]]);
                snode_printScreen(g->node[g->elem2d[i].nodes[1]]);
                snode_printScreen(g->node[g->elem2d[i].nodes[2]]);
                exit(0);
            }
            
            g->elem_error[i] = 0.;
            g->elem2d[i].nedges = 3;
            g->elem2d[i].edges = g->nd_on_TriEdge;
            g->elem2d[i].nnodes_quad = NDONTRIQUAD;
        }
        
        for (i = 0; i < g->nelems1d; i++) {
            get_elem1d_linear_djac_gradPhi(g, &(g->elem1d[i]));
        }
        
        //node2elem2d_map(g);
        
        /* these are for 3d, but they will keep us from if-checking constantly elsewhere */
        
        g->my_nnodes_bed = g->my_nnodes;
        g->my_nnodes_sur = g->my_nnodes;
        g->nnodes_bed = g->nnodes;
        g->nnodes_sur = g->nnodes;
        
        
        if(((g->smpi->partition_flag) == 0) || (flag>0))g->initial_nnodes_bed = g->nnodes;
        if(((g->smpi->partition_flag) == 0) || (flag>0))g->initial_nelems_bed = g->nelems2d;
        if(flag>0) mod->grid->nelems2d_old = mod->grid->nelems2d;
        
        if (((g->smpi->partition_flag) == 0) || (flag>0)) {
            for (i=0; i<MAX_NNODES_ON_ELEM2D; i++) {
                mod->sw->d2->elem_rhs_dacont_extra_terms[i] = (double *) tl_realloc(sizeof(double),g->nelems2d,g->nelems2d_old,mod->sw->d2->elem_rhs_dacont_extra_terms[i]);
            }
            mod->grid->nelems2d_old = mod->grid->nelems2d;
        }
    }
    if (g->ndim == 3) {
        for (i = 0; i < g->nelems2d; i++) {
            assert(g->elem2d[i].nnodes == NDONTRI);
            g->elem2d[i].id = i;
            g->elem2d[i].id_orig = i;
            
            // get 2D element djacs, normals and basis function gradients (elemental constants)
            int inode = UNSET_INT, nnodes2d = g->elem2d[i].nnodes;
            SVECT nds[nnodes2d];
            for (inode=0; inode<nnodes2d; inode++) {
                nds[inode].x = g->node[ g->elem2d[i].nodes[inode] ].x;
                nds[inode].y = g->node[ g->elem2d[i].nodes[inode] ].y;
                nds[inode].z = g->node[ g->elem2d[i].nodes[inode] ].z; // initial displacement should be added here
            }
            
            if (nnodes2d == NDONTRI) {
                get_triangle_linear_djac_nrml_gradPhi(&(g->elem2d[i]), NULL, nds);
            } else {
                g->elem2d[i].nrml = get_elem2d_normals(nds);
            }
            
            if (g->elem2d[i].bflag == 1) { // on bed
                g->elem2d[i].djac = fabs(g->elem2d[i].djac);
                g->elem2d[i].djac3d = fabs(g->elem2d[i].djac3d);
                g->elem2d[i].djac3d_fixed = fabs(g->elem2d[i].djac3d_fixed);
            }
            
            
            if (g->elem2d[i].djac3d < SMALL6 || g->elem2d[i].djac3d_fixed < SMALL6) {
                selem2d_printScreen(&g->elem2d[i]);
                fprintf(stderr, "Improperly numbered triangle number %d :: djac:  %20.10e\n", (i + 1),g->elem2d[i].djac);
                tl_error("ERROR: Improperly numbered triangle number\n");
            }
            
            g->elem2d[i].my_pe = 0;
            if ((g->node[g->elem2d[i].nodes[0]].resident_pe != g->node[g->elem2d[i].nodes[1]].resident_pe) &&
                (g->node[g->elem2d[i].nodes[0]].resident_pe != g->node[g->elem2d[i].nodes[2]].resident_pe) &&
                (g->node[g->elem2d[i].nodes[1]].resident_pe != g->node[g->elem2d[i].nodes[2]].resident_pe)){
                
                for (inode = 0; inode < 3; inode++){
                    if (g->node[g->elem2d[i].nodes[inode]].resident_pe > g->smpi->myid) g->elem2d[i].my_pe++;
                }
                
            }
            else{
                for (inode = 0; inode < 3; inode++){
                    if (g->node[g->elem2d[i].nodes[inode]].resident_pe == g->smpi->myid) g->elem2d[i].my_pe++;
                }
            }
            //selem2d_printScreen(grid->elem2d[i]); exit(-1);
            if (g->elem2d[i].djac < 0.0) {
                fprintf(stderr, "WARNING :: MYID %d file %s line %d  :: Improperly numbered triangle=%d || Nodes %d %d %d Jacobian is: %20.10f \n",g->smpi->myid, __FILE__,__LINE__,(i + 1), g->elem2d[i].nodes[0], g->elem2d[i].nodes[1], g->elem2d[i].nodes[2],g->elem2d[i].djac);
                exit(0);
            }
            
            g->elem2d[i].nedges = 3;
            g->elem2d[i].edges = g->nd_on_TriEdge;
            g->elem2d[i].nnodes_quad = NDONTRIQUAD;
        }
        
        for (i = 0; i < g->nelems3d; i++) {
            assert(g->elem3d[i].nnodes == NDONTET);
            g->elem3d[i].my_pe = 0;
            min_id=g->smpi->npes;
            min_id_wgt=0;
            if ((g->node[g->elem3d[i].nodes[0]].resident_pe != g->node[g->elem3d[i].nodes[1]].resident_pe) &&
                (g->node[g->elem3d[i].nodes[0]].resident_pe != g->node[g->elem3d[i].nodes[2]].resident_pe) &&
                (g->node[g->elem3d[i].nodes[0]].resident_pe != g->node[g->elem3d[i].nodes[3]].resident_pe) &&
                (g->node[g->elem3d[i].nodes[1]].resident_pe != g->node[g->elem3d[i].nodes[2]].resident_pe) &&
                (g->node[g->elem3d[i].nodes[1]].resident_pe != g->node[g->elem3d[i].nodes[3]].resident_pe) &&
                (g->node[g->elem3d[i].nodes[2]].resident_pe != g->node[g->elem3d[i].nodes[3]].resident_pe)){
                
                for (inode = 0; inode < 4; inode++){
                    if (g->node[g->elem3d[i].nodes[inode]].resident_pe > g->smpi->myid) g->elem3d[i].my_pe++;
                    if (g->node[g->elem3d[i].nodes[inode]].resident_pe < min_id) min_id = g->node[g->elem3d[i].nodes[inode]].resident_pe;
                }
                
                if (min_id != g->smpi->myid) g->elem3d[i].my_pe = 0;
                
            }
            else{
                for (inode = 0; inode < 4; inode++){
                    if (g->node[g->elem3d[i].nodes[inode]].resident_pe == g->smpi->myid) g->elem3d[i].my_pe++;
                    if (g->node[g->elem3d[i].nodes[inode]].resident_pe < min_id) min_id = g->node[g->elem3d[i].nodes[inode]].resident_pe;
                }
                for (inode = 0; inode < 4; inode++){
                    if (g->node[g->elem3d[i].nodes[inode]].resident_pe == min_id) min_id_wgt++;
                }
                if(g->elem3d[i].my_pe==2){
                    if((min_id_wgt==2) && (min_id != g->smpi->myid)){
                        g->elem3d[i].my_pe = 0;
                    }
                }
                
            }
            
            if (g->elem3d[i].nnodes == NDONTET) {
                SVECT nds[g->elem3d[i].nnodes];
                for (inode=0; inode<g->elem3d[i].nnodes; inode++) {
                    nds[inode].x = g->node[ g->elem3d[i].nodes[inode] ].x;
                    nds[inode].y = g->node[ g->elem3d[i].nodes[inode] ].y;
                    nds[inode].z = g->node[ g->elem3d[i].nodes[inode] ].z;  // initial displacement should be added here
                }
                
                get_tet_linear_djac_gradPhi(&(g->elem3d[i]), NULL, nds);
                if (g->elem3d[i].djac < 0.0) {
                    fprintf(stderr, "WARNING :: MYID %d file %s line %d  :: Improperly numbered triangle=%d || Nodes %d %d %d Jacobian is: %20.10f \n",g->smpi->myid, __FILE__,__LINE__,(i + 1), g->elem2d[i].nodes[0], g->elem2d[i].nodes[1], g->elem2d[i].nodes[2],g->elem2d[i].djac);
                    exit(0);
                }
                
                g->elem3d[i].nedges = 6;
                g->elem3d[i].edges = g->nd_on_TetEdge;
                g->elem3d[i].nnodes_quad = NDONTETQUAD;
                
            }
        }
        
        /* computes the number material types */
        for (i = 0; i < g->nelems3d; i++) {
            if (g->elem3d[i].mat >= g->nmat) {
                g->nmat = g->elem3d[i].mat + 1;
            }
        }
        
        // determine 3d elements col, bed and surface elements
        // note the original code, this is a loop over nelem3d, not columns, which leads to slightly different residuals
        if (g->type == COLUMNAR) {
            int ie3d = UNSET_INT, icol=UNSET_INT;
            ID_LIST_ITEM *ptr;
            for(icol=0; icol<g->ncolumns; icol++)  {
                ptr = g->column_list[icol];
                while(ptr->next != NULL) {
                    ie3d = ptr->id;
                    g->elem3d[ie3d].icol = icol;
                    g->elem3d[ie3d].elem2d_sur = g->elem2d_sur[icol];
                    g->elem3d[ie3d].elem2d_bed = g->elem2d_bed[icol];
                    ptr = ptr->next;
                }
            }
        }
        
        if (debug.no_hydro == OFF) {
            if (((g->smpi->partition_flag) == 0) || (flag>0)) {
                g->initial_nnodes_bed = g->nnodes_bed;
                g->initial_nelems_bed = g->nelems2d_bed;
                if(mod->flag.SW3_FLOW) {
                    for (i=0; i<MAX_NNODES_ON_ELEM3D; i++) {
                        mod->sw->d3->elem_rhs_supg_dacont[i] = (double *) tl_realloc(sizeof(double),mod->grid->nelems3d,mod->grid->nelems3d_old,mod->sw->d3->elem_rhs_supg_dacont[i]);
                        mod->sw->d3->elem_rhs_supg_cont[i] = (double *) tl_realloc(sizeof(double),mod->grid->nelems3d,mod->grid->nelems3d_old,mod->sw->d3->elem_rhs_supg_cont[i]);
                    }
                }
                
#ifdef _ADH_GROUNDWATER
                if (mod->flag.GW_FLOW) {
                    mod->sgw->elem_3d_data = (ELEMENT_3D_DATA *) tl_realloc(sizeof(ELEMENT_3D_DATA),mod->grid->nelems3d,mod->grid->nelems3d_old,mod->sgw->elem_3d_data);
                    mod->sgw->elem_gw_flux = (SVECT *) tl_realloc(sizeof(SVECT),mod->grid->nelems3d,mod->grid->nelems3d_old,mod->sgw->elem_gw_flux);
                    
                    // cjt :: added to make sure when 3D GW grid is partitition, the surface grid reflects this.
                    for (i=0; i<MAX_NNODES_ON_ELEM2D; i++) {
                        mod->sw->d2->elem_rhs_dacont_extra_terms[i] = (double *) tl_realloc(sizeof(double),mod->grid->nelems2d,mod->grid->nelems2d_old,mod->sw->d2->elem_rhs_dacont_extra_terms[i]);
                    }
                }
#endif

                mod->grid->nelems3d_old = mod->grid->nelems3d;
            }
        }
        
    }
    
    if((g->smpi->partition_flag) == 0 || (flag>0))g->initial_nnodes = g->nnodes;
    if((g->smpi->partition_flag) == 0 || (flag>0))g->orig_initial_nnodes = g->nnodes;
    g->smpi->partition_flag = 1; //set flag for adaption
    g->nnodes_prev = g->nnodes;
    g->nnodes_sur_prev = g->nnodes_sur;
    g->nnodes_bed_prev = g->nnodes_bed;
    if(ierr > 0) exit(-1);
    
    /*update Wind and Wave Series */
    if (mod->flag.WIND) {
        if (mod->flag.WIND_STATION) {
            SSERIES *series = mod->series_wind_head;
            while(series != NULL) {
                
                smeteor_station_realloc(&(series->station), g->nnodes_sur, series->nnodes);
                series->nnodes = g->nnodes_sur;
                series = series->next;
            }
            mod->series_wind_head->nnodes = g->nnodes_sur;
            sseries_set_meteor_stations(mod->series_wind_head, mod->grid, WIND_SERIES);
        }
    }
    if (mod->flag.WAVE) {
        if (mod->flag.WAVE_STATION) {
            SSERIES *series = mod->series_wave_head;
            while(series != NULL) {
                
                smeteor_station_realloc(&(series->station), g->nnodes_sur, series->nnodes);
                series->nnodes = g->nnodes_sur;
                series = series->next;
            }
            if(mod->series_wave_head != NULL) mod->series_wave_head->nnodes = g->nnodes_sur;
            
            sseries_set_meteor_stations(mod->series_wave_head, mod->grid, WAVE_SERIES);
        }
    }
    
    
    // print newly partitioned grid info
    if (g->ndim == 2) {
        sgrid_printScreen(g, mod->io->geo2d.filename);
    } else if (g->ndim == 3) {
        sgrid_printScreen(g, mod->io->geo3d.filename);
    }
    
    
#endif
    
}
