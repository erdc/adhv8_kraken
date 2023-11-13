
#include "global_header.h"

#undef __FUNCT__
#define __FUNCT__ "adh_main_func"

// restructure notes ::
// nnodes2d_surf = nnodes_sur
// my_nnodes2d_surf = my_nnodes_sur
// node_local_map2D_to_3D = nodeID_2d_to_3d_sur
// node_local_map3D_to_2D = nodeID_3d_to_2d_sur

//*******************************************************//
//*******************************************************//
//*******************************************************//

/* The following routine will provide a mechanism to register the mesh
 ** with the mf2 coupler. This rouinte will be call from within the mf2
 ** coupler init routine.
 */

/*
 ** ne_g = number of global elements
 ** ne_l_max = max number of local elements among all procs. Probably don't need for ADH
 ** ne_l = number of local elements
 */

#ifdef _MESSG
void adh_export_func_(MPI_Fint *comm, int *tag2, int *wsid) {
    
    int i, id, nd, surface_nd, rank;
    double *tempadhx, *tempadhy;
    double *rnodecodesadh;
    double mintemp, maxtemp;
    double *buff, *d_ptr;
    int numsend;
    ID_LIST_ITEM *ptr;
    MPI_Comm c_comm;
    
    c_comm = MPI_Comm_f2c(*comm);
    MPI_Comm_rank(c_comm, &rank);
    
    // alias
    int myid = mod_cstorm->grid->smpi->myid;
    int npes = mod_cstorm->grid->smpi->npes;
    SFLAGS flag = mod_cstorm->flag;
    SGRID *grid = mod_cstorm->grid;
    SSW_2D *sw2d = mod_cstorm->sw->d2;
    SSW_3D *sw3d = mod_cstorm->sw->d3;
    
    
    switch (*wsid) {
            
        default:  // case 7
            
            if (flag.SW2_FLOW) {
                
                tempadhx = (double *)tl_alloc(sizeof(double), grid->my_nnodes);
                tempadhy = (double *)tl_alloc(sizeof(double), grid->my_nnodes);
                rnodecodesadh = (double *)tl_alloc(sizeof(double), grid->my_nnodes);
                
                /* calculate and write min and max surge sent to stwave */
                mintemp =  1e10;
                maxtemp = -1e10;
                for (i = 0; i < grid->my_nnodes; i++) {
                    if (sw2d->head[i] < 0.0) {
                        rnodecodesadh[i] = 0.0;
                    } else {
                        rnodecodesadh[i] = 1.0;
                    }
                    tempadhx[i] = sw2d->head[i] + grid->node[i].z;
                    if (tempadhx[i] < mintemp) {
                        mintemp = tempadhx[i];
                    }
                    if (tempadhx[i] > maxtemp) {
                        maxtemp = tempadhx[i];
                    }
                }
                printf("Minimum ADH -> StWave surge :: %lf \n", mintemp);
                printf("Maximum ADH -> StWave surge :: %lf \n", maxtemp);
                MPI_Send(rnodecodesadh, grid->my_nnodes, MPI_DOUBLE, 0, *tag2 - 1, c_comm);
                MPI_Send(tempadhx, grid->my_nnodes, MPI_DOUBLE, 0, *tag2, c_comm);
                printf("P%d: my_nnode=%d \n", rank, grid->my_nnodes);
                fflush(stdout);
                
                /* write data for checking */
                if (0) {
                    FILE *fp = NULL;
                    char name[MAXLINE];
                    int *iptr,k;
                    int it1, it2;
                    double tuse = mod_cstorm->t_prev + mod_cstorm->dt;
                    
                    it1 = (int)(tuse);
                    it2 = (int)(1000.0 * (tuse - (double)(it1)));
                    
                    fp = io_fopen(build_filename2(name, MAXLINE, "surge_expt", "-",it1, "-", myid), "w", TRUE);
                    
                    for (k = 0; k < grid->nnodes; k++) {
                        fprintf(fp," %d %lf %lf %lf \n", k, sw2d->head[k], grid->node[k].z, tempadhx[k]);
                    }
                    
                    fflush(fp);
                    fclose(fp);
                }
                
                /* send winds */
                for (i = 0; i < grid->my_nnodes; i++) {
                    tempadhx[i] = 0.;
                    tempadhy[i] = 0.;
                    if (flag.WIND) {
                        tempadhx[i] = sw2d->winds[i].stress.x;
                        tempadhy[i] = sw2d->winds[i].stress.y;
                    }
                }
                MPI_Send(tempadhx, grid->my_nnodes, MPI_DOUBLE, 0, *tag2 + 1, c_comm);
                MPI_Send(tempadhy, grid->my_nnodes, MPI_DOUBLE, 0, *tag2 + 2, c_comm);
                
                /* send velocity */
                for (i = 0; i < grid->my_nnodes; i++) {
                    tempadhx[i] = sw2d->vel[i].x;
                    tempadhy[i] = sw2d->vel[i].y;
                }
                MPI_Send(tempadhx, grid->my_nnodes, MPI_DOUBLE, 0, *tag2 + 3, c_comm);
                MPI_Send(tempadhy, grid->my_nnodes, MPI_DOUBLE, 0, *tag2 + 4, c_comm);
                
                /* free temp data */
                tempadhx = (double *)tl_free(sizeof(double), grid->my_nnodes, tempadhx);
                tempadhy = (double *)tl_free(sizeof(double), grid->my_nnodes, tempadhy);
                rnodecodesadh = (double *)tl_free(sizeof(double), grid->my_nnodes, rnodecodesadh);
                
            } else if (flag.SW3_FLOW) {
                
                tempadhx = (double *)tl_alloc(sizeof(double),  grid->my_nnodes_sur);
                tempadhy = (double *)tl_alloc(sizeof(double),  grid->my_nnodes_sur);
                rnodecodesadh = (double *)tl_alloc(sizeof(double),  grid->my_nnodes_sur);
                
                /* calculate and write min and max surge sent to stwave */
                mintemp =  1e10;
                maxtemp = -1e10;
                sarray_init_dbl(tempadhx, grid->my_nnodes_sur);
                sarray_init_dbl(tempadhy, grid->my_nnodes_sur);
                sarray_init_value_dbl(rnodecodesadh, grid->my_nnodes_sur,1.);
                for (i = 0; i < grid->nnodes_sur; i++) {
                    ptr = grid->vertical_list[i];
                    nd = ptr->id;
                    surface_nd = grid->nodeID_3d_to_2d_sur[nd];
                    assert(surface_nd != UNSET_INT);
                    if (grid->node[nd].resident_pe == grid->smpi->myid) {
                        tempadhx[surface_nd] = grid->node[nd].z + sw3d->displacement[nd];
                        if (tempadhx[surface_nd] < mintemp) {
                            mintemp = tempadhx[i];
                        }
                        if (tempadhx[surface_nd] > maxtemp) {
                            maxtemp = tempadhx[i];
                        }
                    }
                }
                
                printf("Minimum ADH -> StWave surge :: %lf \n", mintemp);
                printf("Maximum ADH -> StWave surge :: %lf \n", maxtemp);
                MPI_Send(rnodecodesadh,  grid->my_nnodes_sur, MPI_DOUBLE, 0, *tag2 - 1, c_comm);
                MPI_Send(tempadhx,  grid->my_nnodes_sur, MPI_DOUBLE, 0, *tag2, c_comm);
                printf("P%d: my_nnode=%d \n", rank, grid->my_nnodes_sur);
                fflush(stdout);
                
                /* send winds */
                sarray_init_dbl(tempadhx, grid->my_nnodes_sur);
                sarray_init_dbl(tempadhy, grid->my_nnodes_sur);
                for (i = 0; i < grid->nnodes_sur; i++) {
                    ptr = grid->vertical_list[i];
                    nd = ptr->id;
                    surface_nd = grid->nodeID_3d_to_2d_sur[nd];
                    if (grid->node[nd].resident_pe == grid->smpi->myid) {
                        if (flag.WIND) {
                            tempadhx[surface_nd] = sw3d->winds[surface_nd].stress.x;
                            tempadhy[surface_nd] = sw3d->winds[surface_nd].stress.y;
                        }
                    }
                }
                MPI_Send(tempadhx, grid->my_nnodes_sur, MPI_DOUBLE, 0, *tag2 + 1, c_comm);
                MPI_Send(tempadhy, grid->my_nnodes_sur, MPI_DOUBLE, 0, *tag2 + 2, c_comm);
                
                
                /* send velocity */
                sarray_init_dbl(tempadhx, grid->my_nnodes_sur);
                sarray_init_dbl(tempadhy, grid->my_nnodes_sur);
                for (i = 0; i < grid->nnodes_sur; i++) {
                    ptr = grid->vertical_list[i];
                    nd = ptr->id;
                    surface_nd = grid->nodeID_3d_to_2d_sur[nd];
                    if (grid->node[nd].resident_pe == grid->smpi->myid) {
                        tempadhx[surface_nd] = sw3d->vel[nd].x;
                        tempadhy[surface_nd] = sw3d->vel[nd].y;
                    }
                }
                MPI_Send(tempadhx, grid->my_nnodes_sur, MPI_DOUBLE, 0, *tag2 + 3, c_comm);
                MPI_Send(tempadhy, grid->my_nnodes_sur, MPI_DOUBLE, 0, *tag2 + 4, c_comm);
                
                /* free temp data */
                tempadhx = (double *)tl_free(sizeof(double), grid->my_nnodes_sur, tempadhx);
                tempadhy = (double *)tl_free(sizeof(double), grid->my_nnodes_sur, tempadhy);
                rnodecodesadh = (double *)tl_free(sizeof(double), grid->my_nnodes_sur, rnodecodesadh);
            }
            
            break;
        case 15:
            //            /* send data :
            //             1. rnodecodesadh
            //             2. surge
            //             3. winds: winds.x and winds.y
            //             4. velocity: vel.x and vel.y
            //             5. waves: waves.sxx, waves.sxy, waves.syy, waves.heigh, waves.period,
            //             waves.angle, waves.number/speed (?) waves.ediss_break
            //             */
            //            /* numsend = 15; */
            //            numsend = 6;
            //            buff = (double *)tl_alloc(sizeof(double), my_nnode * numsend);
            //            rnodecodesadh = (double *)tl_alloc(sizeof(double), my_nnode);
            //            tempadhx = (double *)tl_alloc(sizeof(double), my_nnode);
            //            tempadhy = (double *)tl_alloc(sizeof(double), my_nnode);
            //
            //
            //            /* calculate and write min and max surge sent to stwave */
            //            mintemp = 1e10;
            //            maxtemp = -1e10;
            //            for (i = 0; i < my_nnode; i++)
            //            {
            //                if (ol_head[i] < 0.0)
            //                {
            //                    rnodecodesadh[i] = 0.0;
            //                }
            //                else
            //                {
            //                    rnodecodesadh[i] = 1.0;
            //                }
            //                tempadhx[i] = ol_head[i] + node[i].z;
            //                if (tempadhx[i] < mintemp)
            //                {
            //                    mintemp = tempadhx[i];
            //                }
            //                if (tempadhx[i] > maxtemp)
            //                {
            //                    maxtemp = tempadhx[i];
            //                }
            //            }
            //            printf("Minimum ADH -> StWave surge :: %lf \n", mintemp);
            //            printf("Maximum ADH -> StWave surge :: %lf \n", maxtemp);
            //
            //            for (i = 0; i < my_nnode; i++)
            //            {
            //                tempadhx[i] = 0.0;
            //                tempadhy[i] = 0.0;
            //            }
            //
            //            /*MPI_Send(rnodecodesadh, my_nnode, MPI_DOUBLE, 0, *tag2 - 1, c_comm);*/
            //            /* pack data */
            //            for (i = 0, d_ptr = buff; i < my_nnode; i++, d_ptr += 6)
            //            {
            //                d_ptr[0] = rnodecodesadh[i];
            //                d_ptr[1] = ol_head[i] + node[i].z;
            //                d_ptr[2] = tempadhx[i];
            //                d_ptr[3] = tempadhy[i];
            //                d_ptr[4] = ol_vel[i].x;
            //                d_ptr[5] = ol_vel[i].y;
            //            }
            //            MPI_Send(buff, my_nnode*numsend, MPI_DOUBLE, 0, *tag2, c_comm);
            break;
    }
    return;
}

#else
void adh_export_func_(void){return;};
#endif

//*******************************************************//
//*******************************************************//
//*******************************************************//

#ifdef _MESSG
void adh_import_func_(MPI_Fint *comm, int *tag1, int *wsid) {
    
    int i, j, rank, numrecv;
    double *buff, *d_ptr;
    double *buff1;
    int ierr_code = MPI_ERR_UNKNOWN;
    MPI_Request  request;
    MPI_Status status;
    MPI_Comm c_comm;
    ID_LIST_ITEM *ptr;
    int nd, surface_nd, ipe;
    
    // alias
    int myid = mod_cstorm->grid->smpi->myid;
    int npes = mod_cstorm->grid->smpi->npes;
    SFLAGS *flag = &(mod_cstorm->flag);
    SGRID *grid = mod_cstorm->grid;
    SSW_2D *sw2d = mod_cstorm->sw->d2;
    SSW_3D *sw3d = mod_cstorm->sw->d3;
    
    c_comm = MPI_Comm_f2c(*comm);
    MPI_Comm_rank(c_comm, &rank);
    
    numrecv = 0;
    switch (*wsid) {
        case 7:
            numrecv = 3;
            flag->CSTORM_WSID = 7; // for adh writing
            
            if (flag->SW2_FLOW) {
                
                buff = (double *)tl_alloc(sizeof(double), numrecv * grid->my_nnodes);
                
                ierr_code = MPI_Irecv(buff, numrecv*grid->my_nnodes, MPI_DOUBLE, 0, *tag1, c_comm, &request);
                MPI_Wait(&request, &status);
                
                for (j = 0, d_ptr = buff; j < grid->my_nnodes; j++, d_ptr += 3) {
                    sw2d->waves[j].rads.xx = d_ptr[0];
                    sw2d->waves[j].rads.yy = d_ptr[1];
                    sw2d->waves[j].rads.xy = d_ptr[2];
                }
                
            } else if (flag->SW3_FLOW) {
                
                buff = (double *)tl_alloc(sizeof(double), numrecv * grid->my_nnodes_sur);
                ierr_code = MPI_Irecv(buff, numrecv*grid->my_nnodes_sur, MPI_DOUBLE, 0, *tag1, c_comm, &request);
                MPI_Wait(&request, &status);
                
                for (j = 0, d_ptr = buff; j < grid->nnodes_sur; j++, d_ptr += 3) {
                    ptr = grid->vertical_list[j];
                    nd = ptr->id;
                    surface_nd = grid->nodeID_3d_to_2d_sur[nd];
                    if (grid->node[nd].resident_pe == grid->smpi->myid) {
                        sw3d->waves[surface_nd].rads.xx = d_ptr[0];
                        sw3d->waves[surface_nd].rads.yy = d_ptr[1];
                        sw3d->waves[surface_nd].rads.xy = d_ptr[2];
                    }
                }
            }
            break;
            
        case 12:
            
            //            cstorm_wsid = 12; // for adh writing
            //
            //            numrecv = 8;
            //            buff = (double *)tl_alloc(sizeof(double), numrecv * my_nnode);
            //
            //            ierr_code = MPI_Irecv(buff, numrecv*my_nnode, MPI_DOUBLE, 0, *tag1, c_comm, &request);
            //            MPI_Wait(&request, &status);
            //
            //            for (j = 0, d_ptr = buff; j < my_nnode; j++, d_ptr += 8)
            //            {
            //                waves[j].Sxx = d_ptr[0];
            //                waves[j].Syy = d_ptr[1];
            //                waves[j].Sxy = d_ptr[2];
            //                waves[j].height = d_ptr[3];
            //                waves[j].period = d_ptr[4];
            //                waves[j].angle = d_ptr[5];
            //                waves[j].number = d_ptr[6];
            //                waves[j].ediss_break = -9810*d_ptr[7];
            //            }
            //            break;
            
        case 14:
            
            numrecv = 2;
            flag->CSTORM_WSID = 14; // for adh writing
            
            if (flag->SW2_FLOW) {
                
                buff = (double *)tl_alloc(sizeof(double), numrecv * grid->my_nnodes);
                
                ierr_code = MPI_Irecv(buff, numrecv*grid->my_nnodes, MPI_DOUBLE, 0, *tag1, c_comm, &request);
                MPI_Wait(&request, &status);
                
                for (j = 0, d_ptr = buff; j < grid->my_nnodes; j++, d_ptr += 2) {
                    sw2d->waves[j].stress.x = d_ptr[0];
                    sw2d->waves[j].stress.y = d_ptr[1];
                }
                
            } else if (flag->SW3_FLOW) {
                
                buff = (double *)tl_alloc(sizeof(double), numrecv * grid->my_nnodes_sur);
                ierr_code = MPI_Irecv(buff, numrecv*grid->my_nnodes_sur, MPI_DOUBLE, 0, *tag1, c_comm, &request);
                MPI_Wait(&request, &status);
                
                for (j = 0, d_ptr = buff; j < grid->my_nnodes_sur; j++, d_ptr += 2) {
                    ptr = grid->vertical_list[j];
                    nd = ptr->id;
                    surface_nd = grid->nodeID_3d_to_2d_sur[nd];
                    if (grid->node[nd].resident_pe == grid->smpi->myid) {
                        sw3d->waves[surface_nd].stress.x = d_ptr[0];
                        sw3d->waves[surface_nd].stress.y = d_ptr[1];
                    }
                }
            }
            break;
            
        default:
            break;
    }
    
    /* update ghost nodes */
    if (numrecv > 0) {
        if (flag->SW2_FLOW) {
            comm_update_swaves(grid, sw2d->waves, flag->CSTORM_WSID);
        } else if (flag->SW3_FLOW) {
            comm_update_swaves(grid, sw3d->waves, flag->CSTORM_WSID);
        }
    }
    
    
    if (1) { // print scatter plot of wave stresses and forces
        FILE *fp;
        messg_barrier(grid->smpi->ADH_COMM);
        char buffer[50];
        snprintf(buffer, 50, "cstorm_wave_force_%f.dat", mod_cstorm->t_prev);
        printf("writing to: %s\n",buffer);
        fp = io_fopen(buffer, "a", TRUE);
        messg_barrier(grid->smpi->ADH_COMM);;
        
        for (ipe=0; ipe<npes; ipe++) {
            if (myid != ipe) continue;
            printf("PE: %d writing wave forces\n",myid);
            
            if (flag->SW3_FLOW) {
                for (j=0; j<grid->nnodes_sur; j++) {
                    ptr = grid->vertical_list[j];
                    nd=ptr->id;
                    surface_nd = grid->nodeID_3d_to_2d_sur[nd];
                    fprintf(fp, "%20.10f %20.10f %20.10f %20.10f %20.10f\n",grid->node[nd].x, grid->node[nd].y, grid->node[nd].z, sw3d->waves[surface_nd].stress.x, sw3d->waves[surface_nd].stress.y);
                }
            } else {
                for (j=0; j<grid->nnodes; j++) {
                    fprintf(fp, "%20.10f %20.10f %20.10f %20.10f %20.10f\n",grid->node[j].x, grid->node[j].y, grid->node[j].z, sw2d->waves[j].stress.x, sw2d->waves[j].stress.y);
                }
            }
            messg_barrier(grid->smpi->ADH_COMM);;
        }
        fflush(fp);
        fclose(fp);
    }
    
    /* free buff */
    if (numrecv > 0) {
        if (flag->SW2_FLOW) {
            buff = (double *)tl_free(sizeof(double), numrecv * grid->my_nnodes, buff);
        } else if (flag->SW3_FLOW) {
            buff = (double *)tl_free(sizeof(double), numrecv * grid->my_nnodes_sur, buff);
        }
    }
    
    return;
}

#else
void adh_import_func_(void);
#endif

//*******************************************************//
//*******************************************************//
//*******************************************************//

#ifdef _MESSG
void adh_cpl_adhadpt_(MPI_Fint *comm, int *tag, int *ADPTCS_FLAG) {
    MPI_Comm c_comm;
    
    c_comm = MPI_Comm_f2c(*comm);
    /* Get the global number of nodes */
    /* need to take out next line ADPT_FLAG_CSTORM = 1; */
    *ADPTCS_FLAG = mod_cstorm->flag.ADAPTED_THE_GRID;
    
    if (mod_cstorm->grid->smpi->myid == 0)  { /* myid is 0 not 1; first argument needs to pass reference */
        MPI_Send(&mod_cstorm->flag.ADAPTED_THE_GRID, 1, MPI_INT, 0, *tag, c_comm);
        printf("P%d: ADH send ADPT_FLAG_STORM=%d to server\n", mod_cstorm->grid->smpi->myid, mod_cstorm->flag.ADAPTED_THE_GRID);
    }
    //ADPT_FLAG_CSTORM = 0;
    return;
}

#else
void adh_cpl_adhadpt_(void);
#endif

//*******************************************************//
//*******************************************************//
//*******************************************************//
// restructure notes ::
// nnodes2d_surf = nnodes_sur
// my_nnodes2d_surf = my_nnodes_sur
// node_local_map2D_to_3D = nodeID_2d_to_3d_sur
// node_local_map3D_to_2D = nodeID_3d_to_2d_sur
#ifdef _MESSG
void adh_cpl_init_(MPI_Fint *comm, int *tag, int *mnp) {
    
    /* nelem2d number of owned elements + ghost elements
     elem2d holds the element ie array with the node info at initial offset
     elem2d  - &elem2d.nodes, and spacing of sizeof(ELEM_2D);
     node is the x array that hold the xyz.  node[0].x,node[0].y, etc */
    /*
     ** This routine will do two things. First is to create a local to global
     ** table. This table will not map back to the original input deck, but will
     ** be relative, based on the current partition.
     ** Secondly, we will send the geometric and connectivity info to the driver.
     */
    
    MPI_Comm c_comm;
    int *owned;
    int buf[2], owner;
    int *ie = NULL, i, num_elems = 0, j, id, nd;
    int global_id, local_id;
    int *node_counts = NULL, *start_index = NULL;
    int *elm_local_to_global = NULL, start_el_global_num=0;
    int *node_local_to_global = NULL, start_node_global_num=0;
    double *geom_coords_x = NULL, *geom_coords_y = NULL;
    int my_nnode_local = 0;
    
    // alias
    int myid = mod_cstorm->grid->smpi->myid;
    int npes = mod_cstorm->grid->smpi->npes;
    SFLAGS flag = mod_cstorm->flag;
    SGRID *grid = mod_cstorm->grid;
    
    if (flag.SW3_FLOW) {
        my_nnode_local = grid->my_nnodes_sur;
    } else {
        my_nnode_local = grid->my_nnodes;
    }
    
    c_comm = MPI_Comm_f2c(*comm);
    buf[0] = grid->nelems2d;
    buf[1] = my_nnode_local;
    
    int nummaxelm = 50000;
    
    // sanity check
    if (messg_comm_rank(grid->smpi->ADH_COMM) != myid) printf("HOUSTON, WE HAVE A PROBLEM.\n");
    
    /* Get the global number of nodes */
    *mnp = my_nnode_local;
    
    int g_node_count;
    MPI_Allreduce(&my_nnode_local, &g_node_count, 1, MPI_INT, MPI_SUM, grid->smpi->ADH_COMM); // total nodes over all PEs
    
    if (flag.SW3_FLOW) assert(g_node_count == grid->macro_nnodes_sur);
    
    if (npes > 1) {
        MPI_Exscan(&buf[1], &start_node_global_num, 1, MPI_INT, MPI_SUM, grid->smpi->ADH_COMM); // get local starting ID in global node array
        //printf("myid: %d start_node_global_num: %d\n",myid,start_node_global_num);
    }
    
    if (messg_comm_rank(grid->smpi->ADH_COMM) == 0) {
        start_node_global_num = 0;
    }
    
    node_counts=(int *)tl_alloc(sizeof(int) , npes); // and local array of all my_nnodes on all PEs
    start_index=(int *)tl_alloc(sizeof(int) , npes);
    MPI_Allgather(&my_nnode_local, 1, MPI_INT, node_counts, 1, MPI_INT, grid->smpi->ADH_COMM);
    
    start_index[0] = 0;
    for (i = 1; i < npes; i++) {
        start_index[i] = start_index[i-1] + node_counts[i-1];
    }
    node_local_to_global = (int *)tl_alloc(sizeof(int) , my_nnode_local);
    
    for (i = 0; i < my_nnode_local; i++) {
        node_local_to_global[i] = ++start_node_global_num; // global node ID for local node #
    }
    
    /* Decide if I am going to be the owner of an element. Look at the ownership of all the nodes. The smallest rank wins */
    owned = (int *)tl_alloc(sizeof(int) , nummaxelm);
    num_elems = 0;
    for (i = 0; i < grid->nelems2d; i++) {
        
        if (flag.SW3_FLOW == ON) {
            if (grid->elem2d[i].bflag != 0) continue;  // not a free surface element
        }
        
        owner = npes;
        for (j = 0; j < 3; j++) {
            //owner = (node_pair[grid->elem2d[i].nodes[j]].sd < owner ? node_pair[grid->elem2d[i].nodes[j]].sd : owner)
            owner = (grid->node[grid->elem2d[i].nodes[j]].resident_pe < owner ? grid->node[grid->elem2d[i].nodes[j]].resident_pe : owner);
        }
        
        if (owner == myid) {
            owned[num_elems] = i;
            num_elems++;
            if ((num_elems % nummaxelm) == 0) {
                owned = realloc(owned, sizeof(int) * nummaxelm * (num_elems / nummaxelm));
            }
        }
    }
    
    if (npes > 1) {
        MPI_Exscan(&num_elems, &start_el_global_num, 1, MPI_INT, MPI_SUM, grid->smpi->ADH_COMM);
    }
    if (messg_comm_rank(grid->smpi->ADH_COMM) == 0) {
        start_el_global_num = 0;
    }
    elm_local_to_global = (int *)tl_alloc(sizeof(int) , num_elems);
    for (i = 0; i < num_elems; i++) {
        /* increment to convert to fotran indexing */
        elm_local_to_global[i] = ++start_el_global_num;
    }
    
    buf[0] = num_elems;
    buf[1] = my_nnode_local;
    MPI_Send(buf, 2, MPI_INT, 0, *tag, c_comm);
    /* Send  global element array */
    MPI_Send(node_local_to_global, my_nnode_local, MPI_INT, 0, *tag + 1, c_comm);
    /* Now we need to send local to global data to the MF_Server */
    /* Send global node number array */
    MPI_Send(elm_local_to_global, num_elems, MPI_INT, 0, *tag + 2, c_comm);
    /*free(elm_local_to_global);*/
    /*free(node_local_to_global);*/
    
    /* Now we need to send ie and xy coords to MF Coupler procs */
    /* Build an IE array */
    ie = (int *)tl_alloc(sizeof(int) , 3 * num_elems);
    
    // cjt :: initilize to catch problems
    int counter=0;
    for (i = 0; i < num_elems; i++) {
        for (j=0; j<3; j++) {
            ie[counter] = UNSET_INT;
            counter++;
        }
    }
    
    int *gnode_map3d_2d;    // a global 3d to surface 2d grid map
    int g_nnode = 0;        // max nnodes (ghost + residential) over all PE X npes
    int lmax_nnode = 0;     // the max # of local nodes (ghost + resid) on all PEs
    if (flag.SW3_FLOW == ON) {
        
        // doing this beause every PE needs all PE maps (cjt)
        MPI_Allreduce(&grid->nnodes, &lmax_nnode, 1, MPI_INT, MPI_MAX, grid->smpi->ADH_COMM); // get max nnodes (ghost + resid) over all PEs
        g_nnode = npes*lmax_nnode; // use max nnodes here for MPI_Allgather and straightforward access to maps
        //if (myid==0) printf("g_nnode: %d\n",g_nnode);
        gnode_map3d_2d  = (int *)tl_alloc(sizeof(int) , g_nnode);
        for (i=0; i<g_nnode; i++) {
            gnode_map3d_2d[i]  = UNSET_INT;
        }
        MPI_Allgather(grid->nodeID_3d_to_2d_sur, grid->nnodes, MPI_INT, gnode_map3d_2d, lmax_nnode, MPI_INT, grid->smpi->ADH_COMM);
        MPI_Barrier(grid->smpi->ADH_COMM);
        
        // test to make sure global maps are being built correctly
        /*
         FILE *fp1, *fp2;
         fp1 = io_fopen("gmap1.txt", "a", TRUE);
         fp2 = io_fopen("gmap2.txt", "w", TRUE);
         for (i=0; i<npes; i++) {
         if (myid==i) {
         printf("PE: %d writing\n",myid);
         for (j=0; j<nnode; j++) {
         fprintf(fp1,"%d\n",node_local_map3D_to_2D[j]);
         }
         for (j=nnode; j<lmax_nnode; j++){
         fprintf(fp1,"%d\n",UNSET_INT);
         }
         }
         fflush(fp1);
         MPI_Barrier(ADH_COMM);
         }
         fflush(fp1);
         fclose(fp1);
         MPI_Barrier(ADH_COMM);
         
         if (myid == 0) {
         for (i=0; i<g_nnode; i++) {
         fprintf(fp2,"%d\n",gnode_map3d_2d[i]);
         }
         }
         fflush(fp2);
         MPI_Barrier(ADH_COMM);
         fclose(fp2);
         */
    }
    
    
    for (i = 0; i < num_elems; i++) {   // loop over all cstorm PE owned elements
        for (j = 0; j < 3; j++) {
            if (flag.SW2_FLOW == ON) {
                
                local_id = grid->elem2d[owned[i]].nodes[j];		        // the local node number on this pe owned by this cstorm element
                owner = grid->node[local_id].resident_pe;                     // owner = node_pair[local_id].sd;  // the adh residential owner
                local_id = grid->node[grid->elem2d[owned[i]].nodes[j]].resident_id; // node_pair[local_id].rnode;  // the adh residential local node number onowning processor
                global_id = start_index[owner] + local_id;	            // the adh global id // Lucas says this is already calculated in AdH
                ie[(i * 3) + j] = global_id + 1;
                
            } else if (flag.SW3_FLOW == ON) {
                
                id = grid->elem2d[owned[i]].nodes[j];                    // local node id on 3d grid
                owner = grid->node[id].resident_pe;                            // owner = node_pair[id].sd;
                local_id = grid->node[grid->elem2d[owned[i]].nodes[j]].resident_id;  // local_id = node_pair[id].rnode;	// note, this local id may be on another PE
                local_id = gnode_map3d_2d[owner*lmax_nnode + local_id];  // local id here could be on another PE, so we need to access other PE maps
                global_id = start_index[owner] + local_id;
                ie[(i * 3) + j] = global_id + 1;                         // global id on 2d grid
                
            }
        }
    }
    
    
    /* Now generate the xy coords array to send to the mf server */
    geom_coords_x = (double *)tl_alloc(sizeof(double) , my_nnode_local);
    geom_coords_y = (double *)tl_alloc(sizeof(double) , my_nnode_local);
    
    if (flag.SW2_FLOW) {
        for (i = 0; i < my_nnode_local; i++) {
            geom_coords_x[i] = grid->node[i].x;
            geom_coords_y[i] = grid->node[i].y;
        }
    } else if (flag.SW3_FLOW) {
        for (i = 0; i < my_nnode_local; i++) {
            geom_coords_x[i] = grid->node[grid->nodeID_2d_to_3d_sur[i]].x;
            geom_coords_y[i] = grid->node[grid->nodeID_2d_to_3d_sur[i]].y;
        }
    }
    
    if (0) {// dt and t_prev not passed
        //        FILE *fp;
        //        char name[MAXLINE];
        //        int *iptr,k;
        //        int it1, it2;
        //        double tuse = mod_cstorm->t_prev + mod_cstorm->dt;
        //
        //        it1 = (int)(tuse);
        //        it2 = (int)(1000.0 * (tuse - (double)(it1)));
        //
        //        fp = io_fopen(build_filename2(name, MAXLINE, "nodes", "-", it1, "-", myid),"w", TRUE);
        //        for (i = 0; i < grid->my_nnodes; i++) {
        //            fprintf(fp," %d %d %lf %lf \n", i, node_local_to_global[i], geom_coords_x[i], geom_coords_y[i]);
        //        }
        //
        //        fflush(fp);
        //        fclose(fp);
    }
    
    MPI_Send(ie, num_elems * 3, MPI_INT, 0, *tag + 3, c_comm);
    MPI_Send(geom_coords_x, my_nnode_local, MPI_DOUBLE, 0, *tag + 4, c_comm);
    MPI_Send(geom_coords_y, my_nnode_local, MPI_DOUBLE, 0, *tag + 5, c_comm);
    
    
    ie = (int *)tl_free(sizeof(int) , 3 * num_elems, ie);
    owned = (int *)tl_free(sizeof(int) , nummaxelm, owned);
    node_counts=(int *)tl_free(sizeof(int) , npes, node_counts);
    start_index=(int *)tl_free(sizeof(int) , npes, start_index);
    geom_coords_x = (double *)tl_free(sizeof(double) , my_nnode_local, geom_coords_x);
    geom_coords_y = (double *)tl_free(sizeof(double) , my_nnode_local, geom_coords_y);
    elm_local_to_global = (int *)tl_free(sizeof(int) , num_elems, elm_local_to_global);
    node_local_to_global = (int *)tl_free(sizeof(int) , my_nnode_local, node_local_to_global);
    
    if (flag.SW3_FLOW) {
        gnode_map3d_2d  = (int *)tl_free(sizeof(int) , g_nnode, gnode_map3d_2d);
    }
    
    return;
}
#else
void adh_cpl_init_(void){return;}
#endif
