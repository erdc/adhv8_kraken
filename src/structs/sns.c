#include "global_header.h"

void sns_alloc_init(SNS **navier_stokes, SGRID *grid, SIO *io, SFLAGS flags) {
#ifdef _DEBUG
    printf("\n*** MYID %d NAVIER STOKES MODEL INITIALIZING ***\n", grid->smpi->myid);
#endif
    
    /* initialize shallow water model */
    (*navier_stokes) = (SNS *) tl_alloc(sizeof(SNS), 1);
    SNS *ns = (*navier_stokes);   // alias
    ns->d2 = NULL;
    ns->d3 = NULL;
    
    if (grid->ndim == 2) {
        sns_2d_alloc_init(&(ns->d2), grid, io, flags);
    }
    if (grid->ndim == 3) {
        sns_3d_alloc_init(&(ns->d3), grid, io, flags);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void sns_free(SNS *ns, SGRID *grid, SFLAGS flags) {
    assert(grid);
    
    /* free sw2d model if used */
    if (ns->d2 != NULL) {
        sns_2d_free(ns->d2, grid, flags);
    }
    
    /* free sw3d model if used */
    if (ns->d3 != NULL) {
        sns_3d_free(ns->d3, grid, flags);
    }
    
    /* free shallow water struct */
    ns = (SNS *) tl_free(sizeof(SNS), 1, ns);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void sns_open_output(SMODEL *mod) {
    
    assert(mod->io);
    int super = strlen(mod->io->sup.filename);    /* whether super file exists */
    
    /* open output for both NS2d and NS3d (head, waves, winds, etc.) */
    //open_output_file(&(mod->io->fout_ns2_head), "ns2d depth file", super);
    //assert(mod->io->fout_ns2_head.fp);
    //print_header(mod, mod->io->fout_ns2_head.fp, PS_FLAG_DEPTH);
    
    if (mod->flag.WIND && mod->file_output.wind == ON) {
        open_output_file(&(mod->io->fout_ns_winds), "sw wind file", super);
        assert(mod->io->fout_ns_winds.fp);
        print_header(mod, mod->io->fout_ns_winds.fp, PS_FLAG_WINDS);
    }
    
    if (mod->flag.WAVE && mod->file_output.wave == ON) {
        open_output_file(&(mod->io->fout_ns_waves), "sw wave file", super);
        assert(mod->io->fout_ns_waves.fp);
        print_header(mod, mod->io->fout_ns_waves.fp, PS_FLAG_WAVES);
    }
    
    /* open specific physics related output */
    if (mod->flag.NS2_FLOW) {
        sns_2d_open_output(mod);
    }
    else if (mod->flag.NS3_FLOW) {
        sns_3d_open_output(mod);
    }
    
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void sns_print_ts(SNS *ns, SIO *io, SGRID *grid, double time, int outfact, SFLAGS flag, SFILE sup, char *proj_name, int it1, int it2, int **ndata, int my_nnode_max, int *my_nnode_ext, int **ndata_sur, int my_nnode_max_sur, int *my_nnode_ext_sur, int flag1, SFILE_OUTPUT file_output) {
    /* convert time to specified output units */
    if (flag.NS2_FLOW) {
        if (ns->d2->winds != NULL && file_output.wind == ON) sns_print_ts_winds(grid, ns->d2->winds, io->fout_ns_winds.fp, time, ndata, my_nnode_max, my_nnode_ext);
        if (ns->d2->waves != NULL && file_output.wave == ON) sns_print_ts_waves(grid, ns->d2->waves, io->fout_ns_waves.fp, time, ndata, my_nnode_max, my_nnode_ext);
        sns_2d_print_ts(ns->d2, io, grid, time, outfact, sup, proj_name, it1, it2, ndata, my_nnode_max, my_nnode_ext, flag1, file_output);
    }
    else if (flag.NS3_FLOW) {
        if (ns->d3->winds != NULL && file_output.wind == ON) sns_print_ts_winds(grid, ns->d3->winds, io->fout_ns_winds.fp, time, ndata_sur, my_nnode_max_sur, my_nnode_ext_sur);
        if (ns->d3->waves != NULL && file_output.wave == ON) sns_print_ts_waves(grid, ns->d3->waves, io->fout_ns_waves.fp, time, ndata_sur, my_nnode_max_sur, my_nnode_ext_sur);
        sns_3d_print_ts(ns->d3, io, grid, time, outfact, sup, proj_name, it1, it2, ndata, my_nnode_max, my_nnode_ext, ndata_sur, my_nnode_max_sur, my_nnode_ext_sur, flag1, file_output);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void sns_print_ts_winds(SGRID *grid, SWIND *winds, FILE *fout, double time, int **ndata, int my_nnode_max, int *my_nnode_ext) {
    int inode = 0;
    FILE *fout1 = NULL;
    SVECT2D *temp = NULL;
    int nnodes = my_nnode_ext[grid->smpi->myid];
    temp = (SVECT2D *)tl_alloc(sizeof(SVECT2D), nnodes);
    if(grid->smpi->myid == 0) fprintf(fout, "TS 0 %15.8e\n", time);
    for(inode=0; inode < nnodes; inode++) {
        temp[inode].x = winds[inode].stress.x;
        temp[inode].y = winds[inode].stress.y;
    }
    svect2d_print_array_MPI(grid, fout, fout1, temp, ndata, my_nnode_max, my_nnode_ext, 0);
    temp = (SVECT2D *)tl_free(sizeof(SVECT2D), nnodes, temp);
    
}

void sns_print_ts_waves(SGRID *grid, SWAVE *waves, FILE *fout, double time, int **ndata, int my_nnode_max, int *my_nnode_ext) {
    int inode = 0;
    FILE *fout1 = NULL;
    SVECT2D *temp = NULL;
    int nnodes = my_nnode_ext[grid->smpi->myid];
    temp = (SVECT2D *)tl_alloc(sizeof(SVECT2D), nnodes);
    if(grid->smpi->myid == 0) fprintf(fout, "TS 0 %15.8e\n", time);
    for(inode=0; inode < nnodes; inode++) {
        temp[inode].x = waves[inode].stress.x;
        temp[inode].y = waves[inode].stress.y;
    }
    svect2d_print_array_MPI(grid, fout, fout1, temp, ndata, my_nnode_max, my_nnode_ext, 0);
    temp = (SVECT2D *)tl_free(sizeof(SVECT2D), nnodes, temp);
    
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void sns_read_hot(SMODEL *mod) {
    
    char line[MAXLINE];           /* the input line */
    char name[MAXLINE];           /* the name of the current data set */
    char msg[MAXLINE];            /* for reporting to screen */
    char *data = NULL;
    int i, itrns = 0, inode = 0, jnode = 0;     /* loop counter */
    double *tmp_data_set;         /* temporary space to store the current data set in */
    SGRID *grid;
    
    assert(mod);
    assert(mod->io);
    
    // these flag determine how to initialize arrays based on what is read in
    int hot_flag_con[mod->ntransport]; sarray_init_int(hot_flag_con, mod->ntransport);
    int hot_flag_pcon[mod->ntransport]; sarray_init_int(hot_flag_pcon, mod->ntransport);
    int flag_depth_found = FALSE, flag_old_depth_found = FALSE, flag_older_depth_found = FALSE;
    int flag_vel_found = FALSE, flag_old_vel_found = FALSE, flag_older_vel_found = FALSE;
    int flag_prs_found = FALSE;
    
    
    SIO io = *(mod->io);          // alias
    SNS_2D *ns2d = mod->ns->d2;   // alias
    SNS_3D *ns3d = mod->ns->d3;   // alias
    grid = mod->grid;
    int nnodes = grid->nnodes;
    int nnodes_surface = grid->nnodes_sur;
    int macro_nnodes = grid->macro_nnodes;
#ifdef _DEBUG
    printf("-- MYID %d reading shallow water hotstart file: %s\n",grid->smpi->myid,io.hot.filename);
#endif
    
    /* allocates space to read the current data set */
    tmp_data_set = (double *) tl_alloc(sizeof(double), 3 * macro_nnodes);
    
    /* loops over the lines in the input hotstart file looking for the file type */
    while ((fgets(line, MAXLINE, io.hot.fp) != NULL) && (strncmp(line, "DATASET", 7) != AGREE));
    
    /* reads the data sets while there are data sets */
    while (read_data_set(io, &io.hot, grid, tmp_data_set, line) == YES) {
        
        /* parse the name from the name line */
        io_save_line(&io, io.hot.fp, io.hot.filename, line);
        if (parse_card(line, &data) != CARD_NAME) {
            io_read_error(io, "Expected to read data set name.", TRUE);
        }
        read_text_field_custom(io, &data, name, MAXLINE, NULL, "data set name", 1, TRUE);
        
        
        /* convert name to uppercase to simplify comparisons */
        convert_to_uppercase(name);
        
        if (strncmp("IOH", name, 3) == AGREE) {
#ifdef _DEBUG
            root_print("------ reading initial depths.");
#endif
            flag_depth_found = TRUE;
            
            if (mod->flag.NS2_FLOW) {
                for (inode = 0; inode < nnodes; inode++)
                    ns2d->head[inode] = tmp_data_set[grid->node[inode].gid];
            }
            else if (mod->flag.NS3_FLOW) {
                for (inode = 0; inode < nnodes; inode++)
                    ns3d->depth[inode] = tmp_data_set[grid->node[inode].gid];
            }
        }
        
        else if (strncmp("IPH", name, 3) == AGREE) {
#ifdef _DEBUG
            root_print("------ reading previous depths.");
#endif
            
            flag_old_depth_found = TRUE;
            
            if (mod->flag.NS2_FLOW) {
                for (inode = 0; inode < nnodes; inode++)
                    ns2d->old_head[inode] = tmp_data_set[grid->node[inode].gid];
            }
            else if (mod->flag.NS3_FLOW) {
                for (inode = 0; inode < nnodes_surface; inode++)
                    ns3d->old_depth[inode] = tmp_data_set[grid->node[inode].gid];
            }
        }
        
        else if (strncmp("ID", name, 3) == AGREE) {
#ifdef _DEBUG
            root_print("------ reading initial displacements.");
#endif
            if (mod->flag.NS3_FLOW) {
                for (inode = 0; inode < nnodes; inode++) {
                    ns3d->displacement[inode] = tmp_data_set[grid->node[inode].gid];
                    ns3d->old_displacement[inode] = ns3d->displacement[inode];
                    ns3d->older_displacement[inode] = ns3d->displacement[inode];
                }
            }
        }
        
        else if (strncmp("IPD", name, 3) == AGREE) {
#ifdef _DEBUG
            root_print("------ reading initial old displacements.");
#endif
            if (mod->flag.NS3_FLOW) {
                for (inode = 0; inode < nnodes; inode++) {
                    ns3d->old_displacement[inode] = tmp_data_set[grid->node[inode].gid];
                }
            }
        }
        
        else if (strncmp("IPPD", name, 4) == AGREE) {
#ifdef _DEBUG
            root_print("------ reading initial older displacements.");
#endif
            if (mod->flag.NS3_FLOW) {
                for (inode = 0; inode < nnodes; inode++) {
                    ns3d->older_displacement[inode] = tmp_data_set[grid->node[inode].gid];
                }
            }
        }
        
        else if (strncmp("IV", name, 3) == AGREE) {
#ifdef _DEBUG
            root_print("------ reading initial velocities.");
#endif
            flag_vel_found = TRUE;
            
            if (mod->flag.NS2_FLOW) {
                for (inode = 0; inode < nnodes; inode++) {
                    jnode = 3 * grid->node[inode].gid;
                    ns2d->vel[inode].x = tmp_data_set[jnode];
                    ns2d->vel[inode].y = tmp_data_set[jnode + 1];
                }
            }
            else if (mod->flag.NS3_FLOW) {
                for (inode = 0; inode < nnodes; inode++) {
                    jnode = 3 * grid->node[inode].gid;
                    ns3d->vel[inode].x = tmp_data_set[jnode];
                    ns3d->vel[inode].y = tmp_data_set[jnode + 1];
                    ns3d->vel[inode].y = tmp_data_set[jnode + 2];
                }
            }
        }
        
        else if (strncmp("IPV", name, 3) == AGREE) {
#ifdef _DEBUG
            root_print("------ reading previous velocities.");
#endif
            flag_old_vel_found = TRUE;
            
            if (mod->flag.NS2_FLOW) {
                for (inode = 0; inode < nnodes; inode++) {
                    jnode = 3 * grid->node[inode].gid;
                    ns2d->old_vel[inode].x = tmp_data_set[jnode];
                    ns2d->old_vel[inode].y = tmp_data_set[jnode + 1];
                }
            }
            else if (mod->flag.NS3_FLOW) {
                for (inode = 0; inode < nnodes; inode++) {
                    jnode = 3 * grid->node[inode].gid;
                    ns3d->old_vel[inode].x = tmp_data_set[jnode];
                    ns3d->old_vel[inode].y = tmp_data_set[jnode + 1];
                    ns3d->old_vel[inode].y = tmp_data_set[jnode + 2];
                }
            }
        }
        
        else if (strncmp("IP", name, 2) == AGREE) {
#ifdef _DEBUG
            root_print("------ reading initial pressures.");
#endif
            flag_prs_found = TRUE;
            
            if (mod->flag.NS3_FLOW) {
                for (inode = 0; inode < nnodes; inode++) {
                    ns3d->prs[inode] = tmp_data_set[grid->node[inode].gid];
                }
            }
        }
        
        else if (strncmp("WAVE_STRESS", name, 11) == AGREE) {
#ifdef _DEBUG
            root_print(" Reading initial wave stresses.");
#endif
            if (mod->flag.WAVE != ON) {
                mod->flag.WAVE = ON;
                if (mod->flag.NS2_FLOW) swave_alloc(&mod->sw->d2->waves, nnodes);
                else if (mod->flag.NS3_FLOW) swave_alloc(&mod->sw->d3->waves, nnodes_surface);
            }
            for (inode = 0; inode < nnodes_surface; inode++) {
                if (mod->flag.NS2_FLOW) jnode = 3 * grid->node[inode].gid;
                else if (mod->flag.NS3_FLOW) jnode = 3 * grid->node[grid->nodeID_2d_to_3d_sur[inode]].global_surf_id;
                /* Removed units for StWave output (cjt) */
                if (mod->flag.NS2_FLOW) {
                    mod->sw->d2->waves[inode].rads.xx = tmp_data_set[jnode]; /* / 1000.0; */
                    mod->sw->d2->waves[inode].rads.yy = tmp_data_set[jnode + 1]; /* / 1000.0; */
                    mod->sw->d2->waves[inode].rads.xy = tmp_data_set[jnode + 2]; /* / 1000.0; */
                } else if (mod->flag.NS3_FLOW) {
                    mod->sw->d3->waves[inode].rads.xx = tmp_data_set[jnode]; /* / 1000.0; */
                    mod->sw->d3->waves[inode].rads.yy = tmp_data_set[jnode + 1]; /* / 1000.0; */
                    mod->sw->d3->waves[inode].rads.xy = tmp_data_set[jnode + 2]; /* / 1000.0; */
                    //printf("wave stress: %20.10f %20.10f %20.10f \n",mod->sw->d3->waves[inode].rads.xx,mod->sw->d3->waves[inode].rads.yy,mod->sw->d3->waves[inode].rads.xy);
                }
            }
        }
        
        
        /* read transport constituent initial conditions */
        else if (strncmp("ICON", name, 4) == AGREE) {
            
            /* get the constituent number from the name line */
            itrns = read_int_field_custom(io, &data, NULL, "constituent ID", 2, TRUE) - 1;
            if (itrns < 0 || itrns >= mod->ntransport) {
                io_read_error(io, "Tried to read concentration for non existent transport " "constituent.", TRUE);
            }
            hot_flag_con[itrns] = ON;
            
            if (mod->con[itrns].type == SLT || mod->con[itrns].type == SND || mod->con[itrns].type == CLA) {
#ifdef _DEBUG
                sprintf(msg, " Reading initital salt or sediment concentrations (constituent " "ID: %d).", itrns + 1);
                root_print(msg);
#endif
                /* if it is sediment it is read in in micromass per unit mass and we need it mass per unit mass */
                for (i = 0; i < nnodes; i++) {
                    mod->con[itrns].concentration[i] = 1.E-6 * tmp_data_set[grid->node[i].gid];
                }
            }
            else {
#ifdef _DEBUG
                sprintf(msg, "------ reading transport concentrations (constituent " "ID: %d).", itrns + 1);
                root_print(msg);
#endif
                for (i = 0; i < nnodes; i++) {
                    mod->con[itrns].concentration[i] = tmp_data_set[grid->node[i].gid];
                }
            }
        }
        
        else if (strncmp("IPC", name, 3) == AGREE) {
            
            /* get the constituent number from the name line */
            itrns = read_int_field_custom(io, &data, NULL, "constituent ID", 2, TRUE) - 1;
            if (itrns < 0 || itrns >= mod->ntransport) {
                io_read_error(io, "Tried to read concentration for non existent transport " "constituent.", TRUE);
            }
            
            hot_flag_pcon[itrns] = ON;
            
            if (mod->con[itrns].type == SLT || mod->con[itrns].type == SND || mod->con[itrns].type == CLA) {
#ifdef _DEBUG
                sprintf(msg, " Reading previous sediment concentrations (constituent " "ID: %d).", itrns + 1);
                root_print(msg);
#endif
                /* if it is sediment it is read in in micromass per unit mass and we need it mass per unit mass */
                for (i = 0; i < nnodes; i++) {
                    mod->con[itrns].old_concentration[i] = 1.E-6 * tmp_data_set[grid->node[i].gid];
                }
            }
            else {
#ifdef _DEBUG
                sprintf(msg, " Reading previous transport concentrations (constituent " "ID: %d).", itrns + 1);
                root_print(msg);
#endif
                for (i = 0; i < nnodes; i++) {
                    mod->con[itrns].old_concentration[i] = tmp_data_set[grid->node[i].gid];
                }
            }
        }
        
    }
    io_save_line(&io, NULL, "", "");
    rewind(io.hot.fp);
    
    /* force the initial condition to meet Dirichlet boundary conditions */
#ifdef _DEBUG
    root_print("------ (WARNING) Dirichlet boundary conditions are being applied to nodes (superseding initial conditions).");
#endif
    int istring = -1, isers = -1;
    for (i = 0; i < nnodes; i++) {
        istring = grid->node[i].string;
        if (istring > NORMAL) {
            if (mod->str_values[istring].ol_flow.bc_flag == BCT_PRS_DIR) {
                isers = mod->str_values[istring].ol_flow.iu_0;
                ns2d->head[i] = sseries_get_value(isers, mod->series_head,0);
            }
            else if (mod->str_values[istring].ol_flow.bc_flag == BCT_VEL_DIR) {
                isers = mod->str_values[istring].ol_flow.ivx;
                ns2d->vel[i].x = sseries_get_value(isers, mod->series_head,0);
                isers = mod->str_values[istring].ol_flow.ivy;
                ns2d->vel[i].y = sseries_get_value(isers, mod->series_head,0);
            }
            else if (mod->str_values[istring].ol_flow.bc_flag == BCT_VEL_PRS_DIR) {
                isers = mod->str_values[istring].ol_flow.ivx;
                ns2d->vel[i].x = sseries_get_value(isers, mod->series_head,0);
                isers = mod->str_values[istring].ol_flow.ivy;
                ns2d->vel[i].y = sseries_get_value(isers, mod->series_head,0);
                isers = mod->str_values[istring].ol_flow.iu_0;
                ns2d->head[i] = sseries_get_value(isers, mod->series_head,0);
            }
            else if (mod->str_values[istring].flow.bc_flag == BCT_VEL_DIR) {
                isers = mod->str_values[istring].flow.ivx;
                ns3d->vel[i].x = sseries_get_value(isers, mod->series_head,0);
                isers = mod->str_values[istring].flow.ivy;
                ns3d->vel[i].y = sseries_get_value(isers, mod->series_head,0);
                isers = mod->str_values[istring].flow.ivz;
                ns3d->vel[i].z = sseries_get_value(isers, mod->series_head,0);
                //printf("vel[%d]:  x=%20.10f y=%20.10f z=%20.10f\n",i,ns3d->vel[i].x,ns3d->vel[i].y,ns3d->vel[i].z);
            }
        }
    }
    
    
    
    //    if (flag_old_depth_found != TRUE) {
    //        // if old values are not given, assume that old and older are same as current
    //        if (mod->flag.NS2_FLOW) {
    //            sarray_copy_dbl(ns2d->old_head, ns2d->head, nnodes);
    //            sarray_copy_dbl(ns2d->older_head, ns2d->head, nnodes);
    //        } else if (mod->flag.NS3_FLOW) {
    //            sarray_copy_dbl(ns3d->old_depth, ns3d->depth, nnodes_surface);
    //        }
    //    } else {
    //        // if old values are given, assume older are the same
    //        if (mod->flag.NS2_FLOW) {
    //            sarray_copy_dbl(ns2d->older_head, ns2d->head, nnodes);
    //        } else if (mod->flag.NS3_FLOW) {
    //        }
    //    }
    
    if (flag_old_vel_found != TRUE) {
        // if old values are not given, assume that old and older are same as current
        if (mod->flag.NS2_FLOW) {
            svect2d_copy_array(ns2d->old_vel, ns2d->vel, nnodes);
            svect2d_copy_array(ns2d->older_vel, ns2d->vel, nnodes);
        } else if (mod->flag.NS3_FLOW) {
            svect_copy_array(ns3d->old_vel, ns3d->vel, nnodes);
            svect_copy_array(ns3d->older_vel, ns3d->vel, nnodes);
        }
    } else {
        // if old values are given, assume older are the same
        if (mod->flag.NS2_FLOW) {
            svect2d_copy_array(ns2d->older_vel, ns2d->vel, nnodes);
        } else if (mod->flag.NS3_FLOW) {
            svect_copy_array(ns3d->older_vel, ns3d->vel, nnodes);
        }
    }
    
    
    for (itrns=0; itrns<mod->ntransport; itrns++) {
        if (hot_flag_pcon[itrns] == OFF) {
            sarray_copy_dbl(mod->con[itrns].old_concentration, mod->con[itrns].concentration, nnodes);
            sarray_copy_dbl(mod->con[itrns].older_concentration, mod->con[itrns].old_concentration, nnodes);
        } else {
            sarray_copy_dbl(mod->con[itrns].older_concentration, mod->con[itrns].old_concentration, nnodes);
        }
    }
    
    /* clean up memory */
    tmp_data_set = (double *) tl_free(sizeof(double), 3 * macro_nnodes, tmp_data_set);
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
