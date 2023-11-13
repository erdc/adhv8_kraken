#include "global_header.h"

void print_smodel_to_file(SMODEL *, char *);

void print_series_to_open_file(FILE *, SSERIES, int);

void print_series_list_to_file(FILE *, int, SSERIES *);

void print_sw2_to_file(SSW_2D *, SGRID *, char *);

void print_sw3_to_file();
//void print_sw3_to_file(SMODEL *, char *);

void print_scon_to_file(SCON *, int, SGRID *, char *);

void print_matrix_to_file(double *, SPARSE_VECT *, int, int, char *);

void print_matrix_profile_to_file(Profile_Info *, char *);

void print_resid_to_file(double *, SGRID *, int, char *);

void print_dble_array_to_file(char *, double *, int);

void print_int_array_to_file(char *, int *, int);

void print_smodel_to_file(SMODEL *mod, char *description){

  SGRID *g = mod->grid; //alias
  FILE *fp;
  int i;
  char filename[MAXLINE];
#ifdef _MESSG
  int myid = g->smpi->myid;
#else
  int myid = 0;
#endif

  fp = io_fopen(build_filename2(filename, MAXLINE, description, "_model_", myid, ".txt", UNSET_INT), "w", TRUE);

    fprintf(fp, "\n****************************************************\n");
    if (mod->flag.SW2_FLOW) {
        fprintf(fp, "MODEL: 2D - SHALLOW WATER");
        if (mod->flag.SW2_TRANSPORT)
            fprintf(fp, " - with transport\n");
    }
    else if (mod->flag.SW3_FLOW) {
        fprintf(fp, "MODEL: 3D - SHALLOW WATER\n");
        if (mod->flag.SW3_TRANSPORT)
            fprintf(fp, " - with transport\n");
    }
    else
        fprintf(fp, "MODEL: UNRECOGNIZED\n");

    fprintf(fp, "\nModel Flags ----------------------------------------\n");
    fprintf(fp, "moving_grid: %d\n", mod->flag.MG);
    fprintf(fp, "adaption: %d\n", mod->flag.GRID_ADAPTION);
    fprintf(fp, "icm: %d\n", mod->flag.ICM);
    fprintf(fp, "nsm: %d\n", mod->flag.NSM);
    fprintf(fp, "wave: %d\n", mod->flag.WAVE);
    fprintf(fp, "wind: %d\n", mod->flag.WIND);
    fprintf(fp, "coriolis: %d\n", mod->flag.CORIOLIS);
    fprintf(fp, "muc: %d \n", mod->flag.MUC);
    fprintf(fp, "conveyance: %d \n", mod->flag.CONVEYANCE);
    fprintf(fp, "ice: %d \n", mod->flag.ICE);
    fprintf(fp, "tidal_flag: %d\n", mod->flag.TIDE);
    fprintf(fp, "o_flag: %d \n", mod->o_flag);

    fprintf(fp, "\nSolver Variables  -----------------------------------\n");
    fprintf(fp, "inc_nonlin; %f\n", mod->inc_nonlin);
    fprintf(fp, "tol_nonlin: %f\n", mod->tol_nonlin);
    fprintf(fp, "max_nonlin_it: %d\n", mod->max_nonlin_it);
    fprintf(fp, "nalloc_inc: %d\n", mod->nalloc_inc);
    fprintf(fp, "nblock: %d\n", mod->nblock);

    fprintf(fp, "\nWetting/Drying Variables  ----------------------------\n");
    fprintf(fp, "drying_lower_limit: %f\n", mod->drying_lower_limit);
    fprintf(fp, "drying_upper_limit: %f\n", mod->drying_upper_limit);
    fprintf(fp, "wd_lower_tol: %f\n", mod->wd_lower_tol);
    fprintf(fp, "wd_upper_tol: %f\n", mod->wd_upper_tol);
    fprintf(fp, "wd_rate_lower_tol: %f\n", mod->wd_rate_lower_tol);
    fprintf(fp, "wd_rate_upper_tol: %f\n", mod->wd_rate_upper_tol);

    fprintf(fp, "\nModel Time Variables  -------------------------------\n");
    fprintf(fp, "t_init: %f\n", mod->t_init);
    fprintf(fp, "t_prev: %f\n", mod->t_prev);
    fprintf(fp, "t_final: %f\n", mod->t_final);
    fprintf(fp, "tau_temporal: %f\n", mod->tau_temporal);
    fprintf(fp, "out_level: %d\n", mod->out_level);

    fprintf(fp, "\nTime-Series: %d  -------------------------------------\n", mod->nseries);
    fprintf(fp, "\nTime-Series: %d  -------------------------------------\n", mod->nseries);
    print_series_to_open_file(fp, *(mod->series_dt), 1);
    print_series_to_open_file(fp, *(mod->series_out), 1);
    if (mod->series_wind_head != NULL) print_series_list_to_file(fp, 1, mod->series_wind_head);
    if (mod->series_wave_head != NULL) print_series_list_to_file(fp, 1, mod->series_wave_head);
    if (mod->series_head != NULL) sseries_printScreen_list(1, mod->series_head);

    fprintf(fp, "\nStrings  ---------------------------------------------\n");
    fprintf(fp, "nstring: %d\n", mod->nstring);

    fclose(fp);
    return;

}

/*************************************************************************************/

void print_series_to_open_file(FILE *fp, SSERIES ts, int full){
  int ientry;
  fprintf(fp, "------------------\n");
  fprintf(fp, "series: id: %d  :: type :: ", ts.id + 1);

  if(ts.type == WIND_SERIES) {
    fprintf(fp, "wind series\n");
  }
  if(ts.type == WAVE_SERIES) {
    fprintf(fp, "wave series\n");
  }
  else if(ts.type == TIME_SERIES) {
    fprintf(fp, "time series\n");
  }
  else if(ts.type == OUTPUT_SERIES) {
    fprintf(fp, "output series\n");
  }
  else if(ts.type == CONSTITUITIVE_SERIES) {
    fprintf(fp, "constituitive series\n");
  }
  else if(ts.type == ANY_SERIES) {
    fprintf(fp, "any series\n");
  }
  else if(ts.type == DT_SERIES) {
    fprintf(fp, "time-step series\n");
  }

  fprintf(fp, "size: %d  ", ts.size);
  fprintf(fp, "infact: %f  outfact: %f  ", ts.infact, ts.outfact);
  fprintf(fp, "nvalues: %d nsurface_nodes: %d\n",ts.nvalues,ts.nnodes);

  if(ts.type == WIND_SERIES) {
      fprintf(fp, "-- station location: {%f,%f} \n", ts.station->x, ts.station->y);
      fprintf(fp, "current time-interpolated value:  wind_stress.x: %10.5f \t wind_stress.y: %10.5f\n", ts.ivalue[0], ts.ivalue[1]);
  }
  else if(ts.type == WAVE_SERIES) {
      fprintf(fp, "-- station location: {%f,%f} \n", ts.station->x, ts.station->y);
      fprintf(fp, "current time-interpolated value:  wave_rads.xx: %10.5f \t wave_rads.xy: %10.5f \t wave_rads.yy: %10.5f\n",ts.ivalue[0], ts.ivalue[1], ts.ivalue[2]);
  }
  else {
    fprintf(fp, "current time-interpolated value: %10.5f\n", ts.ivalue[0]);
  }

  if(full) {
    for(ientry = 0; ientry < ts.size; ientry++) {
      fprintf(fp, "entry: %d \t", ientry);
      fprintf(fp, "time: %10.5f \t", ts.entry[ientry].time);

      if(ts.type == WIND_SERIES) {
          fprintf(fp, "wind_stress.x: %10.5f \t wind_stress.y: %10.5f", ts.entry[ientry].value[0], ts.entry[ientry].value[1]);
      }
      else if(ts.type == WAVE_SERIES) {
          fprintf(fp, "wave_rads.xx: %10.5f \t wave_rads.xy: %10.5f \t wave_rads.yy: %10.5f",ts.entry[ientry].value[0], ts.entry[ientry].value[1], ts.entry[ientry].value[2]);
      }
      else {
          fprintf(fp, "value: %10.5f", ts.entry[ientry].value[0]);
      }

      if(ts.type != OUTPUT_SERIES)
          fprintf(fp, "\t area: %10.5f \t slope: %10.5f\n", ts.entry[ientry].area[0], ts.entry[ientry].slope[0]);
      else
          fprintf(fp, "\n");
    }
  }
}

/***************************************************************************************/

void print_series_list_to_file(FILE *fp, int full, SSERIES *head) {

  SSERIES *ptr = head;
  while(ptr != NULL) {
    print_series_to_open_file(fp, *ptr, full);
    ptr = ptr->next;
  }
  return;

}

/**************************************************************************************/

void print_sw2_to_file(SSW_2D *sw, SGRID *g, char *description){

  FILE *fp;
  int i;
  char filename[MAXLINE];
#ifdef _MESSG
  int myid = g->smpi->myid;
#else
  int myid = 0;
#endif

  fp = io_fopen(build_filename2(filename, MAXLINE, description, "_sw2_", myid, ".txt", UNSET_INT), "w", TRUE);

  fprintf(fp, "i     head      old_head  older_head vel.x      vel.y      old_vel.x   old_vel.y   older_vel.x older_vel.y \n");

  for (i=0;i<g->nnodes;i++){
    fprintf(fp, "%d     %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e %10.5e \n", i, sw->head[i],sw->old_head[i], sw->vel[i].x, sw->vel[i].y, sw->old_vel[i].x, sw->old_vel[i].y, sw->older_vel[i].x, sw->older_vel[i].y);
  }

  fclose(fp);

  return;
}

/**************************************************************************************/

//void print_sw3_to_file(SMODEL *mod, char *description){
void print_sw3_to_file(){

/*printf("%s",description);
  FILE *fp;
  SSW_3D *sw=mod->sw->d3;
  SGRID *g=mod->grid;
  int i;
  char filename[MAXLINE];
  char desc[MAXLINE];
#ifdef _MESSG
  int myid = g->smpi->myid;
#else
  int myid = 0;
#endif
  int time;
 
  time=(mod->t_prev - mod->t_init)/mod->dt;
  sprintf(desc, "%s_%d", description,time);

  fp = io_fopen(build_filename2(filename, MAXLINE, desc, "_sw3_", myid, ".txt", UNSET_INT), "w", TRUE);

  fprintf(fp, "i gid disp old_disp older_disp grid_speed old_grid_speed desity prs prs_plus prs_minus error vert_node_flux velx vely velz old_velx old_vely old_velz older_velx older_vely older_velz tanvecx tanvecy hyd_visc trn_diff dpl_pert \n");

  for (i=0;i<g->nnodes;i++){
    fprintf(fp, "%d %d %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e\n", i, g->node[i].gid, sw->displacement[i],sw->old_displacement[i], sw->older_displacement[i], sw->grid_speed[i], sw->density[i], sw->prs[i], sw->prs_plus[i], sw->prs_minus[i], sw->error[i], sw->vertical_node_flux[i], sw->vel[i].x, sw->vel[i].y, sw->vel[i].z, sw->old_vel[i].x, sw->old_vel[i].y, sw->old_vel[i].z, sw->older_vel[i].x, sw->older_vel[i].y, sw->older_vel[i].z, sw->tanvec[i].x, sw->tanvec[i].y, sw->hyd_viscosity[i], sw->trn_diffusivity[i], sw->dpl_perturbation[i]);
  }

  fclose(fp);

  fp = io_fopen(build_filename2(filename, MAXLINE, description, "_sw3_2D_", myid, ".txt", UNSET_INT), "w", TRUE);

  fprintf(fp, "i gid depth old_depth depth_avg_vel old_depth_avg_vel surface_vel bottom_vel bed_elevation \n");

  for (i=0;i<g->nnodes_bed;i++){
    fprintf(fp, "%d %d %d %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e \n", i, g->node[g->nodeID_2d_to_3d_sur[i]].global_surf_id, g->node[g->nodeID_2d_to_3d_bed[i]].global_bed_id, sw->depth[i], sw->old_depth[i], sw->depth_avg_vel[i].x, sw->depth_avg_vel[i].y, sw->old_depth_avg_vel[i].x, sw->old_depth_avg_vel[i].y, sw->surface_vel[i].x, sw->surface_vel[i].y,sw->bottom_vel[i].x, sw->bottom_vel[i].y,sw->bed_elevation[i]);
  }

  fclose(fp);
*/
  return;
}

/*********************************************************************************/

void print_scon_to_file(SCON *con, int ntransport, SGRID *g, char *description){

  FILE *fp;
  int i, itrns;
  char filename[MAXLINE];
#ifdef _MESSG
  int myid = g->smpi->myid;
#else
  int myid = 0;
#endif

  fp = io_fopen(build_filename2(filename, MAXLINE, description, "_con_", myid, ".txt", UNSET_INT), "w", TRUE);
   
  for (itrns = 0; itrns < ntransport; itrns++) {
    fprintf(fp, "ITRNS #%d \n", itrns);
    for (i=0;i<4;i++) fprintf(fp, "property %d is %f \n", i, con[itrns].property[i]);
    fprintf(fp, "i gid con old_con older_con sink source ndl_decay error mfcf vcf.x vcf.y \n");
    for (i=0;i<g->nnodes;i++){
    fprintf(fp, "%d %d %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e %3.2e \n", i, g->node[i].gid, con[itrns].concentration[i], con[itrns].old_concentration[i], con[itrns].older_concentration[i], con[itrns].sink[i], con[itrns].source[i], con[itrns].nodal_decay_coefficient[i],  con[itrns].error[i], con[itrns].mfcf[i], con[itrns].vcf[i].x, con[itrns].vcf[i].y );
    }
  }

  fclose(fp);
  return;
}

/*****************************************************************************/

void print_matrix_to_file(double *diagonal, SPARSE_VECT * matrix, int nnode, int nsys_sq, char *description){

  FILE *fp;
  int i,j,k;
  char filename[MAXLINE];
#ifdef _MESSG
  int myid, ierr;
  ierr =  MPI_Comm_rank(cstorm_comm, &myid); // cjt :: only works for 1 grid in CSTORM :: MPI_COMM_WORLD,&myid); //only because solver routines don't get whole grid.    Revisit for multigrid
#else
  int myid = 0;
#endif

  fp = io_fopen(build_filename2(filename, MAXLINE, description, "_matrix_", myid, ".txt", UNSET_INT), "w", TRUE);
  
  fprintf(fp, "printing matrix: %s\n",description);
    for (i = 0; i < nnode; i++) {
        for (j = 0; j < nsys_sq; j++) {
            fprintf(fp, " i=%d, diagonal %30.20e \n", i + 1, diagonal[i * nsys_sq + j]);
        }
        fprintf(fp, "\n\n----------row %d--------------\n", (i + 1));
        for (j = 0; j < matrix[i].size; j++) {
            fprintf(fp, "**col %d**\n", (matrix[i].index[j] + 1));
            for (k = j * nsys_sq; k < (j + 1) * nsys_sq; k++)
                fprintf(fp, "%3d %30.20e\n", k - j * nsys_sq, matrix[i].value[k]);
        }
    }
    fclose(fp);
    return;
}

/*****************************************************************************/

void print_matrix_profile_to_file(Profile_Info *profile, char *description){

  FILE *fp;
  char filename[MAXLINE];
#ifdef _MESSG
  int myid, ierr;
  ierr =  MPI_Comm_rank(cstorm_comm, &myid); // cjt :: only works for 1 grid in CSTORM :: MPI_COMM_WORLD,&myid); //only because solver routines don't get whole grid.    Revisit for multigrid
#else
  int myid = 0;
#endif

  fp = io_fopen(build_filename2(filename, MAXLINE, description, "_matrix_profile_", myid, ".txt", UNSET_INT), "w", TRUE);

  int i=0, j=0;
  fprintf(fp, "printing matrix profile: %s \n",description);
  fprintf(fp, "profile size: %d \n",profile->size);
  for (i=0; i<profile->size; i++) {
    fprintf(fp, "rows :: begin: %d end %d size: %d \n",profile->rows[i].begin, profile->rows[i].end, profile->rows[i].size);
    for (j=0; j<profile->rows[i].size; j++) {
      fprintf(fp, "%20.10f \n",profile->rows[i].value[j]);
    }
    fprintf(fp, "cols :: begin: %d end %d size: %d \n",profile->cols[i].begin, profile->cols[i].end, profile->cols[i].size);
    for (j=0; j<profile->cols[i].size; j++) {
      fprintf(fp, "%20.10f \n",profile->cols[i].value[j]);
    }
 }
    fclose(fp);
    return;
}

/*********************************************************************************/

void print_resid_to_file(double *resid, SGRID *g, int nsys, char *description){

  FILE *fp;
  char filename[MAXLINE];
#ifdef _MESSG
  int myid = g->smpi->myid;
#else
  int myid = 0;
#endif

  fp = io_fopen(build_filename2(filename, MAXLINE, description, "_resid_", myid, ".txt", UNSET_INT), "w", TRUE);

    int i, j, k;
    fprintf(fp, "printing residual: %s \n", description);
    for (i = 0; i < g->nnodes; i++) {
        j = i * nsys;
        fprintf(fp, "i=%d, gid=%d final residual = ", i, g->node[i].gid );
        for (k = 0; k < nsys; k++)
            fprintf(fp, " %30.20e ", resid[j + k]);
        fprintf(fp, "\n");
    }
    fclose(fp);
    return;
}

/****************************************************************************************************/

void print_dble_array_to_file(char * description, double *array, int size) {
    int i;
  FILE *fp;
  char filename[MAXLINE];
#ifdef _MESSG
  int myid, ierr;
  ierr =  MPI_Comm_rank(cstorm_comm, &myid); // cjt :: only works for 1 grid in CSTORM :: MPI_COMM_WORLD,&myid); //only because solver routines don't get whole grid.    Revisit for multigrid
#else
  int myid = 0;
#endif

  fp = io_fopen(build_filename2(filename, MAXLINE, description, "", myid, ".txt", UNSET_INT), "w", TRUE);

    fprintf(fp, "printing array: %s sizee: %d \n",description, size);
    for (i=0; i<size; i++) {
        fprintf(fp, "[%d] \t %30.20f \n",i,array[i]);
    }
    fclose(fp);
    return;
}

/***************************************************************************************************/
void print_int_array_to_file(char * description, int *array, int size) {
    int i;
  FILE *fp;
  char filename[MAXLINE];
#ifdef _MESSG
  int myid, ierr;
  ierr =  MPI_Comm_rank(cstorm_comm, &myid); // cjt :: only works for 1 grid in CSTORM :: MPI_COMM_WORLD,&myid); //only because solver routines don't get whole grid.    Revisit for multigrid
#else
  int myid = 0;
#endif

  fp = io_fopen(build_filename2(filename, MAXLINE, description, "", myid, ".txt", UNSET_INT), "w", TRUE);

    fprintf(fp, "printing array: %s sizee: %d \n",description, size);
    for (i=0; i<size; i++) {
        fprintf(fp, "[%d] \t %20d \n",i,array[i]);
    }
    fclose(fp);
    return;
}

