#include "adh.h"

/***************************************************/
/***************************************************/
/***************************************************/


void sfile_open(SFILE *file, char *filebase, char *ext1, char *ext2, char *suffix, char *mode, int stopOnNotFind) {
    
	char filename[MAXLINE];        /* the file name */
    init_adh_file(file);

    // puts the file name together 
    strcpy(filename, filebase);
    if (ext1 != NULL) strcat(filename, ext1);
	if (ext2 != NULL) strcat(filename, ext2);
	strcat(filename, suffix);
    fprintf(stdout,">> opening filename: %s\n",filename);

    // opens a file if it can and returns the pointer to the file 
    if((file->fp = fopen(filename, mode)) != NULL) return;

    // otherwise it writes an error message
    if (stopOnNotFind == YES) {
		fprintf(stderr, "FATAL ERROR :: Can't open file %s mode %s. \n", filename, mode);
    	exit(0);
	} else {
		fprintf(stdout, "WARNING :: file :: %s is not found.  Continuing on though.\n", filename);
	}
}

/* helper functions */
void init_adh_file(SFILE * file)
{
  strcpy(file->filename, "");
  file->fp = NULL;
}

void init_adh_in_file(SFILE_IN * file)
{
  strcpy(file->filename, "");
  file->fp = NULL;
  file->specified = FALSE;
}

//void init_adh_vec_files(ADH_VEC_FILES * vec) {
//  vec->files = NULL;
//  vec->nfiles = 0;
//}

/***************************************************/
/***************************************************/
/***************************************************/

void sio_init(SIO * io, const char *adh_root) {

  char line[MAXLINE];           /* the input line */
  char *data;                   /* the data after the first card is read */

  assert(io);

  /* simulation information */
  strcpy(io->proj_name, "");
  strcpy(io->run_name, "");

  strcpy(io->proj_name, adh_root);


  /* default names of standard input files usi.the project name */
  init_adh_file(&(io->sup));
  sprintf(io->sup.filename, "%s.sup", io->proj_name);

  init_adh_in_file(&(io->bc));
  sprintf(io->bc.filename, "%s.bc", io->proj_name);

  init_adh_in_file(&(io->geo2d));
  sprintf(io->geo2d.filename, "%s.2dm", io->proj_name);

  init_adh_in_file(&(io->geo3d));
  sprintf(io->geo3d.filename, "%s.3dm", io->proj_name);

  init_adh_in_file(&(io->hot));
  sprintf(io->hot.filename, "%s.hot", io->proj_name);

  init_adh_in_file(&(io->face));
  sprintf(io->face.filename, "%s.faces", io->proj_name);
    
  init_adh_in_file(&(io->extrude_bin));
  sprintf(io->extrude_bin.filename, "%s.bin", io->proj_name);
    
  init_adh_in_file(&(io->extrude_node));
  sprintf(io->extrude_node.filename, "%s.node", io->proj_name);

//////////////////////////////////////////////////////////////////////////////////////////////////
// FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. //
#ifdef WINDLIB
  /* initialize wind input file "fort.22". This name is necessary since wind library uses it. */
  init_adh_in_file(&(io->windfile));
  sprintf(io->windfile.filename,"fort.22");
#endif
// ABOVE LINES ADDED BY GAJANAN                                                                 //
//////////////////////////////////////////////////////////////////////////////////////////////////

  /* default names of standard output files */
  init_adh_file(&(io->fout_sw2_head));
  sprintf(io->fout_sw2_head.filename, "%s_dep.dat", io->proj_name);

  init_adh_file(&(io->fout_sw2_old_head));
  sprintf(io->fout_sw2_old_head.filename, "%s_old_dep.dat", io->proj_name);

  init_adh_file(&(io->fout_sw2_vel));
  sprintf(io->fout_sw2_vel.filename, "%s_vel.dat", io->proj_name);

  init_adh_file(&(io->fout_sw2_old_vel));
  sprintf(io->fout_sw2_old_vel.filename, "%s_old_vel.dat", io->proj_name);

  init_adh_file(&(io->fout_sw2_error));
  sprintf(io->fout_sw2_error.filename, "%s_error.dat", io->proj_name);

  init_adh_file(&(io->fout_sw2_error_hydro));
  sprintf(io->fout_sw2_error_hydro.filename, "%s_error_hydro.dat", io->proj_name);

  init_adh_file(&(io->fout_sw3_depth));
  sprintf(io->fout_sw3_depth.filename, "%s_dep.dat", io->proj_name);

  init_adh_file(&(io->fout_sw3_displacement));
  sprintf(io->fout_sw3_displacement.filename, "%s_dpl.dat", io->proj_name);

  init_adh_file(&(io->fout_sw3_old_displacement));
  sprintf(io->fout_sw3_old_displacement.filename, "%s_old_dpl.dat", io->proj_name);

  init_adh_file(&(io->fout_sw3_grid_speed));
  sprintf(io->fout_sw3_grid_speed.filename, "%s_grid_speed.dat", io->proj_name);

  init_adh_file(&(io->fout_sw3_pressure));
  sprintf(io->fout_sw3_pressure.filename, "%s_pressure.dat", io->proj_name);

  init_adh_file(&(io->fout_sw3_vel));
  sprintf(io->fout_sw3_vel.filename, "%s_vel.dat", io->proj_name);

  init_adh_file(&(io->fout_sw3_old_vel));
  sprintf(io->fout_sw3_old_vel.filename, "%s_old_vel.dat", io->proj_name);

  init_adh_file(&(io->fout_sw3_depth_avg_vel));
  sprintf(io->fout_sw3_depth_avg_vel.filename, "%s_depth_avg_vel.dat", io->proj_name);

  init_adh_file(&(io->fout_sw3_surface_vel));
  sprintf(io->fout_sw3_surface_vel.filename, "%s_surface_vel.dat", io->proj_name);

  init_adh_file(&(io->fout_sw3_bottom_vel));
  sprintf(io->fout_sw3_bottom_vel.filename, "%s_bottom_vel.dat", io->proj_name);

  init_adh_file(&(io->fout_sw3_error));
  sprintf(io->fout_sw3_error.filename, "%s_error.dat", io->proj_name);

  init_adh_file(&(io->fout_sw3_error_hydro));
  sprintf(io->fout_sw3_error_hydro.filename, "%s_error_hydro.dat", io->proj_name);

  init_adh_file(&(io->fout_sw3_hyd_viscosity));
  sprintf(io->fout_sw3_hyd_viscosity.filename, "%s_Viscosity.dat", io->proj_name); //GSAVANT

  init_adh_file(&(io->fout_sw3_trn_diffusivity));
  sprintf(io->fout_sw3_trn_diffusivity.filename, "%s_diffusivity.dat", io->proj_name); //GSAVANT

  // navier stokes
  init_adh_file(&(io->fout_ns3_depth));
  sprintf(io->fout_ns3_depth.filename, "%s_dep.dat", io->proj_name);

  init_adh_file(&(io->fout_ns3_displacement));
  sprintf(io->fout_ns3_displacement.filename, "%s_dpl.dat", io->proj_name);

  init_adh_file(&(io->fout_ns3_old_displacement));
  sprintf(io->fout_ns3_old_displacement.filename, "%s_old_dpl.dat", io->proj_name);

  init_adh_file(&(io->fout_ns3_grid_speed));
  sprintf(io->fout_ns3_grid_speed.filename, "%s_grid_speed.dat", io->proj_name);

  init_adh_file(&(io->fout_ns3_pressure));
  sprintf(io->fout_ns3_pressure.filename, "%s_pressure.dat", io->proj_name);

  init_adh_file(&(io->fout_ns3_vel));
  sprintf(io->fout_ns3_vel.filename, "%s_vel.dat", io->proj_name);

  init_adh_file(&(io->fout_ns3_old_vel));
  sprintf(io->fout_ns3_old_vel.filename, "%s_old_vel.dat", io->proj_name);

  init_adh_file(&(io->fout_ns3_depth_avg_vel));
  sprintf(io->fout_ns3_depth_avg_vel.filename, "%s_depth_avg_vel.dat", io->proj_name);

  init_adh_file(&(io->fout_ns3_surface_vel));
  sprintf(io->fout_ns3_surface_vel.filename, "%s_surface_vel.dat", io->proj_name);

  init_adh_file(&(io->fout_ns3_bottom_vel));
  sprintf(io->fout_ns3_bottom_vel.filename, "%s_bottom_vel.dat", io->proj_name);

  init_adh_file(&(io->fout_ns3_error));
  sprintf(io->fout_ns3_error.filename, "%s_error.dat", io->proj_name);

  init_adh_file(&(io->fout_ns3_error_hydro));
  sprintf(io->fout_ns3_error_hydro.filename, "%s_error_hydro.dat", io->proj_name);

  // extrusion output files 
  init_adh_file(&(io->fout_grid));
  sprintf(io->fout_grid.filename, "%s_3d.3dm", io->proj_name);

  init_adh_file(&(io->fout_hot));
  sprintf(io->fout_hot.filename, "%s_3d.hot", io->proj_name);

  init_adh_file(&(io->fout_bc));
  sprintf(io->fout_bc.filename, "%s_3d.bc", io->proj_name);

  init_adh_file(&(io->fout_faces));
  sprintf(io->fout_faces.filename, "%s_3d.faces", io->proj_name);

#ifdef _ADH_GROUNDWATER
  init_adh_file(&(io->fout_gw_phead));
  sprintf(io->fout_gw_phead.filename, "%s_phd.dat", io->proj_name);

  init_adh_file(&(io->fout_gw_thead));
  sprintf(io->fout_gw_thead.filename, "%s_thd.dat", io->proj_name);
  
  init_adh_file(&(io->fout_gw_density));
  sprintf(io->fout_gw_density.filename, "%s_den.dat", io->proj_name);

  init_adh_file(&(io->fout_gw_sat));
  sprintf(io->fout_gw_sat.filename, "%s_sat.dat", io->proj_name);

  init_adh_file(&(io->fout_gw_flx));
  sprintf(io->fout_gw_flx.filename, "%s_flx.dat", io->proj_name);

  init_adh_file(&(io->fout_gw_vel));
  sprintf(io->fout_gw_vel.filename, "%s_gw_vel.dat", io->proj_name);

  init_adh_file(&(io->fout_gw_error));
  sprintf(io->fout_gw_error.filename, "%s_error.dat", io->proj_name);

#endif
  /* check whether superfile exists (try to open it) */
  io->sup.fp = io_fopen(io->sup.filename, "r", FALSE);
  if (!(io->sup.fp)) {
    /* no superfile found so clear its filename and keep the defaults */
    sprintf(io->sup.filename, "");
  }
  else {
    /* read contents of the superfile */
    while (fgets(line, MAXLINE, io->sup.fp)) {
      io_save_line(io, io->sup.fp, io->sup.filename, line);
      if (strip_comments(line) <= 1) {
        continue;
      }
      /* read card */
        // CJT -- switch to string comparisons!
//      switch (parse_card(line, &data)) {
//        case CARD_BC:
//            /* read file path that can include spaces if in quotes */
//            read_text_field_custom(*io, &data, io->bc.filename, MAXLINE, NULL, "boundary condiion file name", 1, TRUE);
//            io->bc.specified = TRUE;
//            break;
//        case CARD_GEO:
//            read_text_field_custom(*io, &data, io->geo2d.filename, MAXLINE, NULL, "geometry file name", 1, TRUE);
//            io->geo2d.specified = TRUE;
//            break;
//        case CARD_GGEO:
//            read_text_field_custom(*io, &data, io->geo3d.filename, MAXLINE, NULL, "geometry file name", 1, TRUE);
//            io->geo3d.specified = TRUE;
//            break;
//        case CARD_HOT:
//            read_text_field_custom(*io, &data, io->hot.filename, MAXLINE, NULL, "hot start file name", 1, TRUE);
//            io->hot.specified = TRUE;
//            break;
//        case CARD_FACE:
//            read_text_field_custom(*io, &data, io->face.filename, MAXLINE, NULL, "boundary face file name", 1, TRUE);
//            io->face.specified = TRUE;
//            break;
//        default:
//            io_read_error(*io, "unrecognized superfile card.", FALSE);
//            break;
//      }
    }
  }
  if (io->sup.fp) {
    fclose(io->sup.fp);
  }
  io->sup.fp = NULL;
  io_save_line(io, NULL, "", "");

  return;

}

/************************************************************************/
/************************************************************************/
/************************************************************************/
void sio_init_transport(SIO *io, int ntransport) {
  int itrns = 0;
  io->fout_con = (SFILE *) tl_alloc(sizeof(SFILE), ntransport);
  io->fout_error_con = (SFILE *) tl_alloc(sizeof(SFILE), ntransport);

  for (itrns=0; itrns<ntransport; itrns++) {
    init_adh_file(&(io->fout_con[itrns]));
    build_filename2(io->fout_con[itrns].filename, MAXLINE, io->proj_name, "_con", itrns + 1, ".dat", UNSET_INT);
    build_filename2(io->fout_error_con[itrns].filename, MAXLINE, io->proj_name, "_error_con", itrns + 1, ".dat", UNSET_INT);
  }
}

/************************************************************************/
/************************************************************************/
/************************************************************************/
#ifdef _SEDIMENT
void sio_init_sediment(SIO *io, int nlayers, int nsed) {
    
    // initialize bed file
    init_adh_file(&(io->fout_bed));
    sprintf(io->fout_bed.filename, "%s_bed.dat", io->proj_name);
    
    // initialize bed load flux file
    init_adh_file(&(io->fout_bed_flux));
    sprintf(io->fout_bed_flux.filename, "%s_bed_flux.dat", io->proj_name);

    // initialize active layer file
    init_adh_file(&(io->fout_active_layer));
    sprintf(io->fout_active_layer.filename, "%s_active_layer.dat", io->proj_name);

    // allocate and initialize array of bed layer files
    io->fout_bed_layer = (SFILE *) tl_alloc(sizeof(SFILE), nlayers);
    int ilayer=0;
    for (ilayer=0; ilayer<nlayers; ilayer++) {
        init_adh_file(&(io->fout_bed_layer[ilayer]));
        build_filename2(io->fout_bed_layer[ilayer].filename, MAXLINE, io->proj_name, "_bed_layer_", ilayer + 1, ".dat", UNSET_INT);
    }

    // allocate and initialize array of grain files
    io->fout_sl_grain = (SFILE *) tl_alloc(sizeof(SFILE), nsed);
    io->fout_bl_grain = (SFILE *) tl_alloc(sizeof(SFILE), nsed);
    int ised = 0;
    for (ised=0; ised<nsed; ised++) {
        init_adh_file(&(io->fout_sl_grain[ised]));
        build_filename2(io->fout_sl_grain[ised].filename, MAXLINE, io->proj_name, "_sl_grain_", ised + 1, ".dat", UNSET_INT);

        init_adh_file(&(io->fout_bl_grain[ised]));
        build_filename2(io->fout_bl_grain[ised].filename, MAXLINE, io->proj_name, "_bl_grain_", ised + 1, ".dat", UNSET_INT);
    }

}
#endif

/************************************************************************/
/************************************************************************/
/************************************************************************/
void sio_init_winds(SIO *io) {
    init_adh_file(&(io->fout_sw_winds));
    sprintf(io->fout_sw_winds.filename, "%s_winds.dat", io->proj_name);
}


/************************************************************************/
/************************************************************************/
/************************************************************************/
void sio_init_waves(SIO *io) {
    init_adh_file(&(io->fout_sw_waves));
    sprintf(io->fout_sw_waves.filename, "%s_waves.dat", io->proj_name);
}

/************************************************************************/
/************************************************************************/
/************************************************************************/
void sio_free(SIO *io, int ntransport, int nlayers, int nsed) {

     if (ntransport>0) {
        io->fout_con = (SFILE *) tl_free(sizeof(SFILE), ntransport, io->fout_con);
        io->fout_error_con = (SFILE *) tl_free(sizeof(SFILE), ntransport, io->fout_error_con);
     }

#ifdef _SEDIMENT
     if (nlayers>0) {
        io->fout_bed_layer = (SFILE *) tl_free(sizeof(SFILE), nlayers, io->fout_bed_layer);
     }
     if (nsed>0) {
         io->fout_bl_grain = (SFILE *) tl_free(sizeof(SFILE), nsed, io->fout_bl_grain);
         io->fout_sl_grain = (SFILE *) tl_free(sizeof(SFILE), nsed, io->fout_sl_grain);
     }
#endif

     io = (SIO *) tl_free(sizeof(SIO), 1, io);
}

/************************************************************************/
/************************************************************************/
/************************************************************************/
