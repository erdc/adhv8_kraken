#ifndef H_SIO_
#define H_SIO_

/*!
   \file type_io.h
   \brief Data Structures Needed by ADH_IO Information 
 */

#define UNSET_INT -3
#define REJECT " \t\v\f\n\r"
#define DELIMITERS " \t\n"
#define COMMENT_CHARACTERS "!#"

#define READ_SUCCESS  1
#define READ_NONE     0
#define READ_FAIL    -1

/* helper structs */
typedef struct {
  char filename[MAXLINE];       /* path of the file */
  FILE *fp;                     /* the file pointer */
} SFILE;

typedef struct {
  char filename[MAXLINE];       /* path of the file */
  FILE *fp;                     /* the file pointer */
  int  specified;               /* whether the file was specified in the super file */
} SFILE_IN;

/* main file information struct */
typedef struct {

  /* simulation information */
  char proj_name[MAXLINE];      /* The Project Name */
  char run_name[MAXLINE];       /* The Running Name */

  /* input files */
  SFILE sup;                 /* super file */
  SFILE_IN bc;               /* boundary condition file */
  SFILE_IN geo2d;            /* 2d geometry file */
  SFILE_IN geo3d;            /* 3d geometry file */
  SFILE_IN hot;              /* hot start file */
  SFILE_IN face;             /* file listing 2d face elements - should be included with 3d SW run */
  SFILE_IN sedfile;          /* sedlib bc input file */
  SFILE_IN nsmfile;          /* WQ NSM input file */

  SFILE_IN extrude_bin;           // for the 2d/3d extrusion utility
  SFILE_IN extrude_node;          // for the 2d/3d extrusion utility
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
// FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. //
#ifdef WINDLIB
  SFILE_IN windfile; /* Unlike other input files however, windfile = "fort.22" necessarily. */
#endif
// ABOVE LINES ADDED BY GAJANAN                                                                 //
//////////////////////////////////////////////////////////////////////////////////////////////////

  /* shallow water files */
  SFILE fout_sw_winds;
  SFILE fout_sw_waves;

  /* sw2d output files */
  SFILE fout_sw2_head;
  SFILE fout_sw2_old_head;
  SFILE fout_sw2_vel;
  SFILE fout_sw2_old_vel;
  SFILE fout_sw2_error;
  SFILE fout_sw2_error_hydro;

  /* sw3d output files */
  SFILE fout_sw3_depth;   /* a 2d grid file */
  SFILE fout_sw3_displacement;
  SFILE fout_sw3_old_displacement;
  SFILE fout_sw3_grid_speed;
  SFILE fout_sw3_pressure;
  SFILE fout_sw3_vel;
  SFILE fout_sw3_old_vel;
  SFILE fout_sw3_error;
  SFILE fout_sw3_error_hydro;
  SFILE fout_sw3_depth_avg_vel;
  SFILE fout_sw3_surface_vel;
  SFILE fout_sw3_bottom_vel;
  SFILE fout_sw3_hyd_viscosity; /* GSAVANT */
  SFILE fout_sw3_trn_diffusivity; /* GSAVANT */

  /* ns wind and wave files */
  SFILE fout_ns_winds;
  SFILE fout_ns_waves;

  /* ns2d output files */
  SFILE fout_ns2_head;
  SFILE fout_ns2_old_head;
  SFILE fout_ns2_vel;
  SFILE fout_ns2_old_vel;
  SFILE fout_ns2_error;
  SFILE fout_ns2_error_hydro;

  /* ns3d output files */
  SFILE fout_ns3_depth;   /* a 2d grid file */
  SFILE fout_ns3_displacement;
  SFILE fout_ns3_old_displacement;
  SFILE fout_ns3_grid_speed;
  SFILE fout_ns3_pressure;
  SFILE fout_ns3_vel;
  SFILE fout_ns3_old_vel;
  SFILE fout_ns3_error;
  SFILE fout_ns3_error_hydro;
  SFILE fout_ns3_depth_avg_vel;
  SFILE fout_ns3_surface_vel;
  SFILE fout_ns3_bottom_vel;

  /* general transport output files */
  SFILE *fout_con;
  SFILE *fout_error_con;

  // for extrusion utility
  SFILE fout_grid;
  SFILE fout_bc;
  SFILE fout_hot;
  SFILE fout_faces;

#ifdef _SEDIMENT
  /* sediment transport output files */
  SFILE fout_bed;               /* the output file for all bed variables */
  SFILE fout_bed_flux;          /* the output file for bed load flux vector */
  SFILE fout_active_layer;      /* the output file for all active layer variables */
  SFILE *fout_bed_layer;        /* an array of nlayer output files for all bed layer variables */
  SFILE *fout_sl_grain;         /* an array of ngrain output files for all suspended load grain variables */
  SFILE *fout_bl_grain;         /* an array of ngrain output files for all bed load grain variables */
#endif

#ifdef _ADH_GROUNDWATER
  SFILE fout_gw_phead;            /* the output file for the pressure heads */
  SFILE fout_gw_thead;            /* the output file for the equivalent freshwater heads */
  SFILE fout_gw_density;          /* the output file for the water densities */
  SFILE fout_gw_sat;              /* the output file for the saturations */
  SFILE fout_gw_flx;              /* the output file for the nodal fluxes */
  SFILE fout_gw_vel;              /* the output file for the velocities */
  SFILE fout_gw_error;            /* the output file for the error */
#endif
  /* other stuff */
  char cur_line[MAXLINE];       /* file line being processed */
  char cur_filename[MAXLINE];   /* name of file being processed */
  double timing_io;             /* Total IO Timing */
  double timing_random;         /* Arbitrary Timing */

} SIO;


/************************************************************/
/* struct methods ----------------------------------------- */

void sfile_open(SFILE *, char *filebase, char *ext1, char *ext2, char *suffix, char *mode, int stopOnNotFind);
void sio_init(SIO *, const char *);
void sio_free(SIO *, int, int, int);
void sio_init_winds(SIO *io);
void sio_init_waves(SIO *io);
void sio_init_transport(SIO *, int);
void sio_init_sediment(SIO *, int, int);
void init_adh_file(SFILE *);
void init_adh_in_file(SFILE_IN *);

/************************************************************/
/************************************************************/


#endif
