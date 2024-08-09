
/* datatype for all model flags */

typedef struct {

  /* Model Flags */
  int NS_FLOW;              /* if TRUE then run Navier-Stokes */
  int NS2_FLOW;             /* if TRUE then run 2D Navier-Stokes */
  int NS3_FLOW;             /* if TRUE then run 3D Navier-Stokes */
  int GW_FLOW;              /* if TRUE then run ground water flow */
  int SW_FLOW;              /* if TRUE then run 2D or 3D shallow water */
  int SW2_FLOW;             /* if TRUE then run 2D shallow water */
  int SW3_FLOW;             /* if TRUE then run 3D shallow water */
  int DIFFUSIVE_WAVE;       /* if TRUE then run the 2D overland flow */

  /* moving grid */
  int MG;

  /* baroclinic :: flag to signify what variables affect density */
  /* baroclinic code = 0  :: no density effects */
  /* baroclinic code = 1  :: salinity only affects density */
  /* baroclinic code = 10 :: temperature only affects density */
  /* baroclinic code = 11 :: salinity and temperature affect density */
  int BAROCLINIC;

  int EOS; 
  /* 0 Linearlized EOS */
  /* 1 Full Equation EOS */

  /* mesh adaption */
  int GRID_ADAPTION;        /* flag that grid *can* be adapted during the simulation */
  int ADAPTED_THE_GRID;     /* flag that grid was adapted during a given time-step */
  int GRID_REFINED;         /* flat that grid was refined */
  int GRID_UNREFINED;       /* flag the grid was unrefined */

  /* time adaption */
  int TIME_ADAPT;           /* was t_adpt_flag */
  int TIME_ADAPT_FAIL;      /* was t_fail_flag */

  /* transport */
  int TRANSPORT;
  int NS2_TRANSPORT;
  int NS3_TRANSPORT;
  int SW2_TRANSPORT;
  int SW3_TRANSPORT;
  int VORTICITY;

  /* sediment */
  int SEDIMENT;
  int SEDLIB;

#ifdef _ADH_GROUNDWATER
  /* groundwater */
  /* taken from share_extern.h */
  int GW_SALINITY;
  int GW_TRANSPORT;
  int GW_REACTION;
  int RAY_TRACING;
  int SOCKETS;
  int METEOROLOGY;
#endif
  /* steady state flag */
  int STEADY_STATE;

  /* external libraries */
  int ICM;
  int NSM;

  /* surface stresses */
  int WAVE;
  int WIND;
  int WAVE_STATION;
  int WIND_STATION;
//////////////////////////////////////////////////////////////////////////////////////////////////
// FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. //
  int WIND_LIBRARY;   /* This flag is switched on by using "OP WNDLIB" */
  int NWS;
// ABOVE LINES ADDED BY GAJANAN                                                                 //
//////////////////////////////////////////////////////////////////////////////////////////////////

  /* coriolis force */
  int CORIOLIS;

  /* manning's unit flag */
  int MUC;

  /* units (MKE or FPS) */
  int UNITS;

  /* conveyance ?? */
  int CONVEYANCE;

  /* flag whether or not to printed the adapted mesh */
  int PRN_ADPT;

  /* flag whether to write output or not */
  int OUTPUT;

  /* flag for solver on initial automatic time determination */
  int SOLVE_ATF;

   /* if TRUE, singular matrix occurred, but don't exit */
  int UMFPACK_FAIL;

  /* ice */
  int ICE;              /* indicates that ice is included (SW) */
  int INS;              /* indicates that ice is designated by circular string */
  int nice_coords;      /* number of ice coordinates */

  /* tides */
  int TIDE;

  /* unrefine this ts */
  int UNREFINE;

  /* CStorm WSID flag */
  int CSTORM_WSID;

#ifdef _ADH_HDF5
  /* Parallel XDMF output */
  int PC_FILE_XDMF;
#endif

  int FLUX_WEIGHTED_NORMALS; /* Gajanan gkc adding. ON/OFF. Default OFF. Triggered using "NB FLXNML" in bc file. */

} SFLAGS;
