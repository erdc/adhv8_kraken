/* An AdH SMODEL is one of any number of models being ran
 * by AdH.  It contains the following:
 *    - one bc, hot and 2/3dm file
 *    - one type of PHYSICS (note sw2d/sw3d is one physics, SW)
 *    - one set of time-series, materials, strings, etc.
 */


#ifndef H_SMODEL_
#define H_SMODEL_

/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {
    
    int id;   /* the model ID for multi-model runs */
    int amICoupled; // monolithic test
    char name[MAXLINE];
    
    SIO *io;  /* file names associated with this model application */
    SFILE_OUTPUT file_output; /*moved here for multi-modle systems*/
    
    /* Shallow Water Physics */
    SSW *sw;
    SNS *ns;
    SCON *con;
    
#ifdef _SEDIMENT
    SSED *sed;
#endif

#ifdef _ADH_GROUNDWATER
    SGW *sgw;
#endif
  
    /* Grids */
    int ndim;                 /* the grid dimension */
    SGRID *grid;              /* the grid */
    
    SFLAGS flag;
    SMAT *mat;                /* materials */
    STR_VALUE *str_values;    /* strings */
    SWEIR_C *weir;            /* weirs*/
    SFLAP_C *flap;            /* flap gates */
    SSLUICE_C *sluice;        /* Sluice Gates */
    
    SSERIES *series_head;
    SSERIES *series_curr;
    SSERIES *series_dt;       /* time-step series */
    SSERIES *series_out;      /* output series */
    SSERIES *series_wind_head, *series_wind_curr;
    SSERIES *series_wave_head, *series_wave_curr;

#ifdef _ADH_GROUNDWATER
    SSERIES *series_gw_psk_head, *series_gw_psk_curr;
#endif
    
    int nseries;              /* the number of series in this model */
    int ntides;               /* the number of tidal series */
    int ntransport;           /* the number of transport consituents */
    int nstring;              /* the total number of user + runtime strings */
    int nstring_rt;           /* the number of runtime created strings (db ovl) */
    int nmat;                 /* the number of materials on mother grid */
    int nice_coords;          /* number of ice coordinates */
    int nweir;                /* numnber of weirs in the simulation */
    int nflap;                /*  number of flap gates in the simulation */
    int nsluice;              /* number of sluice gates in the simulation */
    
    /* ice */
    SVECT *ice_coords;     /* keeps track of center coordinates of ice flos */
    int *ice_nodes;       /* keeps track of the string ids associated with ice coordinates */
    int *ice_elems;       /* keeps track of the elems covered by ice */
    
    double viscosity;
    double manning_units_constant;
    double gravity;
    double density;
    
    int out_level;
    int o_flag;           /* some sort of output flag */
    
    /* wetting and drying variables */
    double drying_lower_limit;
    double drying_upper_limit;
    double wd_lower_tol;
    double wd_upper_tol;
    double wd_rate_lower_tol;
    double wd_rate_upper_tol;
    
    /* solver variables */
    Solver_Info solver_info;
    int max_nonlin_it;
    int nalloc_inc;
    int nblock;
    double inc_nonlin;
    double tol_nonlin;
    
    /* time-steping */
    double t_init;
    double t_prev;
    double t_final;
    double t_final_store; // cjt :: for use in screen write of cstorm progress
    double tau_temporal;
    double dt, old_dt, dt_err;
    int t_adpt_flag;
    int ientry_out;
    int green_ampt;   
    /* petrov coefficient */
    double tau_pg;
    
    /* load perturbution */
    double perturbation;
    
    /* boundary conditions mask */
    int *bc_mask;
    
    /* maximum number of equations solved */
    int max_nsys, nsys;
    int max_nsys_sq, nsys_sq;
    
    /* the current and specificly typed transport IDs */
    int itrns;
    int salinity_id;
    int temperature_id;
    int vorticity_id;  // the vorticity id in the array of con structs
    
    /* Rank and processor count for this model */
    /* Can be different from Grid values and from MPI_COMM_WORLD */
    int npes;
    int myid;
    int proc_flag;
#ifdef _MESSG
    MPI_Comm model_comm;
#endif
    
    /* sediment counts */
    int is_sediment_running; // flag that lets transport routines know that sediment is being applied
    int nsed, nclay, nsand, nlayers, nsfssi, nconti, nsilt;
    int ised; // current grain
    
    /* NSM WQ number */
    int NIO2, NIO3, NIH4, ORN, ORP, DIP, PHO4, Alg, CABOD, DIO, nnsm;
    
    /* mass conservation test variables */
    double initial_grid_mass;
    double grid_mass_error;
    
    /* nonlinear and linear total iteration counts */
    int nonlinear_it_total;
    int nonlinear_it_total_hvel;
    int nonlinear_it_total_wvel;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. //
    /* Wind Library struct: */
#ifdef WINDLIB
    SWINDLIB *windlib;
#endif
    // ABOVE LINES ADDED BY GAJANAN
    
#ifdef _ADH_HDF5
    HDF5 hdf5;
#endif
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for 2D-3D coupling.            //
    /* Mapping from nodes of submodel to the equation numbers: */
    int *fmap;       // alias
    int *fmap_wvel;  // alias
    // ABOVE LINES ADDED BY GAJANAN                                                                 //
    //////////////////////////////////////////////////////////////////////////////////////////////////
    
} SMODEL;

/*********************************************************/
/*********************************************************/
/* struct methods -------------------------------------- */

void smodel_init(SMODEL *, const char *);
void smodel_check(SMODEL *);
void smodel_print(SMODEL *);
void smodel_print_ts(SMODEL *, double);
void smodel_free(SMODEL *);
void smodel_open_input(SMODEL *);
void smodel_close_input(SMODEL * mod);
void smodel_open_output(SMODEL *);
void smodel_read_hot(SMODEL *);
void smodel_defaults(SMODEL *);


#endif
