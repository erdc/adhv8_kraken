#include "global_header.h"

// NOTES:
// -- NO = -3
// -- OFF = 0

void smodel_defaults(SMODEL *mod) {

    mod->sw = NULL;
    mod->ns = NULL;
    mod->con = NULL;
    mod->grid = NULL;
    mod->io = NULL;
    mod->mat = NULL;
    mod->bc_mask = NULL;
    mod->str_values = NULL;
    mod->amICoupled = NO;

#ifdef _SEDIMENT
    mod->sed = NULL;
#endif

#ifdef _ADH_GROUNDWATER
   mod->sgw = NULL;
#endif 
    /* flags **********************************************************/

    mod->flag.NS_FLOW = OFF;
    mod->flag.NS2_FLOW = OFF;
    mod->flag.NS3_FLOW = OFF;
    mod->flag.GW_FLOW = OFF;
    mod->flag.SW_FLOW = OFF;
    mod->flag.SW2_FLOW = OFF;
    mod->flag.SW3_FLOW = OFF;
    mod->flag.DIFFUSIVE_WAVE = OFF;
    
    mod->flag.MG = OFF;
    mod->flag.GRID_ADAPTION = OFF;
    mod->flag.ADAPTED_THE_GRID = NO;
    mod->flag.GRID_REFINED = NO;
    mod->flag.GRID_UNREFINED = NO;
    mod->flag.TRANSPORT = OFF;
    mod->flag.NS2_TRANSPORT = OFF;
    mod->flag.NS3_TRANSPORT = OFF;
    mod->flag.SW2_TRANSPORT = OFF;
    mod->flag.SW3_TRANSPORT = OFF;
    mod->flag.VORTICITY = OFF;
    mod->flag.SEDIMENT = OFF;
    mod->flag.SEDLIB = OFF;
    mod->flag.PRN_ADPT = OFF;
    mod->flag.OUTPUT = OFF;
    mod->flag.ICM = OFF;
    mod->flag.NSM = OFF;
    mod->flag.WAVE = OFF;
    mod->flag.WIND = OFF;
    mod->flag.WAVE_STATION = OFF;
    mod->flag.WIND_STATION = OFF;
    mod->flag.EOS = 0; /* Always default to linearlized EOS */
 
//////////////////////////////////////////////////////////////////////////////////////////////////
// FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. //
#ifdef WINDLIB
    mod->flag.WIND_LIBRARY = OFF;
#endif
// ABOVE LINES ADDED BY GAJANAN                                                                 //
//////////////////////////////////////////////////////////////////////////////////////////////////

    mod->flag.CONVEYANCE = OFF;
    mod->flag.UNITS = OFF;
    mod->flag.ICE = OFF;
    mod->flag.INS = OFF;
    mod->flag.CORIOLIS = OFF;
    mod->flag.TIME_ADAPT = OFF;
    mod->flag.SOLVE_ATF = ON;
    mod->flag.STEADY_STATE = OFF;
    mod->flag.TIME_ADAPT_FAIL = NO;
    mod->flag.UMFPACK_FAIL = NO;
    mod->flag.TIDE = OFF;
    mod->flag.UNREFINE = NO;
    mod->flag.BAROCLINIC = OFF;
    mod->flag.CSTORM_WSID = OFF;

//////////////////////////////////////////////////////////////////////////////////////////////////
// FOLLOWING LINES ADDED BY GAJANAN [ gkc May 2017 ]. These are for XDMF parallel I/O.          //
#ifdef _ADH_HDF5
    mod->flag.PC_FILE_XDMF = OFF;
#endif
// ABOVE LINES ADDED BY GAJANAN                                                                 //
//////////////////////////////////////////////////////////////////////////////////////////////////

    mod->flag.FLUX_WEIGHTED_NORMALS = OFF;

#ifdef _ADH_GROUNDWATER
    mod->flag.GW_SALINITY = OFF;
    mod->flag.GW_TRANSPORT= OFF;
    mod->flag.GW_REACTION = OFF;
    mod->flag.RAY_TRACING = OFF;
    mod->flag.SOCKETS     = OFF;
    mod->flag.METEOROLOGY = OFF;    
#endif    
    /* variables *****************************************************/


    /* wetting and drying variables */
    mod->drying_lower_limit = 0.;
    mod->drying_upper_limit = 0.;
    mod->wd_lower_tol = -0.1;
    mod->wd_upper_tol = 0.1;
    mod->wd_rate_lower_tol = -0.1;
    mod->wd_rate_upper_tol = 0.1;

    mod->ntransport = 0;
    mod->itrns = UNSET_INT;
    mod->salinity_id = UNSET_INT;
    mod->vorticity_id = UNSET_INT;

    /* tidal tidal series */
    mod->ntides = 0;

    /* node string */
    mod->nstring = 0;
    mod->nstring_rt = 0;

    /* series */
    mod->nseries = 0;
    mod->series_head = NULL;
    mod->series_curr = NULL;
    mod->series_dt = NULL;
    mod->series_out = NULL;
    mod->series_wind_head = NULL;
    mod->series_wave_head = NULL;
    mod->series_wind_curr = NULL;
    mod->series_wave_curr = NULL;

#ifdef _ADH_GROUNDWATER
    mod->series_gw_psk_head = NULL; mod->series_gw_psk_curr = NULL;
#endif
    /* solver tolerances */
    mod->max_nonlin_it = 10;
    mod->nalloc_inc = 1;
    mod->nblock = 1;
    mod->inc_nonlin = UNSET_FLT;
    mod->tol_nonlin = UNSET_FLT;

    /* petrov coefficient */
    mod->tau_pg = 0.5;

    /* flag level of output */
    mod->out_level = 1;
    mod->o_flag = 0;   

    /* physics parameters */
    mod->viscosity = 9.8E-7;
    mod->manning_units_constant = 1.0;
    mod->gravity = 9.8;
    mod->density = 1000.;

    /* time-steping */
    mod->t_init = 0.;
    mod->t_prev = 0.;
    mod->t_final = 0.;
    mod->tau_temporal = 0.;
    mod->dt = 0.;
    mod->old_dt = 0.;
    mod->dt_err = 0.;
    mod->ientry_out = 0;       // current entry of the output series
    mod->green_ampt = FALSE;
    /* max system of equations */
    mod->nsys = 0;
    mod->nsys_sq = 0;
    mod->max_nsys = 0;
    mod->max_nsys_sq = 0;
    
    /* load perturbation */
    mod->perturbation  = sqrt(SMALL);
    mod->proc_flag = 0;

    /* MPI */
#ifdef _MESSG
    int ierr, myid, npes;
    //ierr = MPI_Comm_dup(cstorm_comm, &mod->model_comm);
    //if(ierr !=  MPI_SUCCESS) messg_err(ierr);
    myid = messg_comm_rank(mod->model_comm);
    mod->myid =  myid;
    npes = messg_comm_size(mod->model_comm);
    mod->npes = npes;
#else
    mod->myid = 0;
    mod->npes = 1;
#endif

    /* optional packages *********************************************/

    mod->nweir = 0;
    mod->nflap = 0;
	mod->nsluice = 0;
    mod->weir = NULL;
    mod->flap = NULL;
	mod->sluice = NULL;

    /* sediment counts */
    mod->is_sediment_running=FALSE;
    mod->nsed = 0;
    mod->nclay = 0;
    mod->nsand = 0;
    mod->nlayers = 0;
    mod->nsfssi = 0;
    mod->nconti = 0;
    mod->nsilt = 0;
    mod->ised = 0;

    /* global initial water mass */
    mod->initial_grid_mass = 0.;
    mod->grid_mass_error = 0.;

    /* nonlinear and linear total iteration counts */
    mod->nonlinear_it_total = 0;
    mod->nonlinear_it_total_hvel = 0;
    mod->nonlinear_it_total_wvel = 0;


}
