/* standard header files */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <stdlib.h>
#include <stddef.h>
#include <limits.h>
#include <ctype.h>

#ifdef _MPI
#include <mpi.h>
#endif

#include "define.h"
#include "cards.h"
#include "macro.h"

#ifdef _PETSC
#include <petscksp.h>
//#include <petscts.h>
#endif

#include "debug.h"
#include "header_tl_alloc.h"
#include "type.h"

#include "fnctn_structs.h"
#include "friction_ext.h"


#include "fnctn.h"
#include "assert.h"
#include "constants.h"

/*******************************************************/
// CSTORM **********************************************/
//#ifdef _MESSG
//MPI_Comm cstorm_comm;
//#endif
// cjt :: since new AdH does not have global variables ...
//SMODEL *mod_cstorm;
/******************************************************/

/*******************************************************/
// WALL CLOCK TIMINGS *********************************/
double TIME_IN_HVEL_RESID;
double TIME_IN_HVEL_LOAD;
double TIME_IN_HVEL_BODY_RESID;
double TIME_IN_HVEL_BOUNDARY_RESID;

double TIME_IN_WVEL_RESID;
double TIME_IN_WVEL_LOAD;
double TIME_IN_WVEL_BODY_RESID;
double TIME_IN_WVEL_BOUNDARY_RESID;

double TIME_IN_2D_SW_RESID;
double TIME_IN_2D_SW_LOAD;
double TIME_IN_2D_SW_BODY_RESID;
double TIME_IN_2D_SW_BOUNDARY_RESID;

double TIME_IN_2D_TRANSPORT_BODY_RESID;
double TIME_IN_2D_TRANSPORT_BOUNDARY_RESID;
double TIME_IN_3D_TRANSPORT_BODY_RESID;
double TIME_IN_3D_TRANSPORT_BOUNDARY_RESID;

double TIME_IN_GW_BODY_RESID;
double TIME_IN_GW_BOUNDARY_RESID;
double TIME_IN_GW_LOAD;
double TIME_IN_GW_RESID;

double TIME_IN_NS3_BODY_RESID;
double TIME_IN_NS3_BOUNDARY_RESID;
double TIME_IN_NS3_LOAD;
double TIME_IN_NS3_RESID;

/*******************************************************/
// MISCELLANIOUS
//int ROTATE; // replace momentum on boundaries in 3D
//SDEBUG debug;
//double DEBUG_TIME;
//SSCREEN_OUTPUT screen_output;
//STESTCASE test_case_flag;
//
//double ***coupled_normal_flux;
//int **root_ids;
//double total_mass_flux_from_2d;
//double total_mass_flux_into_2d;
//double total_mass_flux_from_3d;
//double total_mass_flux_into_3d;
