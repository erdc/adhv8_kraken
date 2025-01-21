//#ifndef H_SADH_
//#define H_SADH_

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
#include <stdbool.h>
//Mark added
#include <umfpack.h>
#include <scotch.h>

#ifdef _ADH_HDF5
#include <hdf5.h>
#endif

#ifdef _MPI
#include <mpi.h>
#endif

#include "define.h"
#include "macro.h"

#ifdef _PETSC
//Mark changed
#include <petsc.h>
//#include <petscksp.h>
//#include <petscts.h>
#endif

#include "debug.h"
#include "header_tl_alloc.h"

#include "assert.h"
//#include "constants.h"

// STRUCTURES
#include "svect2d.h"
#include "svect.h"
#include "snode.h"
#include "stensor.h"
#include "selem_1d.h"
#include "selem_2d.h"
#include "selem_3d.h"
#include "squad.h"
#include "smodel.h"
#include "slist_items.h"
#include "smpi.h"
#include "smeteor.h"
#include "sflags.h"
#include "sstr_value.h"
#include "sgrid.h"
#include "sarray.h"


//Mark adding
#include "tokens.h"
#include "smat_grid.h"
#include "smat_sw.h"
#include "smat_gw.h"
#include "smat_transport.h"
#include "smat_physics.h"
#include "slin_sys.h"
#include "dofmaps.h"
#include "sdvar.h"
#include "ssw.h"


//
//Mark, not sure if this goes here?
#include "messg.h"

#include "sio.h"
#include "sseries.h"
#include "smodel_super.h"
#include "smodel_design.h"

#include "model_codes.h"
#include "models.h"

//Mark added
#include "residual.h"
#include "jacobian.h"
#include "la.h"
#include "newton.h"
//Mark added physics module
#include "fe.h"
#include "sw2.h"
#include "poisson.h"
#include "heat.h"
//Mark added time loop
#include "time_loop.h"


// FOLDERS

//Mark, does test need to be moved up?
#include "tools.h"
#include "fnctn_xdmf.h"
#include "fr_defs.h"

//Mark added
#include "testla.h"
#include "testresidual.h"
#include "testjacobian.h"
#include "testnewton.h"
#include "testnonlinearnewton.h"
#include "testtimeloop.h"
#include "testsw2_wd.h"
//#endif
