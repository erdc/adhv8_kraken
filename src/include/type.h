
// stand-alone structures --------------------------------------
//#include "stdbool.h"
//#include "svect2d.h"        // dependencies :: smpi.h with _MESSG active
//#include "svect.h"          // dependencies :: smpi.h with _MESSG active
//#include "stensor.h"
//#include "snode.h"
//#include "sdof.h"
//#include "sdebug.h"
//#include "sfile_output.h"
//#include "sscreen_output.h"
//#include "stestcase.h"
//#include "sflags.h"
//#include "stime.h"
//#include "sio.h"

// dependent structures ----------------------------------------
//#include "squad.h"
//#include "selem_1d.h"       // dependencies :: svect2d.h
//#include "selem_2d.h"       // dependencies :: svect.h, svect2d.h
//#include "selem_3d.h"       // dependencies :: svect
//#include "slist_items.h"    // dependencies :: svect.h, svect2d.h
//#include "smpi.h"           // dependencies :: slist_items.h
//#include "sgrid.h"          // dependencies :: selem_1d.h, selem_2d.h, selem_3d.h, snode.h, slist_items.h, smpi.h
//#include "ssolve.h"         // dependencies :: smpi.h with _MESSG active
//#include "sarray.h"         // dependencies :: smpi.h with _MESSG active
//#include "svect2d.h"        // dependencies :: smpi.h with _MESSG active
//#include "svect.h"          // dependencies :: smpi.h with _MESSG active
//#include "smeteor.h"
//#include "sns_2d.h"         // dependencies :: svect2d.h, smeteor.h
//#include "sns_3d.h"         // dependencies :: svect2d.h, svect.h, selem1d.h, selem2d.h
//#include "sns.h"            // dependencies :: sns_2d.h, sns_3d.h
//#include "ssw_2d.h"         // dependencies :: svect2d.h, smeteor.h
//#include "ssw_3d.h"         // dependencies :: svect2d.h, svect.h, selem1d.h, selem2d.h
//#include "ssw.h"            // dependencies :: ssw_2d.h, ssw_3d.h

//#include "scon.h"           // dependencies :: svect2d.h
//#include "smat.h"           // dependencies :: stensor.h
//#include "smpi.h"           // dependencies :: slist_items.h

#ifdef _SEDIMENT
#include "ssediment.h"      // dependencies :: svect2d.h, ssedlib.h
#endif
#ifdef _ADH_GROUNDWATER
#include "sgw.h"
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////
// FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. //
#ifdef WINDLIB
#include "swindlib.h"
#endif
// ABOVE LINES ADDED BY GAJANAN                                                                 //
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
// FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for XDMF parallel I/O.         //
#ifdef _ADH_HDF5
#include "adh_hdf5.h"
#endif
// ABOVE LINES ADDED BY GAJANAN                                                                 //
//////////////////////////////////////////////////////////////////////////////////////////////////

//#include "sseries.h"        // dependencies :: smeteor.h
//#include "sstr_value.h"     // dependencies :: sgrid <--- in function, take out ...
//#include "sstructures_value.h" // depedencies:: WEIR_STRUCT_C and FLAP_STRUCT_C
//#include "smodel.h"         // dependencies :: sio.h

//////////////////////////////////////////////////////////////////////////////////////////////////
// FOLLOWING LINES ADDED BY GAJANAN [ gkc April 2015 ]. These are for 2D-3D coupled models.     //
//#include "sinterface.h" // dependencies :: interface_nodelist.h
//#include "ssuperinterface.h"
// ABOVE LINES ADDED BY GAJANAN                                                                 //
//////////////////////////////////////////////////////////////////////////////////////////////////


