
#include "stdbool.h"
#include "svect2d.h"        // dependencies :: smpi.h with _MESSG active
#include "svect.h"          // dependencies :: smpi.h with _MESSG active
#include "stensor.h"
#include "smeteor.h"
#include "snode.h"
#include "squad.h"
#include "selem_1d.h"       // dependencies :: svect2d.h
#include "selem_2d.h"       // dependencies :: svect.h, svect2d.h
#include "selem_3d.h"       // dependencies :: svect
#include "slist_items.h"    // dependencies :: svect.h, svect2d.h
#include "smpi.h"           // dependencies :: slist_items.h
#include "sgrid.h"          // dependencies :: selem_1d.h, selem_2d.h, selem_3d.h, snode.h, slist_items.h, smpi.h
//#include "ssolve.h"         // dependencies :: smpi.h with _MESSG active
#include "sarray.h"         // dependencies :: smpi.h with _MESSG active
#include "svect2d.h"        // dependencies :: smpi.h with _MESSG active
#include "svect.h"          // dependencies :: smpi.h with _MESSG active
#include "smpi.h"           // dependencies :: slist_items.h

#include "elem_physics.h"
#include "sgrid.h"          // dependencies :: selem_1d.h, selem_2d.h, selem_3d.h,snode.h, slist_items.h, smpi.h
#include "super_model.h"
#include "design_model.h"
