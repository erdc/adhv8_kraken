#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <stdlib.h>
#include <stddef.h>
#include <limits.h>
#include <ctype.h>

#include "debug.h"
#include "header_tl_alloc.h"

#include "assert.h"

#include "stdbool.h"
#include "svect2d.h"        // dependencies :: smpi.h with _MESSG active
#include "svect.h"          // dependencies :: smpi.h with _MESSG active

#include "snode.h"
#include "selem_1d.h"       // dependencies :: svect2d.h
#include "selem_2d.h"       // dependencies :: svect.h, svect2d.h
#include "selem_3d.h"       // dependencies :: svect

#include "elem_physics.h"
#include "sgrid.h"          // dependencies :: selem_1d.h, selem_2d.h, selem_3d.h,snode.h, slist_items.h, smpi.h
#include "super_model.h"
#include "design_model.h"

#define UNSET_INT -3
#define UNSET_FLT -9999999.9
#define one_3  0.333333333333333333333333
#define ON 1            /* on */
#define OFF 0           /* off */
#define YES 1           /* yes */
#define NO -3           /* no */
#define NORMAL -1
#define NOT_QUITE_SMALL FLT_EPSILON
#define SMALL DBL_EPSILON
