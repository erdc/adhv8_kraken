#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>

#include "header_tl_alloc.h"
#include "debug.h"

#include "svect2d.h"

#include "ssedlib.h"
#include "sedlib_prototypes.h"


/* other definitions */
#define NDIM 3          /* the maximum number of dimensions */
#define ON 1            /* on */
#define OFF 0           /* off */
#define YES 1           /* yes */
#define NO -3           /* no */
#define TRUE 1          /* true */
#define FALSE 0         /* false */

/* initialization definitions */
#define UNSET_INT -3        /* an integer variable that has not been set */
#define UNSET_FLT -3.0      /* a float variable that has not been set */
