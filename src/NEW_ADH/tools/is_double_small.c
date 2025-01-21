#include "adh.h"
/* cjt - rewrite to return flag */
int is_double_small( double value ) {
    if ( fabs(value) < SMALL) return YES;
    return NO;
}
