#include "adh.h"
/*!
 \brief Check to see if a number is a NaN
 \param value Number to check
 */
int solv_isnan(
               double value			/* The Number to check */
)
{
#ifdef WINDOWS
    if(value != value)
    {
        return (1);
    }
    else
    {
        return (0);
    }
#else
    return (isnan(value));
#endif
}
