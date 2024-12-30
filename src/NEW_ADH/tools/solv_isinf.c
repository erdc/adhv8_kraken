#include "adh.h"
/*!
 \brief Check to see if a number is an Infinity
 \param value Number to check
 */
int solv_isinf(
               double value			/* The Number to check */
)
{
#ifdef WINDOWS
    if((value * 0.0) != 0.0)
    {
        return (1);
    }
    else
    {
        return (0);
    }
#else
    return (isinf(value));
#endif
}
