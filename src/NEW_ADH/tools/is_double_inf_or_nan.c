#include "adh.h"
/*!
 \brief Check a Double Value for NaN, INF
 
 \param *X Double Array
 \param arraybounds Number of Array Values
 \param *filename File Name of Calling Routine
 \param linenumber Line Number of Calling Routine
 */
void Is_DoubleArray_Inf_or_NaN(
                               double *X,
                               int arraybounds,
                               char *filename,
                               int linenumber
                               )
{
    int ii = 0;			/* Loop Counter */
    for(ii = 0; ii < arraybounds; ii++)
    {
        Is_Double_Inf_or_NaN(X[ii], filename, linenumber);
    }
    return;
}

/*!
 \brief Check a Double Value for NaN, INF
 
 \param X Double Value
 \param *filename File Name of Calling Routine
 \param linenumber Line Number of Calling Routine
 */
void Is_Double_Inf_or_NaN(
                          double X,
                          char *filename,
                          int linenumber
                          )
{
    if(solv_isinf(X) != 0)
    {
        printf("%s:%d Double Value has value INF.\n", filename, linenumber);
        exit(1);
    }
    if(solv_isnan(X) != 0)
    {
        printf("%s:%d Double Value has value NaN.\n", filename, linenumber);
        exit(1);
    }
    return;
}
