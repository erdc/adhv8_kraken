/*!
 \file math_simple.c
 \brief Simple Math Routines
 */
#include "global_header.h"

/*!
 \brief Initializes all entries of a (double) vector to zero
 
 \param n the length of the vectors (int)
 \param v pointer to the vector (double, length \a n)
 */
void solv_init_dbl(
                   int n,			/* the number of degrees of freedom */
                   double *v			/* the vector */
)
{
    int ii = 0;			/* loop counter */
    
    /* zeroes the array */
    for(ii = 0; ii < n; ii++)
    {
        v[ii] = DZERO;
    }
    return;
}

/*!
 \brief Initializes all entries of a (int) vector to zero
 
 \param n the length of the vectors (int)
 \param v pointer to the vector (int, length \a n)
 */
void solv_init_int(
                   int n,			/* the number of degrees of freedom */
                   int *v			/* the vector */
)
{
    int ii = 0;			/* loop counter */
    
    /* zeroes the array */
    for(ii = 0; ii < n; ii++)
    {
        v[ii] = 0;
    }
    return;
}

/*!
 \brief Function to allows us to sort list of ints
 
 \param *p1 Pointer to First Entry
 \param *p2 Pointer to Second Entry
 \return 1 If P1>P2, -1 If P1<P2, 0 If P1==P2
 */
int intcompare(
               const void *p1,
               const void *p2
               )
{
    int i = *((int *)p1);
    int j = *((int *)p2);
    
    if(i > j)
    {
        return (1);
    }
    if(i < j)
    {
        return (-1);
    }
    return (0);
}

/*!
 \brief Function to allows us to sort list of doubles
 
 \param *p1 Pointer to First Entry
 \param *p2 Pointer to Second Entry
 \return 1 If P1>P2, -1 If P1<P2, 0 If P1==P2
 */
int dblcompare(
               const void *p1,
               const void *p2
               )
{
    double i = *((double *)p1);
    double j = *((double *)p2);
    
    if(i > j)
    {
        return (1);
    }
    if(i < j)
    {
        return (-1);
    }
    return (0);
}

/*!
 \brief Check for very small values
 
 \param gamma pointer to double, possibly changed on output
 */

void check_for_small_dbl(
                         double *gamma			/* the value */
)
{
#if 1
    if(fabs(*gamma) < SMALL)
    {
        *gamma = SMALL;
    }
#endif
    return;
}

/* cjt - rewrite to return flag */
int is_double_small( double value ) {
    if ( fabs(value) < SMALL) return YES;
    return NO;
}
int compare_double( double v1, double v2) {
    if ( fabs(v1 - v2) < SMALL) return YES;
    return NO;
}

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


/*!
 \brief Check a Double Value for NaN, INF but do not exit!
 \param *X Double Array
 \param arraybounds Number of Array Values
 \param *filename File Name of Calling Routine
 \param linenumber Line Number of Calling Routine
 */
int Is_DoubleArray_Inf_or_NaN_noExit(double *X, int arraybounds, char *filename, int linenumber) {
    int flag = NO;
    int ii = 0;           /* Loop Counter */
    for(ii = 0; ii < arraybounds; ii++) {
        flag = Is_Double_Inf_or_NaN_noExit(X[ii], filename, linenumber);
    }
    return flag;
}

/*!
 \brief Check a Double Value for NaN, INF but do not exit!
 \param *X Double Array
 \param arraybounds Number of Array Values
 \param *filename File Name of Calling Routine
 \param linenumber Line Number of Calling Routine
 */
int Is_vectorArray_Inf_or_NaN_noExit(SVECT *X, int arraybounds, char *filename, int linenumber) {
    int flag = NO;
    int ii = 0;           /* Loop Counter */
    for(ii = 0; ii < arraybounds; ii++) {
        flag = Is_Double_Inf_or_NaN_noExit(X[ii].x, filename, linenumber);
        if (flag != 0) return flag;
        flag = Is_Double_Inf_or_NaN_noExit(X[ii].y, filename, linenumber);
        if (flag != 0) return flag;
        flag = Is_Double_Inf_or_NaN_noExit(X[ii].z, filename, linenumber);
        if (flag != 0) return flag;
    }
    return flag;
}

/*!
 \brief Check a Vector for NaN, INF but do not exit!
 \param *X Double Array
 \param arraybounds Number of Array Values
 \param *filename File Name of Calling Routine
 \param linenumber Line Number of Calling Routine
 */
int Is_Vector_Inf_or_NaN_noExit(SVECT X, char *filename, int linenumber) {
    int flag = NO;
    int ii = 0;           /* Loop Counter */
    flag = Is_Double_Inf_or_NaN_noExit(X.x, filename, linenumber);
    if (flag != 0) return flag;
    flag = Is_Double_Inf_or_NaN_noExit(X.y, filename, linenumber);
    if (flag != 0) return flag;
    flag = Is_Double_Inf_or_NaN_noExit(X.z, filename, linenumber);
    if (flag != 0) return flag;
    return flag;
}


/*!
 \brief Check a Double Value for NaN, INF but do not exit
 \param X Double Value
 \param *filename File Name of Calling Routine
 \param linenumber Line Number of Calling Routine
 */
int Is_Double_Inf_or_NaN_noExit(double X, char *filename, int linenumber) {
    int flag = NO;
    if(solv_isinf(X) != 0) {
        printf("%s:%d Double Value has value INF.\n", filename, linenumber);
        flag = 1;
    }
    if(solv_isnan(X) != 0) {
        printf("%s:%d Double Value has value NaN.\n", filename, linenumber);
        flag = 2;
    }
    return flag;
}

/*!
 \brief Returns the \f$l_2\f$ norm of a (double) vector
 
 \param n the length of the vectors (int)
 \param v1 pointer to the vector (double, length n)
 */
double solv_l2_norm_scaled(
                           int n,			/* the number of nodes */
                           double *v1,			/* the vector */
                           double denominator		/* Scale Factor */
)
{
    double value = DZERO;		/* the partial sum */
    double result = DZERO;	/* the result */
    
    /* Check Array First */
#ifdef _DEBUG_OJE
    Is_DoubleArray_Inf_or_NaN(v1, n, __FILE__, __LINE__);
#endif
    
    /* Perform Local Dot Product */
    value = solv_dot(n, v1, v1);
#ifdef _DEBUG_OJE
    Is_Double_Inf_or_NaN(value, __FILE__, __LINE__);
#endif
    
    /* Sum over Processors */
    /*#ifdef _MPI
     result = messg_dsum(value);
     #else*/
    result = value;
    /*#endif*/
    
#ifdef _DEBUG_OJE
    Is_Double_Inf_or_NaN(result, __FILE__, __LINE__);
#endif
    
    /* Divide By some Factor */
    result /= (denominator);
#ifdef _DEBUG_OJE
    Is_Double_Inf_or_NaN(result, __FILE__, __LINE__);
#endif
    
    /* Get Square Root, as this is a l_2 norm */
    value = sqrt(result);
#ifdef _DEBUG_OJE
    Is_Double_Inf_or_NaN(value, __FILE__, __LINE__);
#endif
    
    /* Return the l_2 norm */
    return (value);
}

/*!
 \brief Returns the \f$l_2\f$ norm of a (double) vector
 
 \param n the length of the vectors (int)
 \param v1 pointer to the vector (double, length n)
 */
double solv_l2_norm(
                    int n,			/* the number of degrees of freedom */
                    double *v1			/* the vector */
)
{
    double value = DZERO;		/* the partial sum */
    double result = DZERO;	/* the result */
    
    /* Check Array First */
#ifdef _DEBUG_OJE
    Is_DoubleArray_Inf_or_NaN(v1, n, __FILE__, __LINE__);
#endif
    
    /* Perform Local Dot Product */
    value = solv_dot(n, v1, v1);
#ifdef _DEBUG_OJE
    Is_Double_Inf_or_NaN(value, __FILE__, __LINE__);
#endif
    
    /* Sum over Processors */
    /*#ifdef _MPI
     result = messg_dsum(value);
     #else*/
    result = value;
    /*#endif*/
    
#ifdef _DEBUG_OJE
    Is_Double_Inf_or_NaN(result, __FILE__, __LINE__);
#endif
    
    /* Get Square Root, as this is a l_2 norm */
    value = sqrt(result);
#ifdef _DEBUG_OJE
    Is_Double_Inf_or_NaN(value, __FILE__, __LINE__);
#endif
    
    /* Return the l_2 norm */
    return (value);
}

/*!
 \brief Returns the \f$l_\infty\f (infinity) norm of a vector (maximum absolute value)
 
 \param n the length of the vectors (int)
 \param v1 pointer to the vector (double, length n)
 */
double solv_infty_norm(
                       int n,			/* the number of degrees of freedom */
                       double *v1			/* the vector */
#ifdef _MESSG
                 , MPI_Comm ADH_COMM
#endif
                       )
{
    double value = DZERO;		/* the partial sum */
    double result = DZERO;	/* the result */
    int ii = 0;			/* loop counter */
    
#ifdef _DEBUG_OJE
    /* Do a Quick Check in Debug Mode */
    Is_DoubleArray_Inf_or_NaN(v1, n, __FILE__, __LINE__);
#endif
    
    /* computes the value for this processor */
    for(ii = 0; ii < n; ii++)
    {
        value = MAX(value, fabs(v1[ii]));
    }
#ifdef _DEBUG_OJE
    Is_Double_Inf_or_NaN(value, __FILE__, __LINE__);
#endif
    
    /* finds maximum with other processors for the parallel run */
#ifdef _MESSG
    result = messg_dmax(value, ADH_COMM); 
#else
    result = value;
#endif
#ifdef _DEBUG_OJE
    Is_Double_Inf_or_NaN(result, __FILE__, __LINE__);
#endif
    
    /* returns the maximum */
    return (result);
}

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
