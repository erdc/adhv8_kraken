/*!
 \file solv_blas_tools.c
 \brief Interface to blas routines, with native C option
 
 Accepts compile option _BLAS (must also link to library).
 See http://www.netlib.org/blas/ \n
 Otherwise, uses native C code (slower). \n
 Not all of these functions are true BLAS interfaces, but all
 do provide very basic computations.
 */
#include "global_header.h"
/*!
 \brief Copies one (double) vector into another, \f$ (y = x) \f$
 
 Entries of vector \f$ y \f$ are overwritten.
 
 \param n the length of the vectors (int)
 \param x pointer to the first vector (double, length n)
 \param y pointer to the second vector (double, length n)
 */
void solv_copy(
               int n,			/* the number of degrees of freedom */
               double *x,			/* the original */
               double *y			/* the target */
)
{
#ifdef _BLAS
    int IONE = 1;
    dcopy_(&n, x, &IONE, y, &IONE);
#else
    memcpy(y, x, ((size_t) (n)) * sizeof(double));
#endif
    return;
}

/*!
 \brief Performs daxpy operation for (double) vectors, \f$y = \alpha x  + y\f$
 
 \param n the length of the vectors (int)
 \param alpha scalar multiplier (double)
 \param x pointer to the first vector (double, length \a n)
 \param y pointer to the second vector (double, length \a n), changed on output
 */
void solv_daxpy(
                int n,			/* the number of degrees of freedom */
                double alpha,			/* the scalar */
                double *x,			/* x vector */
                double *y			/* y vector */
)
{
#ifdef _BLAS
    int IONE = 1;
    
    daxpy_(&n, &alpha, x, &IONE, y, &IONE);
#else
    int i = 0;			/* loop counter */
    
    /* adds the arrays */
    for(i = 0; i < n; i++)
        y[i] += alpha * x[i];
#endif
    return;
}

/*!
 \brief Returns the inner product, \f$(x \cdot y)\f$, of two (double) vectors
 
 Sums among all processors.
 
 \param n the length of the vectors (int)
 \param x pointer to the first vector (double, length \a n)
 \param y pointer to the second vector (double, length \a n)
 */
double solv_dot(
                int n,			/* the number of degrees of freedom */
                double *x,			/* the first vector */
                double *y			/* the second vector */
)
{
    double value = DZERO;		/* the partial sum */
    double result = DZERO;	/* the result */
#ifdef _BLAS
    int IONE = 1;
    
    value = ddot_(&n, x, &IONE, y, &IONE);
#else
    int i = 0;			/* loop counter */
    
    for(i = 0; i < n; i++)
        value += x[i] * y[i];
#endif
    
    /* sums with the other processors for the parallel run */
//#ifdef _MESSG
//    result = messg_dsum(value, cstorm_comm); // cjt :: only works for one grid in CSTORM!!!!  MPI_COMM_WORLD);
//#else
//    result = value;
//#endif
  
    /* sums with the other processors for the parallel run
     * done after function returns due to multi-model runs*/
    result = value;

    /* returns the sum */
    return (result);
}

/*!
 \brief Multiplies all entries of a (double) vector by a scalar
 
 \param n the length of the vectors (int)
 \param alpha scalar multiplier (double)
 \param x pointer to the vector (double, length n)
 */
void solv_scal(
               int n,			/* the number of degrees of freedom */
               double alpha,			/* the scalar */
               double *x			/* the vector */
)
{
#ifdef _BLAS
    int IONE = 1;
    
    dscal_(&n, &alpha, x, &IONE);
#else
    int i = 0;			/* loop counter */
    
    /* multiplies by a scalar */
    for(i = 0; i < n; i++)
        x[i] *= alpha;
#endif
    return;
}
