#ifndef H_MESSG_
#define H_MESSG_

//MPI related routines
double messg_dmax(double x      /* the value */
#ifdef _MESSG
                  , MPI_Comm ADH_COMM
#endif
  );

int messg_imin(int x            /* the value */
#ifdef _MESSG
               , MPI_Comm ADH_COMM
#endif
  );

int messg_imax(int x            /* the value */
#ifdef _MESSG
               , MPI_Comm ADH_COMM
#endif
  );

double messg_dsum(double x);


void comm_update_double(double *vec, int size_v, int npe, int rank);


#endif
