#include "global_header.h"

/* finds the maximum norm across processors and returns */
/* the value and corresponding node number and coordinates */
void messg_max_norm_loc(double *mn_value,   /* max norm value */
                        int *mn_node,   /* max norm node */
                        SVECT * mn_coord /* max norm node coordinates */
#ifdef _MESSG
                        ,MPI_Comm ADH_COMM, int myid
#endif
                        
  )
{
#ifdef _MESSG
  int ierr_code;                /* the error code from an mpi call */
  struct
  {
    double value;               /* data value */
    int rank;                   /* data processor */
  } mn_in, mn_out;

  mn_in.rank = myid;
  mn_in.value = *mn_value;

  ierr_code = MPI_Allreduce(&mn_in, &mn_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, ADH_COMM);
  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }

  *mn_value = mn_out.value;

  ierr_code = MPI_Bcast(mn_node, 1, MPI_INT, mn_out.rank, ADH_COMM);
  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }
  ierr_code = MPI_Bcast(mn_coord, 3, MPI_DOUBLE, mn_out.rank, ADH_COMM);
  if (ierr_code != MPI_SUCCESS)
    {
      messg_err(ierr_code);
    }

#endif
  return;
}
