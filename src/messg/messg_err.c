#include "global_header.h"

/*!
   \brief Message Passing Error, Report it, then Exit
 */
void messg_err(int ierr         /* the error code */
  )
{
#ifdef _MESSG
  char err_string[MPI_MAX_ERROR_STRING];    /* the string with the error description */
  int err_string_size = 0;      /* the size of the printed string */
  int ierr_code;  /* the error code from an mpi call */

  /* Get the error from MPI */
  ierr_code = MPI_Error_string(ierr, err_string, &err_string_size);
  if(ierr_code != MPI_SUCCESS) fprintf(stderr, "Error in MPI_Error_string  \n");
  fprintf(stderr, "Message passing error:  \n");
  fputs(err_string, stderr);
  fprintf(stderr, "\n");
#endif
 // messg_finalize();
  exit(1);
}
