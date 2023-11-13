#include "global_header.h"

/* Write text to standard output on root process only */
void root_print(
  char *text  /* the text message to be printed */
  )
{
#ifdef _MESSG
  int myid, ierr_code;
  ierr_code = MPI_Comm_rank(cstorm_comm, &myid); // cjt :: only works for 1 grid in CSTORM!!! MPI_COMM_WORLD, &myid);
#else
  int myid = 0;
#endif

  if (myid <= 0)
    {
      fprintf(stdout, "%s\n", text);
    }
  return;
}
