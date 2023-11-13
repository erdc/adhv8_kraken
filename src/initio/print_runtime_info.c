#include "global_header.h"

/*!
   \brief Print ADH runtime information to screen
 */
void print_runtime_info(SIO io)
{
  struct tm *tmp = NULL;
  time_t t = time(NULL);
#ifdef _PRE_ADH
  char exec_name[] = "pre-";
#else
  char exec_name[] = "";
#endif

  printf("\n -------------------\n Runtime Information\n -------------------\n");
  
  tmp = localtime(&t);  /* or gmtime, if you want GMT^H^H^HUTC */
  printf(" %sAdH execution Date/Time: %04d.%02d.%02d / %02d:%02d:%02d\n", exec_name,
    tmp->tm_year + 1900, tmp->tm_mon + 1, tmp->tm_mday, tmp->tm_hour, tmp->tm_min, tmp->tm_sec);
  printf(" Launching %sAdH with project name: %s\n", exec_name, io.proj_name);
  printf(" Launching %sAdH with run name: %s\n", exec_name,strlen(io.run_name) ? io.run_name : "(unspecified)");
#ifdef _MESSG
  int npes, ierr_code;
  ierr_code = MPI_Comm_size(cstorm_comm, &npes); // cjt :: only works for 1 grid in CSTORM :: MPI_COMM_WORLD, &npes);
  printf(" %sAdH was launched with %d processor%s\n", exec_name, npes, (npes > 1) ? "s" : "");
#endif

  fflush(NULL);
  return;
}
