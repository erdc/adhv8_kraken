/* This routine scales the time step for iteration failure */
#include "adh.h"
void tc_scale(double * dt
#ifdef _MESSG
              , MPI_Comm ADH_COMM
#endif
              ) {

  /* scales the time step */
  (*dt) *= DT_REDUCE_FACTOR;
#ifdef _MESSG
  /* makes sure everyone is on the same dt */
  *dt = messg_dmin(*dt, ADH_COMM);
#endif

  /* stop if the time step is too small */
  if(*dt <= SMALL)
    tl_error("Time step is ridiculous.\n");
}
