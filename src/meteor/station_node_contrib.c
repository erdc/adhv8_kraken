/* JNT 7-16-2002  -- cjt updated */
/* this routine calcultes the weighting factors for the wind shear stresses */

#include "global_header.h"

void station_node_contrib(SSERIES *head, double *x, double *y, int nstations)
{
  int i, n, k = 0;              /* variables used for looping */
  double **R;                   /* the distance from the wind station to the node */
  double *sumR;                 /* the sum of the distances from the stations to a node */
  double *sum_contrib;
  double node_contrib_sum;

#ifdef _DEBUG
  if (debug.winds) {
//    if (myid==0) printf("\n \n");
    printf("Checking memory before tl_eval_windcontrib ::"); tl_check_all_pickets(__FILE__, __LINE__);
    //messg_barrier();
  }
#endif

  SSERIES *series = head;

  int nnodes = head->nnodes;
  assert(nnodes>0);

  /* if only 1 station, set weights to 1 and bail out */
  if (nstations == 1) {
    for (i = 0; i < nnodes; i++) {
      series->station->node_contrib[i] = 1.;
    }
    return;
  }

  /* multiple stations continue forward */
  sumR = (double *) tl_alloc(sizeof(double), nnodes);
  sum_contrib = (double *) tl_alloc(sizeof(double), nnodes);
  for (i = 0; i < nnodes; i++) {
	sumR[i] = 0.;
	sum_contrib[i] = 0.;
  }

  R = (double **) tl_alloc(sizeof(double *), nstations);
  for (n = 0; n < nstations; n++) {
    R[n] = (double *) tl_alloc(sizeof(double), nnodes);
    for (i = 0; i < nnodes; i++) {
      R[n][i] = 0.;
    }
  }

  for (i = 0; i < nnodes; i++) {

    series = head; sumR[i] = 0.; k = 0;
    while (series != NULL) {
        R[k][i] = pow(pow(series->station->x - x[i], 2.0) + pow(series->station->y - y[i], 2.0), 1. / 2.);
        sumR[i] += R[k][i];
        k++;
        series = series->next;
    }
    
    // sanity check
    if (k != nstations) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        printf("k: %d \t nstations: %d\n",k,nstations);
        tl_error(">> number of meterologic stations is no correct in interpolation!");
    }
    k = 0; node_contrib_sum = 0.0;

    series = head; k = 0;
    while (series != NULL) {
        if (sumR[i] < SMALL) {
            series->station->node_contrib[i] = 1.0;
        } else {
            series->station->node_contrib[i] = (sumR[i] - R[k][i]) / sumR[i];
        }
        node_contrib_sum += series->station->node_contrib[i];
        k++;
        series = series->next;
    } 

    series = head;
    while (series != NULL) {
        if (node_contrib_sum < SMALL) {
            series->station->node_contrib[i] = 1.0;
        } else {
            series->station->node_contrib[i] /= node_contrib_sum;
        }
        series = series->next;
    }

  }

#ifdef _DEBUG
  if (debug.winds) {
    for (i = 0; i < nnodes; i++) {
        series = head;
        while (series != NULL) {
            Is_DoubleArray_Inf_or_NaN(&(series->station->node_contrib[i]), 1, __FILE__, __LINE__);
            series = series->next;
        }
    }
  }
#endif

  /* deallocate */
  series = NULL;
  sumR = (double *) tl_free(sizeof(double), nnodes, sumR);
  sum_contrib = (double *) tl_free(sizeof(double), nnodes, sum_contrib);
  for (n = 0; n < nstations; n++) {
    R[n] = (double *) tl_free(sizeof(double), nnodes, R[n]);
  }
  R = (double **) tl_free(sizeof(double *), nstations, R);


#ifdef _DEBUG
  if (debug.winds) {
    printf("Checking memory after tl_eval_windcontrib ::"); tl_check_all_pickets(__FILE__, __LINE__);
    //messg_barrier();
  }
#endif

}
