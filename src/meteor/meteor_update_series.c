/* updates wind vector components when given through a bc time-series (cjt) */
#include "global_header.h"

void meteor_update_series(int nnodes, double t_prev, double density, SSERIES *series_head, double **values)
{

    int i, k, ixy, interval, indx;
    double winds_x, winds_y;
    SSERIES *series;

#ifdef _DEBUG
    if (debug.meteor) {
        sseries_printScreen(*series_head, 1);
       // if (myid == 0)
            printf("\n \n");
        printf("Checking memory before wind series update :: ");
        tl_check_all_pickets(__FILE__, __LINE__);
        //messg_barrier();
    }
#endif

    int type = series_head->type;

    for (i = 0; i < nnodes; i++) {

        for (indx = 0; indx < series_head->nvalues; indx++) {
            values[indx][i] = 0.;
        }

        series = series_head;
        while (series != NULL) {

            interval = sseries_get_interval(*series, t_prev);
            for (indx = 0; indx < series->nvalues; indx++) {
                series->ivalue[indx] = tc_eval_series(*series, interval, t_prev, indx);
            }
            if (density > 0.0) {
                for (indx = 0; indx < series->nvalues; indx++) {
                    values[indx][i] += series->ivalue[indx] * series->station->node_contrib[i];
                }
            }
#ifdef _DEBUG
            if (debug.meteor) {
                for (indx = 0; indx < series->nvalues; indx++) {
                    Is_DoubleArray_Inf_or_NaN(&values[indx][i], 1, __FILE__, __LINE__);
                    Is_DoubleArray_Inf_or_NaN(&values[indx][i], 1, __FILE__, __LINE__);
                }
            }
#endif
            series = series->next;
        }

        
    }

    series = NULL;

#ifdef _DEBUG
    if (debug.meteor) {
        printf("Checking memory after winds xy-series update :: ");
        tl_check_all_pickets(__FILE__, __LINE__);
        //messg_barrier();
    }
#endif

}
