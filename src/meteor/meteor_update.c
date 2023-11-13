#include "global_header.h"

/******************************************************************/
/******************************************************************/
void values_alloc(int nnodes, int nvalues, double ***values) {
  int i = 0, j = 0;    
  (*values) = (double **) tl_alloc(sizeof(double *), nvalues);
  for (i = 0; i < nvalues; i++) {
    (*values)[i] = (double *) tl_alloc(sizeof(double), nnodes);
    for (j = 0; j < nnodes; j++) {
        (*values)[i][j] = 0.;
    }
  }
}

void values_free(int nnodes, int nvalues, double **values) {
    int i = 0;
    if (values != NULL) {
        for (i = 0; i < nvalues; i++) {
            values[i] = (double *) tl_free(sizeof(double), nnodes, values[i]);
        }
        values = (double **)  tl_free(sizeof(double *), nvalues, values);
    }
}


/******************************************************************/
/******************************************************************/

void meteor_update(SMODEL *mod, double time) {
    int i = 0;
    double **values = NULL;

    int nnodes = mod->grid->nnodes_sur;

    if (mod->flag.WIND_STATION) {  
        values_alloc(nnodes, mod->series_wind_head->nvalues, &values);
        meteor_update_series(nnodes, time, mod->density, mod->series_wind_head, values);
        if (mod->flag.SW2_FLOW) {
            for (i=0; i< nnodes; i++) {
                mod->sw->d2->winds[i].stress.x = values[0][i];
                mod->sw->d2->winds[i].stress.y = values[1][i];
            }
        } else if (mod->flag.SW3_FLOW) {
            for (i=0; i< nnodes; i++) {
                mod->sw->d3->winds[i].stress.x = values[0][i];
                mod->sw->d3->winds[i].stress.y = values[1][i];
            }
        }
        values_free(nnodes, mod->series_wind_head->nvalues, values);
    }

    if (mod->flag.WAVE_STATION) {
        values_alloc(nnodes, mod->series_wave_head->nvalues, &values);
        meteor_update_series(nnodes, time, mod->density, mod->series_wave_head, values);
        if (mod->flag.SW2_FLOW) {
            for (i=0; i<nnodes; i++) {
                mod->sw->d2->waves[i].rads.xx = values[0][i];
                mod->sw->d2->waves[i].rads.xy = values[1][i];
                mod->sw->d2->waves[i].rads.yy = values[2][i];
            }
        } else if (mod->flag.SW3_FLOW) {
            for (i=0; i<nnodes; i++) {
                mod->sw->d3->waves[i].rads.xx = values[0][i];
                mod->sw->d3->waves[i].rads.xy = values[1][i];
                mod->sw->d3->waves[i].rads.yy = values[2][i];
            }
        }

        values_free(nnodes, mod->series_wave_head->nvalues, values);
    }
//////////////////////////////////////////////////////////////////////////////////////////////////
// FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. //
#ifdef WINDLIB
    if (mod->flag.WIND_LIBRARY==ON) {
        values_alloc(nnodes, 2, &values);      // Necessarily "2" because of Wind_X and Wind_Y.
        swindlib_update(mod, mod->windlib, time, nnodes, values);
        if (mod->flag.SW2_FLOW) {
            for (i=0; i< nnodes; i++) {
                mod->sw->d2->winds[i].stress.x = values[0][i];
                mod->sw->d2->winds[i].stress.y = values[1][i];
            }
        } else if (mod->flag.SW3_FLOW) {
            for (i=0; i< nnodes; i++) {
                mod->sw->d3->winds[i].stress.x = values[0][i];
                mod->sw->d3->winds[i].stress.y = values[1][i];
            }
        }
        values_free(nnodes, 2, values);
    }
#endif
// ABOVE LINES ADDED BY GAJANAN                                                                 //
//////////////////////////////////////////////////////////////////////////////////////////////////
}
