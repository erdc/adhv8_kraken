#include "adh.h"
void ssw_alloc_init(SSW *sw, int nnodes, int nelems) {

        //sw->series_wind_head = NULL;
        //sw->series_wave_head = NULL;
        //sw->series_wind_curr = NULL;
        //sw->series_wave_curr = NULL;


        // wetting and drying variables
        sw->drying_lower_limit = 0.;
        sw->drying_upper_limit = 0.;
        sw->wd_lower_tol = -0.1;
        sw->wd_upper_tol = 0.1;
        sw->wd_rate_lower_tol = -0.1;
        sw->wd_rate_upper_tol = 0.1;

        // surface water parameters
        sw->viscosity = 9.8E-7;
        sw->manning_units_constant = 1.0;
        sw->density = 1000.;
        sw->tau_pg = 0.5;

        // hydraulic structures
        //sw->nweir = 0;
        //sw->nflap = 0;
        //sw->nsluice = 0;
        //sw->weir = NULL;
        //sw->flap = NULL;
        //sw->sluice = NULL;

        sw->nsw_nodes = nnodes;
        sw->nd = (SSW_NODE *) tl_alloc(sizeof(SSW_NODE),nnodes);

        sw->nsw_elem = nelems;
        sw->elem = (SSW_NODE *) tl_alloc(sizeof(SSW_ELEM),nelems);


}