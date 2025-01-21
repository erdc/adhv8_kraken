#include "adh.h"
void ssw_alloc_init(SSW *sw) {

        //sw->series_wind_head = NULL;
        //sw->series_wave_head = NULL;
        //sw->series_wind_curr = NULL;
        //sw->series_wave_curr = NULL;
        //must allocate the memory
        
        printf("allocated ssw pointer\n");
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
        sw->elem_rhs_realloc = 0;

        // hydraulic structures
        //sw->nweir = 0;
        //sw->nflap = 0;
        //sw->nsluice = 0;
        //sw->weir = NULL;
        //sw->flap = NULL;
        //sw->sluice = NULL;

        //sw->nsw_nodes = nnodes;
        //sw->nd = (SSW_NODE *) tl_alloc(sizeof(SSW_NODE),nnodes);

        //sw->nsw_elem = nelems;
        //sw->elem = (SSW_NODE *) tl_alloc(sizeof(SSW_ELEM),nelems);

        //set default vals
        //nodal variable types
        // vertical grid speeds
        sw->GS = UNSET_INT; // 0
        sw->GS_OLD = UNSET_INT; //1
        // bed displacement from sediment transport
        sw->BED_DPL = UNSET_INT;// 2
        sw->BED_DPL_OLD = UNSET_INT;// 3
        sw->BED_DPL_OLDER = UNSET_INT;// 4
        // pressure and pressure pertubations
        sw->PRS = UNSET_INT;// 5
        sw->PRS_PLUS = UNSET_INT;// 6
        sw->PRS_MINUS = UNSET_INT;// 7
        //other dependent variables
        sw->DENSITY = UNSET_INT;// 8            // density
        sw->DPL_PERTURBATION = UNSET_INT;// 9   // displacement perturbation
        sw->CONT_ERROR = UNSET_INT;// 10              // nodal continuity errors
        sw->VERTICAL_NODE_FLUX = UNSET_INT;// 11
        sw->HYD_VISCOSITY = UNSET_INT;// 12
        sw->TRN = UNSET_INT;// 13
        //must deconstruct into double fields
        sw->GRAD_BED_X = UNSET_INT;// 14
        sw->GRAD_BED_Y = UNSET_INT;// 15          // bed gradient
        //tangent vector
        sw->TAN_VEC_X = UNSET_INT;// 16
        sw->TAN_VEC_Y = UNSET_INT;// 17
        //SVECT2D tanvec;
        sw->BED_ELEVATION = UNSET_INT;// 18 // node.z + bed displacement from sediment
        sw->SURFACE_VEL_X = UNSET_INT;// 19    //SVECT2D surface_vel;
        sw->SURFACE_VEL_Y = UNSET_INT;// 20
        sw->BOTTOM_VEL_X = UNSET_INT;// 21
        sw->BOTTOM_VEL_Y = UNSET_INT;// 22   //SVECT2D bottom_vel;
        // wind and wave stresses
        //not just doubles, need to think about this
        //SWIND winds;
        //SWAVE waves;
        sw->WIND_STRESS_X = UNSET_INT;// 23
        sw->WIND_STRESS_Y = UNSET_INT;// 24
        sw->WAVE_RADS_XX = UNSET_INT;// 25
        sw->WAVE_RADS_XY = UNSET_INT;// 26
        sw->WAVE_RADS_YY = UNSET_INT;// 27
        sw->WAVE_STRESS_X = UNSET_INT;// 28
        sw->WAVE_STRESS_Y = UNSET_INT;// 29

        //elemental variable types (doubles), none so far

        //elemental variable types (integers)
        sw->WD_FLAG = UNSET_INT;// 0


}
