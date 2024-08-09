
#include "adh.h"

/***********************************************************/
/***********************************************************/
/***********************************************************/

void smeteor_file_init(SMETEOR_FILE *meteor) {
    meteor->n = 0;
    meteor->dt = 0;
    meteor->tprev = 0;
    meteor->tnext = 0;
    meteor->fin = NULL;
}

void smeteor_station_init(SMETEOR_STATION *station, int nnodes) {
    station->x = UNSET_FLT;
    station->y = UNSET_FLT;
    int i = 0;
    for (i = 0; i <nnodes; i++) {
        station->node_contrib[i] = 0.;
    }
}

void swind_init(SWIND *wind) {
    svect2d_init(&(wind->stress));
}

void swave_init(SWAVE *wave) {
    stensor2d_init(&(wave->rads));
    svect2d_init(&(wave->stress));
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void smeteor_station_alloc(SMETEOR_STATION **station, int nnodes) {
    (*station) = (SMETEOR_STATION *) tl_alloc(sizeof(SMETEOR_STATION), 1);
    (*station)->node_contrib = (double *) tl_alloc(sizeof(double), nnodes);
    smeteor_station_init(*station, nnodes);
}

void swave_alloc(SWAVE **wave, int nnodes) {
    int i = 0;
    SWAVE *lwave; // alias
    (*wave) = (SWAVE *) tl_alloc(sizeof(SWAVE), nnodes);
    lwave = *wave;
    for (i = 0; i < nnodes; i++) {
        swave_init(&lwave[i]);
    }
}

void swind_alloc(SWIND **wind, int nnodes) {
    int i = 0;
    SWIND *lwind; // alias
    (*wind) = (SWIND *) tl_alloc(sizeof(SWIND), nnodes);
    lwind = *wind;
    for (i = 0; i < nnodes; i++) {
        swind_init(&lwind[i]);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

void smeteor_station_realloc(SMETEOR_STATION **station, int nnodes_new, int nnodes_old) {
    (*station)->node_contrib = (double *) tl_realloc(sizeof(double), nnodes_new, nnodes_old, (*station)->node_contrib);
    int i;
    for (i = nnodes_old; i < nnodes_new; i++) {
        (*station)->node_contrib[i] = 0.;
    }
}

void swave_realloc(SWAVE **wave, int nnodes_old, int nnodes_new) {
    
    int i = 0;
    SWAVE *lwave; // alias
    (*wave) = (SWAVE *) tl_realloc(sizeof(SWAVE), nnodes_new, nnodes_old, (*wave));
    lwave = *wave;
    for (i = nnodes_old; i < nnodes_new; i++) {
        swave_init(&lwave[i]);
    }
}

void swind_realloc(SWIND **wind, int nnodes_old, int nnodes_new) {
    int i = 0;
    SWIND *lwind; // alias
    (*wind) = (SWIND *) tl_realloc(sizeof(SWIND), nnodes_new, nnodes_old, (*wind));
    lwind = *wind;
    for (i = nnodes_old; i < nnodes_new; i++) {
        swind_init(&lwind[i]);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/

//void smeteor_station_renumber(SMETEOR_STATION *station, int max_nnode, int *new_numbers, int *order_tmp, double *dtmp){
//    node_renumber_double(max_nnode, station->node_contrib, dtmp, new_numbers, order_tmp);
//}
//
//void swave_renumber(SWAVE *wave, int max_nnode, int *new_number, int *order_tmp) {
//    SWAVE *tmp_wave;
//    int i;
//    tmp_wave = (SWAVE *)tl_alloc(sizeof(SWAVE), max_nnode);
//    
//    for (i = 0; i < max_nnode; i++) {
//        
//        if(new_number[i] != i) {
//            tmp_wave[new_number[i]].rads.xx = wave[new_number[i]].rads.xx;
//            tmp_wave[new_number[i]].rads.xy = wave[new_number[i]].rads.xy;
//            tmp_wave[new_number[i]].rads.yy = wave[new_number[i]].rads.yy;
//            tmp_wave[new_number[i]].stress.x = wave[new_number[i]].stress.x;
//            tmp_wave[new_number[i]].stress.y = wave[new_number[i]].stress.y;
//            if (order_tmp[i]==0) {
//                wave[new_number[i]].rads.xx = wave[i].rads.xx;
//                wave[new_number[i]].rads.xy = wave[i].rads.xy;
//                wave[new_number[i]].rads.yy = wave[i].rads.yy;
//                wave[new_number[i]].stress.x = wave[i].stress.x;
//                wave[new_number[i]].stress.y = wave[i].stress.y;
//            }
//            else {
//                wave[new_number[i]].rads.xx = tmp_wave[i].rads.xx;
//                wave[new_number[i]].rads.xy = tmp_wave[i].rads.xy;
//                wave[new_number[i]].rads.yy = tmp_wave[i].rads.yy;
//                wave[new_number[i]].stress.x = tmp_wave[i].stress.x;
//                wave[new_number[i]].stress.y = tmp_wave[i].stress.y;
//            }
//        }
//    }
//    
//    tmp_wave = (SWAVE *)tl_free(sizeof(SWAVE), max_nnode, tmp_wave);
//    
//    return;
//}
//
//void swind_renumber(SWIND *wind, int max_nnode, int *new_numbers, int *order_tmp, SVECT2D *vtmp){
//    
//    node_renumber_vect2d(max_nnode, &(wind->stress), vtmp, new_numbers, order_tmp);
//    
//}
//
//

/***********************************************************/
/***********************************************************/
/***********************************************************/

void smeteor_station_copy(SMETEOR_STATION *dest, SMETEOR_STATION *src, int nnodes) {
    
    
    if(src == NULL) {
        dest = NULL;
        printf("SOURCE IS NULL (file:line) %s:%d\n", __FILE__, __LINE__);
        return;
    }
    dest->x = src->x;
    dest->y = src->y;
    
    if(src->node_contrib == NULL) {
        dest->node_contrib = NULL;
        printf("station node contribution IS NULL in copy(file:line) %s:%d\n", __FILE__, __LINE__);
        return;
    }
    
    int i=0;
    for (i=0; i<nnodes; i++) {
        dest->node_contrib[i] = src->node_contrib[i];
    }
    
    
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
void smeteor_station_free(SMETEOR_STATION *station, int nnodes) {
    
    // free node contribution
    if (station->node_contrib != NULL) {
        station->node_contrib = (double *) tl_free(sizeof(double), nnodes, station->node_contrib);
        station->node_contrib = NULL;
    }
    
    // free the meteorologic station
    station = (SMETEOR_STATION *) tl_free(sizeof(SMETEOR_STATION), 1, station);
    
}

void swave_free(SWAVE *wave, int nnodes) {
    if (wave != NULL) {
        wave = (SWAVE *) tl_free(sizeof(SWAVE), nnodes, wave);
    }
    wave = NULL;
}

void swind_free(SWIND *wind, int nnodes) {
    if (wind != NULL) {
        wind = (SWIND *) tl_free(sizeof(SWIND), nnodes, wind);
    }
    wind = NULL;
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
void swave_elem2d_local_init(SWAVE *waves) {
    int i;
    for (i=0; i<NDONTRI; i++) {
        swave_init(&waves[i]);
    }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
void swind_elem2d_local_init(SWIND *winds) {
    int i;
    for (i=0; i<NDONTRI; i++) {
        swind_init(&winds[i]);
    }
}

///***********************************************************/
///***********************************************************/
///***********************************************************/
//void swind_elem2d_get_local(SWIND *winds, SWIND *elem_winds, int *nodeIDs) {
//    int i;
//    for (i=0; i<NDONTRI; i++) {
//        elem_winds[i].stress.x = winds[nodeIDs[i]].stress.x;
//        elem_winds[i].stress.y = winds[nodeIDs[i]].stress.y;
//    }
//}
//
///***********************************************************/
///***********************************************************/
///***********************************************************/
//SVECT2D swind_elem2d_get_local_stress(int ndim, SWIND *winds, SELEM_2D *elem2d, int *grid2d_nodeIDs, double elem_depth_avg, double gravity, double density, double windatt, int wind_flag, int WIND_LIBRARY, int NWS) {
//    int i, nnodes2d = elem2d->nnodes;
//    
//    SWIND elem_winds[nnodes2d]; swind_elem2d_local_init(elem_winds);
//    SVECT2D wind_stress;  svect2d_init(&wind_stress);
//    
//    assert(ndim == 2 || ndim == 3);
//    swind_elem2d_get_local(winds, elem_winds, grid2d_nodeIDs);
//    
//    double shear_x, shear_y, xws, yws, wind_depth;
//    wind_depth = elem_depth_avg;
//    for (i = 0; i < nnodes2d; i++) {
//        shear_x = 0.0; shear_y = 0.0; xws = elem_winds[i].stress.x; yws = elem_winds[i].stress.y;
//        if (density > 0.0) {
//            //printf("wind_flag: %d\n",mod->mat[imat].sw->wind_flag);
//            if      (wind_flag == 0) { shear_x = xws; shear_y = yws; }
//            else if (wind_flag == 1) winds_proc_wu(gravity, xws, yws, &shear_x, &shear_y);
//            else if (wind_flag == 2) winds_proc_teeter(gravity, wind_depth, xws, yws, &shear_x, &shear_y);
//        }
//        //////////////////////////////////////////////////////////////////////////////////////////////////
//        // FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. //
//#ifdef WINDLIB
//        if (density > 0.0 && WIND_LIBRARY == ON) {
//            switch (NWS){
//                case 1:
//                    shear_x = xws; shear_y = yws;
//                    break;
//                case 2:
//                    shear_x = xws; shear_y = yws;
//                    break;
//                case 3:
//                    winds_proc_teeter(gravity, wind_depth, xws, yws, &shear_x, &shear_y);
//                    break;
//                case 4:
//                    winds_proc_teeter(gravity, wind_depth, xws, yws, &shear_x, &shear_y);
//                    break;
//                case 5:
//                    winds_proc_teeter(gravity, wind_depth, xws, yws, &shear_x, &shear_y);
//                    break;
//                case 6:
//                    winds_proc_teeter(gravity, wind_depth, xws, yws, &shear_x, &shear_y);
//                    break;
//                case 7:
//                    shear_x = xws; shear_y = yws;
//                    break;
//                case 8:
//                    winds_proc_teeter(gravity, wind_depth, xws, yws, &shear_x, &shear_y);
//                    break;
//                default:
//                    break;
//            }
//        }
//#endif
//        // ABOVE LINES ADDED BY GAJANAN                                                                 //
//        //////////////////////////////////////////////////////////////////////////////////////////////////
//        
//        //printf("wind shear: %30.20f %30.20f\n",shear_x,shear_y);
//        wind_stress.x += shear_x;
//        wind_stress.y += shear_y;
//    }
//    /* average over triangle */
//    wind_stress.x *= pow(windatt,2.0) / (nnodes2d * density);
//    wind_stress.y *= pow(windatt,2.0) / (nnodes2d * density);
//    //printf("wind stress: %30.20f %30.20f\n",wind_stress.x,wind_stress.y);
//    
//    return wind_stress;
//}
//
///***********************************************************/
///***********************************************************/
///***********************************************************/
//SVECT2D swave_elem2d_get_local_stress(int ndim, int CSTORM_WSID, SWAVE *waves, SELEM_2D elem2d, int *grid2d_nodeIDs) {
//    int i;
//    
//    SWAVE elem_waves[NDONTRI]; swave_elem2d_local_init(elem_waves);
//    SVECT2D wave_stress;  svect2d_init(&wave_stress);
//    
//    assert(ndim == 2 || ndim == 3);
//    swave_elem2d_get_local(waves, elem_waves, grid2d_nodeIDs);
//    
//    if (CSTORM_WSID == OFF || CSTORM_WSID == 7) {
//        for (i = 0; i<NDONTRI; i++) {
//            wave_stress.x -= (elem2d.grad_shp[i].x * elem_waves[i].rads.xx + elem2d.grad_shp[i].y * elem_waves[i].rads.xy);
//            wave_stress.y -= (elem2d.grad_shp[i].x * elem_waves[i].rads.xy + elem2d.grad_shp[i].y * elem_waves[i].rads.yy);
//        }
//    } else {
//        wave_stress.x = (1./3.) * (elem_waves[0].stress.x + elem_waves[1].stress.x + elem_waves[2].stress.x);
//        wave_stress.y = (1./3.) * (elem_waves[0].stress.y + elem_waves[1].stress.y + elem_waves[2].stress.y);
//        
//    }
//
//    // store stresses for node writing
//	for (i = 0; i<NDONTRI; i++) {
//    	assert( grid2d_nodeIDs[i] != UNSET_INT);
//        waves[ grid2d_nodeIDs[i] ].stress.x = wave_stress.x;
//        waves[ grid2d_nodeIDs[i] ].stress.y = wave_stress.y;
//	}
//    
//    return wave_stress;
//}
//
///***********************************************************/
///***********************************************************/
///***********************************************************/
//void swave_elem2d_get_local(SWAVE *waves, SWAVE *elem_waves, int *nodeIDs) {
//    int i;
//    for (i=0; i<NDONTRI; i++) {
//        elem_waves[i].rads.xx  = waves[nodeIDs[i]].rads.xx;
//        elem_waves[i].rads.yy  = waves[nodeIDs[i]].rads.yy;
//        elem_waves[i].rads.xy  = waves[nodeIDs[i]].rads.xy;
//        elem_waves[i].stress.x = waves[nodeIDs[i]].stress.x;
//        elem_waves[i].stress.y = waves[nodeIDs[i]].stress.y;
//    }
//}

