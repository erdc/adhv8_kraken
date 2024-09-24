/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_SERIES.c  This file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file for reading a model series
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] *sm (SUPER_MODEL *)  a double pointer to an array of AdH supermodels
 * @param[in] **token (CHAR) a BC file line string token
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void read_bc_SERIES(SMODEL_SUPER *mod, char **token, , FILE *fp) {
    
    int i, j, k;
    int ientry = 0, ivalue = 0, nentries = 0, uout = UNSET_INT, uin = UNSET_INT, isize = 0, infact = 0;
    double newtime = 0;
    
    int id = get_id(token,mod->nseries,"SERIES ID Error\n");
    int size = get_next_token_int(token);
    
    SSERIES *series;
    
    // note: the grid nodes is for node contribution, always for 2d
    sseries_alloc(&series, size, type, mod->grid->nnodes_sur);
    series->id = id;
    series->nnodes = mod->grid->nnodes_sur;
    
    if (type == WAVE_SERIES || type == WIND_SERIES) {
        //smeteor_station_alloc(&(series->station), mod->grid->nnodes_sur);
        series->station->x = get_next_token_dbl(token);
        series->station->y = get_next_token_dbl(token);
    };
    
    if (type == OUTPUT_SERIES) {
        mod->flag.OUTPUT = ON;
    }
    
    /* input units */
    uin = get_next_token_int(token);
    if (uin < SECONDS || uin > WEEKS) {
        io_read_error(info, "Invalid time units. Use a value between 0 and 4.", TRUE);
    }
    series->infact = tc_conversion_factor(uin, TO);
    
    /* output units */
    uout = get_next_token_int(token);
    if (uout < SECONDS || uout > WEEKS) {
        io_read_error(info, "Invalid time units. Use a value between 0 and 4.", TRUE);
    }
    series->outfact = tc_conversion_factor(uout, FROM);
    
    /* read series from bc file */
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    read = getline(&line, &len, fp));
    get_token(line,&token); sscanf((*token), "%lf", &series->entry[0].time);
    series->entry[0].time *= series->infact;
    for (ivalue=0; ivalue<series->nvalues; ivalue++) {
        series->entry[0].value[ivalue] = get_next_token_dbl(token);
    }
    
    for (ientry = 1; ientry < series->size; ientry++) {
        read = getline(&line, &len, fp));
        get_token(line,&token); sscanf((*token), "%lf", &series->entry[0].time);
        series->entry[0].time *= series->infact;
        if (series->entry[ientry].time < series->entry[ientry - 1].time) {
            io_read_error(info, "series1 series is not in the correct order.", TRUE);
        }
        
        if (series->entry[ientry].time - series->entry[ientry - 1].time < SMALL) {
            io_read_error(info, "series1 series has same time for two values.", TRUE);
        }
        
        for (ivalue=0; ivalue<series->nvalues; ivalue++) {
            series->entry[ientry].value[ivalue] = get_next_token_dbl(token);
        }
    }
    
    for (ivalue=0; ivalue<series->nvalues; ivalue++) {
        /* calculate slopes */
        for (ientry = 1; ientry < series->size; ientry++) {
            series->entry[ientry - 1].slope[ivalue] = (series->entry[ientry].value[ivalue] - series->entry[ientry - 1].value[ivalue]) / (series->entry[ientry].time - series->entry[ientry - 1].time);
        }
        /* calculate areas */
        for (ientry = 0; ientry < series->size - 1; ientry++) {
            series->entry[ientry].area[ivalue] = tc_trap_area(series->entry[ientry].value[ivalue], series->entry[ientry + 1].value[ivalue], series->entry[ientry].time, series->entry[ientry + 1].time);
        }
    }
    
    // sanity check the series
    sseries_check(*series);
    
    // print the series to screen
    //sseries_sers_print(*series, 1);
    
    if (type == DT_SERIES) {
        mod->series_dt = series;
    }
    else if (type == OUTPUT_SERIES) {
        mod->series_out = series;
    }
    else if (type == WIND_SERIES) { // add to wind series linked list
        mod->flag.WIND = ON;
        mod->flag.WIND_STATION = ON;
        sseries_add(series, &(mod->series_wind_head), &(mod->series_wind_curr), TRUE);
        sseries_free(series);
    }
    else if (type == WAVE_SERIES) { // add to wave series linked list
        //sseries_printScreen(*series, 1);
        
        mod->flag.WAVE = ON;
        mod->flag.WAVE_STATION = ON;
        sseries_add(series, &(mod->series_wave_head), &(mod->series_wave_curr), TRUE);
        sseries_free(series);
        
    }
#ifdef _ADH_GROUNDWATER
    else if (type == CONSTITUITIVE_SERIES) { /* add to relativer permeability list */
        sseries_add(series, &(mod->series_gw_psk_head), &(mod->series_gw_psk_curr), TRUE);
        sseries_free(series);
    }
#endif
    else {
        // add to linked list (which will allocate another, so we can delete this one)
        sseries_add(series, &(mod->series_head), &(mod->series_curr), TRUE);
        sseries_free(series);
    }
}
