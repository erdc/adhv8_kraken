/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  read_bc_AWRITE This file collects methods to read an AdH SuperModel boundary condition input file       */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
static int DEBUG = OFF;
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Reads a SuperModel boundary condition file for automatic writing
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

void read_bc_AWRITE(SMODEL_SUPER *mod, char **token, FILE *fp) {
    
    int isers, ientry, uin, uout, infact, outfact;
    
    // series id
    int id = get_next_token_int(token) - 1;
    if (isers > mod->nseries - 1 || isers < 0) {
        io_read_error(io, "The specified series does not exist.", TRUE);
    }
    
    // series size
    int awrite_nentries = get_next_token_int(token);
    double *start=NULL; start = (double *) tl_alloc(sizeof(double), awrite_nentries);
    double *end=NULL;   end   = (double *) tl_alloc(sizeof(double), awrite_nentries);
    double *dt=NULL;    dt    = (double *) tl_alloc(sizeof(double), awrite_nentries);
    for (ientry = 0; ientry < awrite_nentries; ientry++) {
        start[ientry] = 0;
        end[ientry] = 0;
        dt[ientry] = 0;
    }
    
    // series output units
    uout = get_next_token_int(token);
    outfact = tc_conversion_factor(uout, FROM);
    
    // read awrite entries
    char *line = NULL;  size_t len = 0; ssize_t read;
    for (ientry = 0; ientry < awrite_nentries; ientry++) {
        read = getline(&line, &len, fp);
        cdum = get_token(line,&token);
        start[ientry] = get_next_token_dbl(&token);
        end[ientry]   = get_next_token_dbl(&token);
        dt[ientry]    = get_next_token_dbl(&token);
        
        // series input units
        uin = get_next_token_int(token);
        infact = tc_conversion_factor(uin, TO);
        start[ientry] *= infact;
        end[ientry]   *= infact;
        dt[ientry]    *= infact;
        
        // if one entries end time is greater than the next start time, reset the next start time
        if (ientry > 0 && (end[ientry - 1] >= start[ientry])) {
            start[ientry] = end[ientry - 1] + dt[ientry];
        }
        
        // check input units
        if (uin < SECONDS || uin > WEEKS) {
            io_read_error(info, "Invalid time units. Use a value between 0 and 4.", TRUE);
        }
        
        // calculate size of automatically built output time-series
        isize += (end[ientry] - start[ientry]) / dt[ientry] + 1;
        if (fmod(isize, 1) != 0) {
            isize += 1 - fmod(isize, 1);
        }
    }
    
    // add new series
    SSERIES *series;
    sseries_alloc(&series, (int) isize, OUTPUT_SERIES, UNSET_INT);
    series->id = id;
    series->size = (int) isize;
    series->outfact = outfact;
    series->infact = infact;
    
    /* create time-series from awrite */
    int awrite_final_entry = awrite_nentries - 1;
    double awrite_final_start = start[awrite_final_entry];
    double awrite_final_end = end[awrite_final_entry];
    int flag = 1;                 /* am I at the begginning of an awrite entry */
    int awrite_ientry = 0;
    double awrite_start, awrite_end, awrite_dt;
    
    for (ientry = 0; ientry < series->size; ientry++) {
        if (flag) {
            awrite_start = start[awrite_ientry];
            awrite_end = end[awrite_ientry];
            awrite_dt = dt[awrite_ientry];
            newtime = awrite_start;
            flag = 0;
        } else {
            newtime = series->entry[ientry - 1].time + awrite_dt;
        }
        if (newtime > awrite_final_end) {
            series->entry[ientry].time = awrite_final_end;
        } else if (newtime >= awrite_end) {
            series->entry[ientry].time = awrite_end;
            
            if (ientry < series->size - 1) {
                awrite_ientry++;
                flag = 1;
            }
        } else {
            series->entry[ientry].time = newtime;
            
        }
        //printf("entry: %d \t time: %f \n", ientry, series->entry[ientry].time);
    }
    
    /* calculate slopes */
    for (ientry = 1; ientry < series->size; ientry++) {
        series->entry[ientry - 1].slope[0] = (series->entry[ientry].value[0] - series->entry[ientry - 1].value[0]) / (series->entry[ientry].time - series->entry[ientry - 1].time);
    }
    
    /* calculate areas */
    for (ientry = 0; ientry < series->size - 1; ientry++) {
        series->entry[ientry].area[0] = tc_trap_area(series->entry[ientry].value[0], series->entry[ientry + 1].value[0], series->entry[ientry].time, series->entry[ientry + 1].time);
    }
    
    uin = uout = 0;
    series->infact = tc_conversion_factor(uin, TO);
    
    //sseries_sers_print(series, 1);
    
    /* check series-series */
    sseries_check(*series);
    
    /* define as AdH output series  */
    mod->series_out = series;
    
    /* SET MODEL FLAGS */
    if (mod->flag.OUTPUT == ON)
        fprintf(stderr, "WARNING: Overwriting previous output series with output series %d!\n", series->id);
    
    /* assign model variables */
    mod->flag.OUTPUT = ON;
    
    // free variables
    start = (double *) tl_free(sizeof(double), awrite_nentries, start);
    end = (double *) tl_free(sizeof(double), awrite_nentries, end);
    dt = (double *) tl_free(sizeof(double), awrite_nentries, dt);
    
}
