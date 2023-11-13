#include "global_header.h"
void read_bc_AWRITE(SMODEL *mod, char *data) {
  
  char line[MAXLINE];           /* the input line */
  char *subdata = NULL;         /* the data after the second card   is read */
  char *subsubdata = NULL;

  int i, j, k;
  int id = 0, ientry = 0, awrite_nentries = 0;
  int uout = UNSET_INT, uin = UNSET_INT, isize = 0, infact = 0, outfact = 0;
  double newtime = 0;
  double *start=NULL; // awrite entry data
  double *end=NULL;   // awrite entry data
  double *dt=NULL;    // awrite entry data

  SIO info = *(mod->io);   // alias
 
  /*********************************************************************/
  /* begin reading in the data for this series *********************** */

  /* series id */
  id = sseries_read_id(info, &data, mod->nseries);

  /* series size */
  awrite_nentries = read_int_field(info, &data);
  start = (double *) tl_alloc(sizeof(double), awrite_nentries);
  end   = (double *) tl_alloc(sizeof(double), awrite_nentries);
  dt    = (double *) tl_alloc(sizeof(double), awrite_nentries);
  for (ientry = 0; ientry < awrite_nentries; ientry++) {
      start[ientry] = 0;
      end[ientry] = 0;
      dt[ientry] = 0;
  }

  /* series outfact */
  uout = read_int_field(info, &data);
  outfact = tc_conversion_factor(uout, FROM);

  /* loop over awrite entries */
  for (ientry = 0; ientry < awrite_nentries; ientry++) {
    fgets(line, MAXLINE, info.bc.fp);
    data = &line[0];
    start[ientry] = read_dbl_field(info, &data);
    end[ientry]   = read_dbl_field(info, &data);
    dt[ientry]    = read_dbl_field(info, &data);

    /* series infact */
    uin = read_int_field(info, &data);
    infact = tc_conversion_factor(uin, TO);
    start[ientry] *= infact;
    end[ientry]   *= infact;
    dt[ientry]    *= infact;

    /* if one entries end time is greater than the next start time, reset the next start time */
    if (ientry > 0 && (end[ientry - 1] >= start[ientry])) {
      start[ientry] = end[ientry - 1] + dt[ientry];
    }

    /* check input units */
    if (uin < SECONDS || uin > WEEKS) {
      io_read_error(info, "Invalid time units. Use a value between 0 and 4.", TRUE);
    }

    /* calculate size of automatically built output time-series */
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
    }
    else {
      newtime = series->entry[ientry - 1].time + awrite_dt;
    }
    if (newtime > awrite_final_end) {
      series->entry[ientry].time = awrite_final_end;
    }
    else if (newtime >= awrite_end) {
      series->entry[ientry].time = awrite_end;

      if (ientry < series->size - 1) {
        awrite_ientry++;
        flag = 1;
      }
    }
    else {
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
