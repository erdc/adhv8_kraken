/* allocates and initializes with defaults the ts-sers */

// NOTE: THE SMODEL_SUPER_SANITY SHOULD CHECK FOR DUPLICATED SERIES!!
// do this by looking at series ids

#include "adh.h"
#include <stdbool.h>
#define BIG -9999999.0

/***********************************************************/
/***********************************************************/
/* series entry initialization */
void sseries_entry_init(SSERIES_ENTRY *entry) {
  entry->time = 0.0;
  entry->area = NULL;
  entry->slope = NULL;
  entry->value = NULL;
}

/***********************************************************/
/***********************************************************/
/* series initialization */
void sseries_init(SSERIES *series) {
  if(series == NULL) {
    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
    tl_error(">> trying to initialize a series that is not allocated!");
  }

  // default value for ts-sers
  series->id = -1;
  series->size = 0;
  series->type = -1;
  series->mat_id = -1;
  series->nvalues = 0;
  series->nnodes = UNSET_INT;
  series->infact = 1.0;
  series->outfact = 1.0;
  series->tol = XY_DEFAULT_TOL;
  series->ivalue = NULL;
  series->entry = NULL;
  series->station = NULL;
}

/***********************************************************/
/***********************************************************/
/* series allocation */
void sseries_alloc(SSERIES **series, int nentries, int type, int nnodes) {
  int i = 0, j = 0;

  // allocate the time-series
  (*series) = (SSERIES *) tl_alloc(sizeof(SSERIES), 1);
  sseries_init(*series);

  // fill in what we know from call arguements
  (*series)->size = nentries;
  (*series)->type = type;

  // allocate the entries of the time-series
  (*series)->nvalues = 1;
  if(type == WIND_SERIES) {
    (*series)->nvalues = 2;	// tauX, tauY
  }
  else if (type == WAVE_SERIES) {
    (*series)->nvalues = 3;	// Sxx, Sxy, Syy
  }

  // allocate interpolated value
  (*series)->ivalue = (double *)tl_alloc(sizeof(double), (*series)->nvalues);
  for (i=0; i<(*series)->nvalues; i++) {
      (*series)->ivalue[i] = 0.;
  }

  // allocate entries
  (*series)->entry = (SSERIES_ENTRY *) tl_alloc(sizeof(SSERIES_ENTRY), nentries);

  // if this is a meteor station (node contibutions done later)
  if (type == WIND_SERIES || type == WAVE_SERIES) {
    smeteor_station_alloc(&((*series)->station), nnodes);
  }

  // allocate value in entry and initialize entry
  for (i = 0; i < nentries; i++) {
    sseries_entry_init(&((*series)->entry[i]));
    (*series)->entry[i].area  = (double *)tl_alloc(sizeof(double), (*series)->nvalues);
    (*series)->entry[i].slope = (double *)tl_alloc(sizeof(double), (*series)->nvalues);
    (*series)->entry[i].value = (double *)tl_alloc(sizeof(double), (*series)->nvalues);
    for(j = 0; j < (*series)->nvalues; j++) {
         (*series)->entry[i].area[j] = 0.;
         (*series)->entry[i].slope[j] = 0.;
         (*series)->entry[i].value[j] = 0.;
    }
  }
}

/***********************************************************/
/***********************************************************/
/* series deallocation */
void sseries_free(SSERIES *series) {
  int i = 0;

  // free value in entry
  for(i = 0; i < series->size; i++) {
      if (series->entry[i].value != NULL) {
          series->entry[i].area  = (double *)tl_free(sizeof(double), series->nvalues, series->entry[i].area);
          series->entry[i].slope = (double *)tl_free(sizeof(double), series->nvalues, series->entry[i].slope);
          series->entry[i].value = (double *)tl_free(sizeof(double), series->nvalues, series->entry[i].value);
      }
  }

  // free interpolated value
  if(series->ivalue != NULL) {
    series->ivalue = (double *)tl_free(sizeof(double), series->nvalues, series->ivalue);
  }

  // free station
  if(series->station != NULL) {
    smeteor_station_free(series->station, series->nnodes);
  }

  // free entries
  if(series->entry != NULL) {
    series->entry = (SSERIES_ENTRY *) tl_free(sizeof(SSERIES_ENTRY), series->size, series->entry);
  }

  // now free time-series
  if(series != NULL) {
    series = (SSERIES *) tl_free(sizeof(SSERIES), 1, series);
  }

  series = NULL;
}

/***********************************************************/
/***********************************************************/
/* perform deep copy of series */
void sseries_cpy(SSERIES *dest, SSERIES *source) {

  int i = 0, ientry = 0;

  if(source == NULL) {
    dest = NULL;
    printf("SOURCE IS NULL (file:line) %s:%d\n", __FILE__, __LINE__);
    return;
  }
  dest->id = source->id;
  dest->size = source->size;
  dest->type = source->type;
  dest->mat_id = source->mat_id;
  dest->nvalues = source->nvalues;
  dest->nnodes = source->nnodes;
  dest->infact = source->infact;
  dest->outfact = source->outfact;
  dest->tol = source->tol;
  dest->next = source->next;
  dest->prev = source->prev;

  // later put in deep copy of station
  if (source->station != NULL) {
    smeteor_station_copy(dest->station,source->station, source->nnodes);
  }

  for(i = 0; i < source->nvalues; i++) {
    dest->ivalue[i] = source->ivalue[i];
  }

  for(ientry = 0; ientry < source->size; ientry++) {
    sseries_copy_entry(&(dest->entry[ientry]), source->entry[ientry], source->nvalues);
  }
}

/***********************************************************/
/***********************************************************/
/* deep copy of series entry */
void sseries_copy_entry(SSERIES_ENTRY *dest, SSERIES_ENTRY source, int nvalues) {
  dest->time = source.time;
  int i = 0;
  for(i = 0; i < nvalues; i++) {
      dest->area[i] = source.area[i];
      dest->slope[i] = source.slope[i];
      dest->value[i] = source.value[i];
  }
}

/***********************************************************/
/***********************************************************/
/* sanity check new ts-series */
void sseries_check(SSERIES series) {

  int ientry;

  // cjt :: since series is passed by copy, it is never null
  //if(&series == NULL) {
  //  printf("\nSERIES %d ERROR\n", series.id);	//sseries_print(series);
  //  printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
  //  tl_error(">> trying to initialize a series that is not allocated!");
  //}

  if(series.size <= 0) {
    printf("\nSERIES %d ERROR\n", series.id);
    sseries_printScreen(series, 0);
    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
    tl_error(">> series must have at least one entry!");
  }

  if(series.type < ANY_SERIES || series.type > WAVE_SERIES) {
    printf("\nSERIES %d ERROR\n", series.id);
    sseries_printScreen(series, 0);
    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
    tl_error(">> series type is not supported!");
  }

/*  if(series.infact < SECONDS || series.infact > WEEKS) {
    printf("\nSERIES %d ERROR\n", series.id);
    sseries_printScreen(series, 0);
    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
    printf("series.infact is %f but must be between %d and %d \n", series.infact, SECONDS, WEEKS);
    tl_error(">> series input time units are invalid!");
  }

  if(series.outfact < SECONDS || series.outfact > WEEKS) {
    printf("\nSERIES %d ERROR\n", series.id);
    sseries_printScreen(series, 0);
    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
    tl_error(">> series output time units are invalid!");
  }
*/
  if(series.nvalues < 1 || series.nvalues > 3) {
    printf("\nSERIES %d ERROR\n", series.id);
    sseries_printScreen(series, 0);
    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
    tl_error(">> number of value should be between 1 and 3");
  }

  if(series.size == 1)
    return;

  for(ientry = 1; ientry < series.size; ientry++) {

    if(series.entry[ientry].time < series.entry[ientry - 1].time) {
      printf("\nERROR: SERIES %d :: ENTRY: %d\n", series.id, ientry);
      sseries_printScreen(series, 0);
      printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
      tl_error(">> series entries are not in the correct order!");
    }
    if(series.entry[ientry].time - series.entry[ientry - 1].time < SMALL) {
      printf("\nERROR: SERIES %d :: ENTRY: %d\n", series.id, ientry);
      sseries_printScreen(series, 0);
      printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
      tl_error(">> series has the same time for two entries!");
    }

  }
}

/***********************************************************/
/***********************************************************/
/* check the series list */
void sseries_checklist(int nseries, SSERIES *head) {

  int sers_flag[nseries];
  SSERIES *ptr = head;

  int isers = 0;

  /* check for missing and duplicate ids */
  for(isers = 0; isers < nseries; isers++)
    sers_flag[isers] = 0;
  while(ptr != NULL) {
    sers_flag[ptr->id] += 1;
    ptr = ptr->next;
  }

  for(isers = 0; isers < nseries; isers++) {
    if(sers_flag[isers] > 1) {
      printf("\nSERIES %d DUPLICATED\n", isers);
      printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
      tl_error(">> series duplicated!");
    }
    /*
    if(sers_flag[isers] == 0) {
      printf("\nSERIES %d NOT INCLUDED\n", isers);
      printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
      tl_error(">> -series skipped!");
    }
    */
  }
}

/***********************************************************/
/***********************************************************/
/* print an ts-series, full means actually print all entries as well */
void sseries_printScreen(SSERIES ts, int full){
  int ientry;
  printf("------------------\n");
  printf("series: id: %d  :: type :: ", ts.id + 1);

  if(ts.type == WIND_SERIES) {
    printf("wind series\n");
  }
  if(ts.type == WAVE_SERIES) {
    printf("wave series\n");
  }
  else if(ts.type == TIME_SERIES) {
    printf("time series\n");
  }
  else if(ts.type == OUTPUT_SERIES) {
    printf("output series\n");
  }
  else if(ts.type == CONSTITUITIVE_SERIES) {
    printf("constituitive series\n");
  }
  else if(ts.type == ANY_SERIES) {
    printf("any series\n");
  }
  else if(ts.type == DT_SERIES) {
    printf("time-step series\n");
  }

  printf("size: %d  ", ts.size);
  printf("infact: %f  outfact: %f  ", ts.infact, ts.outfact);
  printf("nvalues: %d nsurface_nodes: %d\n",ts.nvalues,ts.nnodes);

  if(ts.type == WIND_SERIES) {
      printf("-- station location: {%f,%f} \n", ts.station->x, ts.station->y);
      printf("current time-interpolated value:  wind_stress.x: %10.5f \t wind_stress.y: %10.5f\n", ts.ivalue[0], ts.ivalue[1]);
  }
  else if(ts.type == WAVE_SERIES) {
      printf("-- station location: {%f,%f} \n", ts.station->x, ts.station->y);
      printf("current time-interpolated value:  wave_rads.xx: %10.5f \t wave_rads.xy: %10.5f \t wave_rads.yy: %10.5f\n",ts.ivalue[0], ts.ivalue[1], ts.ivalue[2]);
  }
  else {
    printf("current time-interpolated value: %10.5f\n", ts.ivalue[0]);
  }

  if(full) {
    for(ientry = 0; ientry < ts.size; ientry++) {
      printf("entry: %d \t", ientry);
      printf("time: %10.5f \t", ts.entry[ientry].time);

      if(ts.type == WIND_SERIES) {
          printf("wind_stress.x: %10.5f \t wind_stress.y: %10.5f", ts.entry[ientry].value[0], ts.entry[ientry].value[1]);
      }
      else if(ts.type == WAVE_SERIES) {
          printf("wave_rads.xx: %10.5f \t wave_rads.xy: %10.5f \t wave_rads.yy: %10.5f",ts.entry[ientry].value[0], ts.entry[ientry].value[1], ts.entry[ientry].value[2]);
      }
      else {
          printf("value: %10.5f", ts.entry[ientry].value[0]);
      }

      if(ts.type != OUTPUT_SERIES)
          printf("\t area: %10.5f \t slope: %10.5f\n", ts.entry[ientry].area[0], ts.entry[ientry].slope[0]);
      else
          printf("\n");
    }
  }
}

/***********************************************************/
/***********************************************************/
/* print the entire list of ts-series, including entries */
void sseries_printScreen_list(int full, SSERIES *head) {

  SSERIES *ptr = head;
  while(ptr != NULL) {
    sseries_printScreen(*ptr, full);
    ptr = ptr->next;
  }
  return;
}

/***********************************************************/
/***********************************************************/
/* create ts linked list given the first ts-series */
SSERIES *sseries_create_list(SSERIES *series, SSERIES **head, SSERIES **curr) {

#ifdef _ts_VERBOSE
  printf("\n creating series linked list \n");
#endif

  /* allocate new ts-series node */
  SSERIES *ptr;
  sseries_alloc(&ptr, series->size, series->type, series->nnodes);  

  if(ptr == NULL) {
    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
    tl_error(">> series linked list creation failed.");
  }
  sseries_cpy(ptr, series);
  ptr->next = NULL;
  *head = *curr = ptr;
  return ptr;
}

/***********************************************************/
/***********************************************************/
/* add an ts-series to the linked list */
SSERIES *sseries_add(SSERIES *series, SSERIES **head, SSERIES **curr, bool add_to_end) {

  if(*head == NULL) {
    return (sseries_create_list(series, head, curr));
  }

#ifdef _ts_VERBOSE
  if(add_to_end)
    printf("\n Adding node to end of the series\n");
  else
    printf("\n Adding node to beginning of the series list\n");
#endif

  /* allocate new ts-series node */
  SSERIES *ptr;
  sseries_alloc(&ptr, series->size, series->type, series->nnodes);
  if(ptr == NULL) {
    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
    tl_error(">> series linked list addition failed.");
  }
  sseries_cpy(ptr, series);
  ptr->next = NULL;
  if(add_to_end) {
    (*curr)->next = ptr;
    *curr = ptr;
  }
  else {
    ptr->next = *head;
    *head = ptr;
  }
  return ptr;
}

/***********************************************************/
/***********************************************************/
/* return the series given id */
SSERIES *sseries_search(int id, SSERIES ** prev, SSERIES *head) {

  SSERIES *ptr = head;
  SSERIES *tmp = NULL;
  bool found = false;

#ifdef _ts_VERBOSE
  printf("\n Searching the series list for id: [%d] \n", id);
#endif

  while(ptr != NULL) {
    if(ptr->id == id) {
      found = true;
      break;
    }
    else {
      tmp = ptr;
      ptr = ptr->next;
    }
  }

  if(true == found) {
    if(prev)
      *prev = tmp;
    return ptr;
  }
  else {
    return NULL;
  }
}

/***********************************************************/
/***********************************************************/
/* delete a given series from the list */
int sseries_delete(int id, SSERIES **head, SSERIES **curr) {

  SSERIES *prev = NULL;
  SSERIES *del = NULL;

//#ifdef _ts_VERBOSE
  printf("\n Deleting series [%d] from list\n", id);
//#endif

  del = sseries_search(id, &prev, *head);
  if(del == NULL) {
    return -1;
  }
  else {
    if(prev != NULL)
      prev->next = del->next;

    if(del == *curr) {
      *curr = prev;
    }
    else if(del == *head) {
      *head = del->next;
    }
  }

  sseries_free(del);
  return 0;
}

/***********************************************************/
/***********************************************************/
/* free the entire list of series */
void sseries_free_list(SSERIES **t_head) {

  SSERIES *head = *(t_head);

#ifdef _ts_VERBOSE
  printf("... freeing ts-series list ...\n");
#endif

  SSERIES *temp = NULL;
  while(head != NULL) {
    temp = head->next;
    sseries_free(head);
    head = temp;
  }

  return;
}

/***********************************************************/
/***********************************************************/
/* get a value of a series */
double sseries_get_value(int id, SSERIES *head, int index) {

  assert(head);
  SSERIES *ts = head;
  while(ts != NULL) {
    if(ts->id == id) {
      return ts->ivalue[index];
    }
    ts = ts->next;
  }

  // should never get here
  printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
  printf("Series %d is bad || type: %d \n",id+1,ts->type);
  tl_error(">> could not find series id in get value.");

  return UNSET_FLT;
}

/***********************************************************/
/***********************************************************/
/* set the interpolated values of a series to a given time*/
void sseries_set_values(SSERIES *series, double time) {
    int i;
    int interval = sseries_get_interval(*series, time);
    for (i=0; i<series->nvalues; i++) {
        series->ivalue[i] = tc_eval_series(*series, interval, time, i); 
    }
}

/***********************************************************/
/***********************************************************/
/* set the interpolated values of a series to a given time for all time-series*/
void sseries_setall_ts_values(SSERIES *head, double time) {
    SSERIES *series = head;
    while (series != NULL) {
        if (series->type == TIME_SERIES) {
            sseries_set_values(series, time);
        }
        series = series->next;
    }
    series = NULL;
}

/***********************************************************/
/***********************************************************/
/* set the interpolated values of a series as the average of two times*/
void sseries_set_values_average(SSERIES *series, double t1, double t2) {

    int i;                      /* loop counter */
    int interval1, interval2;   /* the left hand endpoints of the intervals  containing t1  and t2 */
    int min_trap;               /* right endpoint of interval containing t1 - first data    point contained in (t1,t2) */
    int max_trap;               /* right endpoint of interval containing t2 -  one more     endpoint than the last contained in (t1,t2) */
    double tot_area;            /* total area from t to t+dt */
    double value1, value2;      /* the values of the two points */
    
    /* calculate where t1 and t2 are in the series series */
    interval1 = sseries_get_interval(*series, t1);
    interval2 = sseries_get_interval(*series, t2);
    min_trap = interval1 + 1;
    max_trap = interval2 + 1;

    int index;
    for (index = 0; index < series->nvalues; index++) {
        value1 = tc_eval_series(*series, interval1, t1, index);
        value2 = tc_eval_series(*series, interval2, t2, index);

        /* initializes the total area under the time series to zero */
        tot_area = 0.0;

        /* initializes the total area under the time series to zero */
        tot_area = 0.0;

        /* splits into two cases */
        /* both times are in the same interval */
        if (interval1 == interval2) { /* add the area of the partial trapezoid between the two points */
            tot_area += tc_trap_area(value1, value2, t1, t2);
        }
        else { /* they are in different intervals */
            /* calculate the area full trapezoids contribute */
            for (i = min_trap; i < (max_trap - 1); i++)
                tot_area += series->entry[i].area[index];

            /* add the area for the first partial trapezoid */
            tot_area += tc_trap_area(value1, series->entry[interval1 + 1].value[index], t1, series->entry[interval1 + 1].time);

            /* add the area for the second partial trapezoid */
            tot_area += tc_trap_area(series->entry[interval2].value[index], value2, series->entry[interval2].time, t2);
        }

        /* sets the value to the average value over the interval (t1,t2) */
        series->ivalue[index] = tot_area / (t2 - t1);
    }
}

/***********************************************************/
/***********************************************************/
/* set the interpolated values of a series as average  time for all time-series*/
void sseries_setall_ts_valuesAVG(SSERIES *head, double t1, double t2) {
    SSERIES *series = head;
    while (series != NULL) {
        if (series->type == TIME_SERIES) {
            sseries_set_values_average(series, t1, t2);
        }
        series = series->next;
    }
    series = NULL;
}

/***********************************************************/
/***********************************************************/
/* extra a series from the linked list */
SSERIES *sseries_extract(SSERIES *head, int id) {
    SSERIES *series;
    series = sseries_search(id, NULL, head);
    if (series == NULL) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> dt series not found!");
    } 

    // allocated new series
    SSERIES *series_new;
    sseries_alloc(&series_new, series->size, series->type, series->nnodes);
 
    // copy list series into stand-alone
    sseries_cpy(series_new, series);

    // delete list series
    sseries_delete(id, &(head), &(head));
    series = NULL;

    /* not sure why this is here ...   
    series = head;
    while (series!= NULL) {
        printf("series id: %d \n",series->id);
        series = series->next;
    }
    */
    series = NULL;

    return series_new;
}

/***********************************************************/
/***********************************************************/
int sseries_read_id(SIO io, char **data, int nseries) {
  /* gets and checks the series number */
  int isers = read_int_field_custom(io, data, NULL, "series ID", UNSET_INT, TRUE);
  if (isers > nseries || isers < 1) {
    io_read_error(io, "The specified series does not exist.", TRUE);
  }
  isers--;

  return (isers);
}

/***********************************************************/
/***********************************************************/
int sseries_set_type(SMODEL_SUPER *mod, char **subdata, int type) {
    SSERIES *sers = NULL;
    int iseries = sseries_read_id(*(mod->io), subdata, mod->nseries);
    sers = sseries_search(iseries, NULL, mod->series_head);
    if (sers == NULL) {
        printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
        tl_error(">> cannot find series id! (Check series numbers in BC strings - AWRITE is not included in list) \n");
    }
    sers->type = type;
    return iseries;
}

/***********************************************************/
/***********************************************************/
/* find the time-series interval the given time falls in */
int sseries_get_interval(SSERIES series, double time) {

  int i;                        /* loop counter */
  int ilast;                    /* the last acceptable i */
  double min_t, max_t;          /* the minimum and maximum t */

  /* set the minimum and maximum t in the series */
  min_t = series.entry[0].time * (1.0 - SMALL);
  ilast = series.size - 1;
  max_t = series.entry[ilast].time * (1.0 + SMALL);

  /* check if t1 and t2 are in the interval */
  if (time < min_t || time > max_t) {
    printf("Error (file:line) %s:%d\n", __FILE__, __LINE__);
    fprintf(stderr, "Series:  %d\n", series.id + 1);
    fprintf(stderr, "Time:  %14.6e\n", time);
    tl_error("Time is not in the range of the series.");
  }

  /* finds the interval containing the time */
  i = 0;
  while (i < ilast && time >= series.entry[i].time * (1.0 - SMALL)) {
    i++;
  }
  i--;

  return (i);
}

/***********************************************************/
/***********************************************************/
int series_get_type_count(SSERIES *head, int type) {
    int count = 0;   
    SSERIES *series = head;
    while (series != NULL) {
        if (series->type == type) {
            count++;
        }
        series = series->next;
    } 
    return count;
}
/***********************************************************/
/***********************************************************/
void sseries_set_meteor_stations(SSERIES *head, SGRID *grid, int type) {
    
// CJT -- fix this later!

    int i=UNSET_INT, id3d=UNSET_INT;

    // calculate nodal contributions for interpolation
    int nstations =  series_get_type_count(head, type);

    double *x = (double *) tl_alloc(sizeof(double), grid->nnodes_bed);
    double *y = (double *) tl_alloc(sizeof(double), grid->nnodes_bed);
    if (grid->ndim == 2) {
        for (i=0; i<grid->nnodes_bed; i++) {
            x[i] = grid->node[i].x;
            y[i] = grid->node[i].y;
        }
    } else if (grid->ndim == 3) {
        for (i=0; i<grid->nnodes_bed; i++) {
            id3d = grid->nodeID_2d_to_3d_sur[i];
            x[i] = grid->node[id3d].x;
            y[i] = grid->node[id3d].y;
        }
    }

    //station_node_contrib(head, x, y, nstations);


    x = (double *) tl_free(sizeof(double), grid->nnodes_bed, x);
    y = (double *) tl_free(sizeof(double), grid->nnodes_bed, y);
}

/***********************************************************/
/***********************************************************/
// this is for a full list and more general - making it different than set_value
void series_list_update_values(SSERIES *series_head, double t_prev) {
        int interval, indx;

        SSERIES *series;
        series = series_head;
        while (series != NULL) {
            interval = sseries_get_interval(*series, t_prev);
            for (indx = 0; indx < series->nvalues; indx++) {
                series->ivalue[indx] = tc_eval_series(*series, interval, t_prev, indx);
            }
            series = series->next;
        }
        series = NULL;
}

