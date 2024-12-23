#ifndef H_TIME_LOOP_
#define H_TIME_LOOP_

/***********************************************************/
/***********************************************************/
/***********************************************************/
int update_dt(SMODEL_DESIGN *dm);
int time_loop(SMODEL_DESIGN *dm);
int advance_time(SMODEL_DESIGN *dm, int nsuper);
void set_forward_steppers(int (*forward_stepper[N_TIME_STEPPERS]) (SMODEL_SUPER*));
/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif
