/* This routine decides whether or not this is the last step and sets dt for the next step */

#include "adh.h"

int tc_end(SMODEL_DESIGN *mod)
{
    int interval;               /* the interval of the new time step */
    double dt_new;              /* the new delta t */

    /* check if we are done then return a YES */
    mod->t_prev += mod->dt;
    //mod->t_prev = messg_dmin(t_prev);

//    if (STEADY_STATE) {
//        if (myid <= 0)
//            printf("\n***STD: old_spatial_residual %f tol_nonlin %f\n", old_spatial_residual, tol_nonlin);
//        if (old_spatial_residual < tol_nonlin) {
//            if (myid <= 0)
//                printf("SOLUTION HAS REACHED STEADY STATE\n");
//            return (YES + STEADY_STATE);
//        }
//        if (t_prev >= tfinal) {
//            if (myid <= 0)
//                printf("WARNING: SOLUTION HAS REACHED TFINAL BEFORE ACHIEVING STEADY STATE\n");
//            return (YES);
//        }
//        return (NO);
//    }
//    if (AUTO_TIME_FIND) {
//        if (t_prev >= tfinal) {
//            return (YES);
//        }
//
//        interval = tc_find_interval(dt_sers, t_prev);
//        max_stdt = tc_eval_series(dt_sers, interval, t_prev);
//
//        /* if we are within one time step of the end, then take exactly the right amount */
//        if (t_prev + dt >= tfinal) {
//            dt = tfinal - t_prev;
//        }
//        /*if we are within one and a little time step of the end, then take half the time */
//        else if (t_prev + 2.0 * dt >= tfinal) {
//            dt = 0.5 * (tfinal - t_prev);
//        }
//
//        return (NO);
//    }

    //Mark, need to adapt this to new way
//    if (mod->t_prev >= mod->t_final) {
//        return (YES);
//    }//

//    /* if we failed on the last step then do not adjust the time step 
//       if we did not fail, then relax the time step */
//    if (mod->flag.TIME_ADAPT_FAIL == NO)
//        mod->dt *= DT_ENLARGE_FACTOR;//

//    /* initializes the new dt to the time series provided by the user */
//    interval = sseries_get_interval(*(mod->series_dt), mod->t_prev);
//    dt_new = tc_eval_series(*(mod->series_dt), interval, mod->t_prev, 0);//

//    /* compares to the time step indicated by the error indicator */
//    if (mod->t_adpt_flag == ON)
//        if (dt_new > mod->dt_err)
//            dt_new = mod->dt_err;//

//    /* compares to the current increment */
//    if (dt_new < mod->dt)
//        mod->dt = dt_new;//

//    /* if we are within one time step of the end, then take exactly the right amount */
//    if (mod->t_prev + mod->dt >= mod->t_final)
//        mod->dt = mod->t_final - mod->t_prev;
//    /* if we are within one and a little time step of the end, then take half the time */
//    else if (mod->t_prev + 2.0 * mod->dt >= mod->t_final)
//        mod->dt = 0.5 * (mod->t_final - mod->t_prev);//

//    /* makes sure everyone is on the same dt */
//#ifdef _MESSG
//    mod->dt = messg_dmin(mod->dt,mod->grid->smpi->ADH_COMM);
//#endif

    /* returns a no */
    return (NO);
}
