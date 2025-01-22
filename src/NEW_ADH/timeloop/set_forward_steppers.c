/*! \file  set_forward_steppers.c This file has sets global array of function pointers to time
 * stepping algorithms */
#include  "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function sets the global array of pointers to time step routines
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] forward_stepper (int (* [N_TIME_STEPPERS])(SMODEL_SUPER *)) - pointer
 *  to functions that advance time of a super model by one time step
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void set_forward_steppers(int (*forward_stepper[N_TIME_STEPPERS]) (SMODEL_SUPER*)){
	forward_stepper[0] = fe_newton;
}
