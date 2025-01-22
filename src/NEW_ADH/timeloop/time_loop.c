/*! \file  time_loop.c This file has functions responsible for running an array of supermodels in time */
#include "adh.h"
//could move to global file? or maybe better to keep here
static int (*forward_stepper[N_TIME_STEPPERS]) (SMODEL_SUPER*);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function performs the time loop on a design model (array of super models)
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] dm (SMODEL_DESIGN*) - the design model where the supermodels are
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int time_loop(SMODEL_DESIGN *dm){
	int TIME_STEP_WORKED = TRUE, LOOP_INCOMPLETE = TRUE;
	int ierr = YES;
	int ts=0;
	//assign function pointers, maybe move to a global file?
	set_forward_steppers(forward_stepper);
	//aliasing for convenience
	int nsuper = dm->nSuperModels;
	//the outer time loop
	do {
			ts+=1;
			//update t_prev and change dt if timeseries says to
			LOOP_INCOMPLETE = update_dt(dm);
			//attempt to advance the whole design model by one time step
			TIME_STEP_WORKED = advance_time(dm, nsuper);
			//printf("TIME_STEP_WORKED = %d, LOOP_INCOMPLETE = %d\n",TIME_STEP_WORKED,LOOP_INCOMPLETE);

			//output step

			//calculate adaption

			//fluxes for external coupling 
			

	}while(LOOP_INCOMPLETE && TIME_STEP_WORKED);


	return ierr;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function advances each super model by calling the forward stepper of the
 *             super model
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] dm (SMODEL_DESIGN*) - the design model where the supermodels are
 *  @param [in] nsuper (int) - number of super models to call
 *  \returns - potentially an error code of whether the time step was succesfull or not
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int advance_time(SMODEL_DESIGN *dm, int nsuper){
	int TS_SUCCESS = TRUE;
	int nsubsteps;
	int isub;
	int isuper;
	SMODEL_SUPER *sm;
	//inner loop, loop over each supermodel
	//for now just simple time lagging from super to super
	//maybe in future can do Strang splitting
	for(isuper = 0; isuper<nsuper; isuper++){
		sm = &(dm->superModel[isuper]);
		nsubsteps = sm->nsubsteps;
		//advance one super model the proper number of steps forward
		//we restrict super models to have time steps as integer fractions
		//of the design model time step
		for (isub=0; isub<nsubsteps; isub++){
			//update sol_old and sol_older
			// prep solutions
 			smodel_super_prep_sol(sm);
			//update dirichlet data
			smodel_super_update_dirichlet_data(sm);

			//forward step in time
			TS_SUCCESS = smodel_super_forward_step(sm, forward_stepper[sm->forward_step]);
			//check if time step succeeded
			if (!TS_SUCCESS){
				return TS_SUCCESS;
			}
		}

		//exchange info between supermodels, need to think how to do this
		//exchange_info(dm,isuper);

	}
	return TS_SUCCESS;

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function updates time and potentially the dt based on time series input
 *             super model
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *  @param[in,out] dm (SMODEL_DESIGN*) - the design model where the supermodels are
 *  \returns - False if final time is passed. true if final time is not passed.
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int update_dt(SMODEL_DESIGN *dm){
	//advances time
	dm->t_prev+=dm->dt;

	//need to add a check and swap to new dt if necessary
	//dm->dt_old = dm->dt

	if(dm->t_prev >= dm->t_final){
		return FALSE;
	}

	return TRUE;
}
