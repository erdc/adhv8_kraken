#include "adh.h"


int time_loop(SMODEL_DESIGN *dm){
	int TIME_STEP_WORKED = TRUE, LOOP_INCOMPLETE = TRUE;
	int ierr = YES;
	int ts=0;

	//aliasing for convenience
	int nsuper = dm->nSuperModels;
	//the outer time loop
	do {
			ts+=1;
			printf("DESIGN MODEL ON TIME STEP %d\n",ts);
			//update t_prev and change dt if timeseries says to
			LOOP_INCOMPLETE = update_dt(dm);

			//attempt to advance the whole design model by one time step
			TIME_STEP_WORKED = advance_time(dm, nsuper);

			printf("TIME_STEP_WORKED = %d, LOOP_INCOMPLETE = %d\n",TIME_STEP_WORKED,LOOP_INCOMPLETE);

			//output step

			//calculate adaption

			//fluxes for external coupling 
			

	}while(LOOP_INCOMPLETE && TIME_STEP_WORKED);


	return ierr;
}


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
			TS_SUCCESS = sm->forward_step(sm);
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