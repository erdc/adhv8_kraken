/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  smodel_super.c This file collects methods of the SUPER_MODEL structure */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
int smodel_super_forward_step(SMODEL_SUPER* sm, int (*ts_fnctn)(SMODEL_SUPER*)){
	//a simple function wrapper
	return ts_fnctn(sm);
}
