#include "adh.h"
void smodel_super_prep_sol(SMODEL_SUPER *sm){
	int idof;
	int ndof = *(sm->ndofs);
    for (idof = 0; idof < ndof; idof++) {
      sm->sol_older[idof] = sm->sol_old[idof];
      sm->sol_old[idof] = sm->sol[idof];
    }
}
