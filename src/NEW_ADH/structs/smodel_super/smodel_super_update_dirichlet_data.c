#include "adh.h"
void smodel_super_update_dirichlet_data(SMODEL_SUPER *sm){
	//this is where all sseries stuff should go
	//use it to update sm->dirichlet_data
	double t = *(sm->t_prev);
	//for now just a simple hard code to get things working
	if (sm->elem2d_physics_mat[0].model[0].fe_resid == HEAT){
		double alpha = 3.0;
		double beta = 1.2;
		//u_{exact} = 1 + x^2 + \alpha y^2 + \beta t
		//printf("Updating heat boundary for T=%f\n",t);
		double x_coord, y_coord;
		int nnodes = sm->grid->nnodes;
		SGRID *grid = sm->grid;
		int id;
		for (int i=0; i<nnodes; i++){
			//mark the boundary only
			x_coord = grid->node[i].x;
			y_coord = grid->node[i].y;
			//id = grid->node[i].id;
			if ( is_near(x_coord,0.0) || is_near(x_coord,1.0) || is_near(y_coord,0.0) || is_near(y_coord,1.0) ){
				sm->dirichlet_data[i*3+1] = 1.0 + x_coord*x_coord + alpha * y_coord*y_coord + beta*t;
			}
		}	
	}
}
