/*! \file timeloop_test.c This file tests the timeloop via heat equation */
#include "adh.h"
static double NEWTON_TEST_TOL = 1e-7;
static int NEWTON_TEST_NX = 100;//150;
static int NEWTON_TEST_NY = 100;//150;
static double alpha = 3;
static double beta = 1.2;
static void compute_exact_solution_heat(double *u_exact, int ndof, SGRID *grid, double t);
static void permute_array(double *arr,int *p, int n);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function tests the time looper and Newton solver
 *  using a heat equation with analytic
 *  solution
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int timeloop_test(int argc, char **argv) {
	//create a grid
	SGRID *grid;
	grid = (SGRID *) tl_alloc(sizeof(SGRID), 1);
	//let's just do something simple
	//3x3 triangular element grid
	double xmin = 0.0;
	double xmax = 1.0;
	double ymin = 0.0;
	double ymax = 1.0;
	int npx = NEWTON_TEST_NX;
	int npy = NEWTON_TEST_NY;
	double theta = 0.0;
	double dz = 1.0;
	double a0 = -5.0;
	double ax = 0.0;
	double ax2 = 0.0;
	double ay = 0.0;
	double ay2 = 0.0;
	double axy = 0.0;
	double ax2y = 0.0;
	double axy2 = 0.0;
	double ax2y2 = 0.0;
	int flag3d =0;
	int err_code=-1;
    *grid = create_rectangular_grid(xmin, xmax, ymin, ymax, npx, npy,
 	theta, dz, a0, ax, ax2, ay, ay2, axy,
    ax2y, axy2, ax2y2, flag3d );
    int nnodes;
    nnodes = grid->nnodes;
    sgrid_reorder(grid,2);
	//print coordinates
//  for(int local_index =0; local_index<grid.nnodes; local_index++){
//		printf("Node %d: (x,y) = {%f,%f}\n",grid.node[local_index].gid,grid.node[local_index].x,grid.node[local_index].y);
//	}
//	//print connectivity
//	for(int local_index =0; local_index<grid.nelems2d; local_index++){
//		printf("Element %d: (nd1,nd2,nd3) = {%d,%d,%d}\n",local_index ,grid.elem2d[local_index].nodes[0], grid.elem2d[local_index].nodes[1], grid.elem2d[local_index].nodes[2]);
//	}


	//create a supermodel with a grid in it
	//SMODEL_SUPER sm;
	//sm.grid = &grid;
    SMODEL_DESIGN dm;
	//specify elemental physics and other properties in super model
	double dt = 2.0/20.0;
	double t0 = 0.0;
	double tf = 2.0;

	char elemVarCode[4]; 
	strcpy(&elemVarCode[0],"2");//SW2D
	strcpy(&elemVarCode[1],"0"); //GW
	strcpy(&elemVarCode[2],"0"); //Transport
	printf("GRID NELEMS2D = %d\n",grid->nelems2d);
	//smodel_super_no_read_simple(&sm, dt, t0, tf, 0 , 1, 0, elemVarCode);
	smodel_design_no_read_simple(&dm, dt, t0, tf,0, 1, 0, elemVarCode, grid);
	printf("NDOFS %d\n",dm.ndofs[0]);

	//OVER WRITE TO HEAT
	dm.superModel[0].elem2d_physics_mat[0].model[0].fe_resid = HEAT;
	//printf("Supermodel no read complete\n");


	//allocate linear system
	//doesn't currently work, need to go back and fix
	//fe_allocate_initialize_linear_system(&sm);
//	sm.cols_diag = NULL;
//	sm.vals_diag = NULL;
//	sm.cols_off_diag=NULL;
//	sm.vals_off_diag=NULL;
//	sm.bc_mask = NULL;
//	sm.nnz_diag_old=0;
//	sm.nnz_off_diag_old=0;
//	printf("Calling sparsity split CSR\n");
//	create_sparsity_split_CSR(&sm, sm.grid);
	//Screen_print_CSR(sm.indptr_diag, sm.cols_diag, sm.vals_diag, sm.ndofs);

	//do we want to stor nnz? it is stored in sm->indptr[nrows]
    //do we want to store local_size = local_range[1]-local_range[0]
//    printf("NNZ = %d = %d\n",sm.indptr_diag[sm.my_ndofs], sm.nnz_diag);
//    sarray_init_dbl(sm.vals_diag, sm.indptr_diag[sm.my_ndofs]);//

//	//manually set up some solver parameters
//	sm.macro_ndofs = sm.ndofs;
//	sm.dsol = (double*) tl_alloc(sizeof(double), sm.ndofs);
//	
//	sm.dirichlet_data = (double*) tl_alloc(sizeof(double), sm.ndofs);
//	sm.scale_vect = (double*) tl_alloc(sizeof(double), sm.ndofs);
////    for(int i=0;i<sm.ndofs;i++){
////    	sm.scale_vect[i] = 1.0;
////    }
//	sm.tol_nonlin = 1e-5;
//	sm.inc_nonlin = 1e-3;
//	sm.max_nonlin_linesearch_cuts = 5;
//	sm.it_count_nonlin_failed = 0;
//	sm.max_nonlin_it = 20;
//	sm.LINEAR_PROBLEM = NO;
//	sm.force_nonlin_it = NO;
//	sm.force_nonlin_it = NO;
//	sm.nonlinear_it_total = 0;
//	sm.nghost=0;
//	sm.ghosts = NULL;
//	sm.local_size = sm.ndofs;
//	sm.sol_old = (double*) tl_alloc(sizeof(double), sm.ndofs);
//	sm.sol_older = (double*) tl_alloc(sizeof(double), sm.ndofs);

 	SMODEL_SUPER *sm;
	sm = &(dm.superModel[0]);

	printf("SETTING UP BCMASK\n");

	// intialize dirichlet and old sol (initial guess)
	for (int local_index=0; local_index<dm.ndofs[0]; local_index++){

		dm.superModel[0].dirichlet_data[local_index] = 0.0;
		dm.superModel[0].sol_old[local_index] = 20.0;
		dm.superModel[0].sol[local_index] = 20.0;
		dm.superModel[0].lin_sys->dsol[local_index] = 0.0;
		dm.superModel[0].bc_mask[local_index] = YES;
	}

	//overwrite intial condition
	double x_coord, y_coord;
	int id;
	for (int i=0; i<nnodes; i++){
		//mark the boundary only
		x_coord = grid->node[i].x;
		y_coord = grid->node[i].y;
		//id = grid->node[i].id;
		id=i;
		//need to set IC
		dm.superModel[0].sol[id*3+1] = 1 + x_coord*x_coord + alpha * y_coord*y_coord;
		if ( is_near(x_coord,xmin) || is_near(x_coord,xmax) || is_near(y_coord,ymin) || is_near(y_coord,ymax) ){
			continue;
		}else{
			dm.superModel[0].bc_mask[id*3+1]=NO;
		}
		//printf("Dirichlet data node[%d] = %f\n", i, sm.dirichlet_data[i*3+1]);
	}
	printf("BCMASK COMPLETE\n");
	//set up bc mask

//	for (int i=0;i<sm.ndofs;i++){
//		printf("sm bc mask[%d] = %d\n",i,sm.bc_mask[i]);
//	}


	//see if it works
	//apply_Dirichlet_BC(&sm);
	//see if something is happening within newton loop or something else
//	initialize_system(&sm);
//	assemble_residual(&sm,sm.grid);
//	assemble_jacobian(&sm,sm.grid);
//	apply_Dirichlet_BC(&sm);
//	int status;
//	status = prep_umfpack(sm.indptr_diag,sm.cols_diag,sm.vals_diag, sm.dsol, sm.residual, sm.local_size);
//	status = solve_umfpack(sm.dsol, sm.indptr_diag, sm.cols_diag, sm.vals_diag, sm.residual, sm.local_size);
//	increment_function(&sm);
	//Screen_print_CSR(sm.indptr_diag, sm.cols_diag, sm.vals_diag, sm.ndofs);

	//set forward_step and call timeloop
	time_loop(&dm); 

	//compare with analytic solution
	//it is a scalar
	double *u_exact;
	u_exact = (double*) tl_alloc(sizeof(double),nnodes);
	compute_exact_solution_heat(u_exact, nnodes, grid,tf);
//	for(int i=0;i<nnodes;i++){
//		printf("Exact solution[%d] = %f\n",i,u_exact[i]);
//	}
	
	//extract second variable here
	double *uh;
	uh = (double*) tl_alloc(sizeof(double), nnodes);
	//create temporary integer array for nodes
	int *nodes;
	nodes = (int*) tl_alloc(sizeof(int), nnodes);
	for(int i=0;i<nnodes;i++){
		nodes[i] = i;
	}

	//global_to_local_dbl_cg_2(uh, sm.sol, nodes, nnodes, PERTURB_U, sm.node_physics_mat, sm.node_physics_mat_id);
	//global_to_local_dbl_cg(uh, sm->sol, nodes, nnodes, PERTURB_U, sm->dof_map_local, sm->node_physics_mat, sm->node_physics_mat_id);
	global_to_local_dbl_cg(uh, sm->sol, nodes, nnodes, PERTURB_U, sm->dof_map_local, sm->node_physics_mat);
//not needed anymore since nodes are reordered
//if (grid->inv_per_node!=NULL){
//	permute_array(uh,grid->inv_per_node,nnodes);
//}
//	printf("Final solution:\n");
//	for(int i=0; i<nnodes;i++){
//		printf("sol[%d] = %f, exact sol[%d] = %f\n",i,uh[i],i,u_exact[i]);
//	} 

	//compute L2 and Linf error
	double l2_err =  l2_error(uh, u_exact, nnodes);
	double linf_err =  linf_error(uh, u_exact, nnodes);

	printf("Final errors: %6.4e , %6.4e\n", l2_err,linf_err);

	//plot grid in h5?
//    strcpy(sm.grid->filename, "residtest");
//    init_hdf5_file(sm.grid);
//    printf("hdf5 initialized\n");
//    sgrid_write_hdf5(sm.grid);
//    printf("hdf5 written\n");
//    sgrid_write_xdmf(sm.grid);
//    printf("xmf written\n");

	//return -1 if failed, 0 if good
	
	if(l2_err < NEWTON_TEST_TOL && linf_err < NEWTON_TEST_TOL){
		err_code=0;
	}
	//printf("Final error code %d\n",err_code);

	//free memory
	u_exact = (double *) tl_free(sizeof(double), nnodes, u_exact);
	uh = (double *) tl_free(sizeof(double), nnodes, uh);
	nodes = (int *) tl_free(sizeof(int), nnodes, nodes);


	smodel_design_free(&dm);


	return err_code;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Computes analytic solution used to calculate error
 *  solution
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void compute_exact_solution_heat(double *u_exact, int ndof, SGRID *grid, double t){

	//test case comes from
	//https://jsdokken.com/dolfinx-tutorial/chapter2/heat_code.html

	//problem is du/dt - /\u + f = 0 on Omega
	//f = beta - 2 - 2*alpha
	// beta = 1.2
	// alpha = 3.0
	//u_D = u_exact on dOmega
	//u_exact = 1 + x^2 + alpha*y^2+beta*t

	//works for cg only at the moment, would need cell by cell loop for dg
	for(int i =0; i<ndof ; i++){
		u_exact[i] = 1.0 + grid->node[i].x*grid->node[i].x + alpha*grid->node[i].y*grid->node[i].y + beta*t;
	}



}
void permute_array(double *arr,int *p, int n){
	double *temp;
	temp = (double *) tl_alloc(sizeof(double),n);
	for(int i =0;i<n;i++){
		temp[p[i]] = arr[i];
	}
	// Copy permuted elements back to the original array
    for (int i = 0; i < n; i++) {
        arr[i] = temp[i];
    }
	temp = (double *) tl_free(sizeof(double),n,temp);
}
