/*! \file  jacobian_test.c This file tests the assembly of a Jacobian matrix */
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     This function tests assembly of the Jacobian matrix
 *  \author    Count Corey J. Trahan
 *  \author    Mark Loveland
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int jacobian_test(int argc, char **argv) {
	//create a grid
	SGRID grid;
	//let's just do something simple
	//3x3 triangular element grid
	double xmin = 0.0;
	double xmax = 2.0;
	double ymin = 0.0;
	double ymax = 2.0;
	int npx = 3;
	int npy = 3;
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
    grid = create_rectangular_grid(xmin, xmax, ymin, ymax, npx, npy,
 	theta, dz, a0, ax, ax2, ay, ay2, axy,
    ax2y, axy2, ax2y2, flag3d );
	//print coordinates
  for(int local_index =0; local_index<grid.nnodes; local_index++){
		printf("Node %d: (x,y) = {%f,%f}\n",grid.node[local_index].gid,grid.node[local_index].x,grid.node[local_index].y);
	}
//	//print connectivity
//	for(int local_index =0; local_index<grid.nelems2d; local_index++){
//		printf("Element %d: (nd1,nd2,nd3) = {%d,%d,%d}\n",local_index ,grid.elem2d[local_index].nodes[0], grid.elem2d[local_index].nodes[1], grid.elem2d[local_index].nodes[2]);
//	}
	//create a supermodel with a grid in it
	SMODEL_SUPER sm;
	sm.grid = &grid;
	//specify elemental physics and other properties in super model
	double dt = 1.0;
	double t0 = 0.0;
	double tf = 1.0;
	char elemVarCode[4]; 
	strcpy(&elemVarCode[0],"2");//SW2D
	strcpy(&elemVarCode[1],"0"); //GW
	strcpy(&elemVarCode[2],"0"); //Transport
	smodel_super_no_read_simple(&sm, dt, t0, tf, 0 , 1, 0, elemVarCode);
	printf("Supermodel no read complete\n");
	//allocate linear system
	//doesn't currently work, need to go back and fix
	//fe_allocate_initialize_linear_system(&sm);
	sm.cols_diag = NULL;
	sm.vals_diag = NULL;
	create_sparsity_split_CSR(&sm, sm.grid);
	//do we want to stor nnz? it is stored in sm->indptr[nrows]
    //do we want to store local_size = local_range[1]-local_range[0]
    sarray_init_dbl(sm.vals_diag, sm.indptr_diag[sm.my_ndofs-1]);
    //sarray_init_dbl(sm.vals_off_diag, sm.indptr_off_diag[sm.my_ndofs-1]);
	//assemble a residual and check correctness
	assemble_residual(&sm, sm.grid);
	//print final residual
	for(int local_index=0;local_index<grid.nnodes;local_index++){
		printf("Node %d: (x,y) = {%f,%f}, Residual = {%f,%f,%f}\n",grid.node[local_index].gid,grid.node[local_index].x,grid.node[local_index].y,sm.residual[local_index*3],sm.residual[local_index*3+1],sm.residual[local_index*3+2]);
	}
	//assemble a jacobian and check correcntess
	assemble_jacobian(&sm);
	Screen_print_CSR(sm.indptr_diag, sm.cols_diag, sm.vals_diag, sm.ndofs);
	//plot grid in h5?
    strcpy(sm.grid->filename, "residtest");
    init_hdf5_file(sm.grid);
    printf("hdf5 initialized\n");
    sgrid_write_hdf5(sm.grid);
    printf("hdf5 written\n");
    sgrid_write_xdmf(sm.grid);
    printf("xmf written\n");
	return 0;
}