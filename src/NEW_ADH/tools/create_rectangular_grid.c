#include "adh.h"


SGRID create_rectangular_grid(double xmin, double xmax, double ymin, double ymax, int npx, int npy,
 double theta, int dz, double a0, double ax, double ax2, double ay, double ay2, double axy,
 double ax2y, double axy2, double ax2y2, int flag3d ) {

	// MPI
//#ifdef _MESSG
//    debug_initialize(MPI_COMM_WORLD);
//#else
//    debug_initialize();
//#endif
    	
	//if (argc != 20 || strcmp(argv[1], "-help") == 0) { // Help option
//		printf("flumeBuilder [root_file_name] [xmin] [xmax] [ymin] [ymax] [npx] [npy] [theta] [dz] [a0] [ax] [ax2] [ay] [ay2] [axy] [ax2y] [axy2] [ax2y2] [flag_3D]\n");
//		printf("[root_file_name]: name for output\n");
//                printf("[xmin]: unrotated minimum x (length of flume)\n");
//                printf("[xmax]: unrotated maximum x (length of flume)\n");
//                printf("[ymin]: unrotated minimum y (width of flume)\n");
//                printf("[ymax]: unrotated maximum y (width of flume)\n");
//                printf("[npx]: number of nodes spanning length\n");
//                printf("[npy]: number of nodes spanning width\n");
//                printf("[theta]: rotation applied to flume nodes in xy plane in degrees\n");
//                printf("[dz]: distance by which to extrude the mesh if flag_3D is on, must be more than 0\n");
//                printf("[a0]: coefficient of degree 0\n");
//                printf("[ax]: coefficient\n");
//                printf("[ax2]: coefficient\n");
//                printf("[ay]: coefficient\n");
//                printf("[ay2]: coefficient\n");
//                printf("[axy]: coefficient\n");
//                printf("[ax2y]: coefficient\n");
//                printf("[axy2]: coefficient\n");
//                printf("[ax2y2]: coefficient\n");
//                printf("for a given node, z-value is\nh(x, y) = a0 + ax*x + ax2*x^2 + ay*y + ay2*y^2 + axy*x*y + ax2y*x^2*y + axy2*x*y^2 + ax2y2*x^2*y^2\n");
//                printf("[flag_3D]: whether to vertically extrude the mesh\n");
//		exit(0); // Not sure if exit(); is necessary or not
//	}


 	assert(npx > 1); // More than 1
	assert(npy > 1);
	// if (npx > 2) assert(npx % 2 != 0); // 2 or odd
	// if (npy > 2) assert(npy % 2 != 0);
	
	theta *= 3.141592653589793 / 180.; // Convert to radians
	
	assert(dz > 0);

	
	// Create upwards reference vector
	SVECT ref_vector;
	ref_vector.x = 0;
	ref_vector.y = 0;
	ref_vector.z = 1; // Should be positive


	// Start writing nodes to SGRID object
	double dx = (xmax - xmin)/(double)(npx-1);
	double dy = (ymax - ymin)/(double)(npy-1);


	//create a grid to return
	SGRID grid;

	grid.elem2d = NULL;
	grid.inv_per_node = NULL;

	//for now serial only, maybe add scotch in later


	//set some parameters, hard coded for triangles for now
	grid.nnodes = npx*npy;
	grid.max_nnodes = grid.nnodes;
	grid.nelems3d = 0;
	grid.nelems2d = (npx-1)*(npy-1)*2;
	grid.nelems1d = 0;
	grid.max_nelems1d = 0;
	grid.max_nelems2d = (npx-1)*(npy-1)*2;
	grid.max_nelems3d = 0;
	grid.nTets=0;
	grid.nPrisms=0;
	grid.nTris= grid.nelems2d;
	grid.nQuads=0;
	grid.nnodes_sur = grid.nnodes;      // number of ghost + residential surface nodes
    grid.nnodes_bed = 0;      // number of ghost + residential bed nodes

    grid.inv_per_node = (int*) tl_alloc(sizeof(int),grid.nnodes);

	// HPC totals, serial only for now
    grid.macro_nnodes = grid.nnodes;          // total # of nodes in the global mesh
    grid.macro_nnodes_sur = grid.nnodes;      // total # of surface nodes
    grid.macro_nnodes_bed = 0;      // total # of bed nodes
    grid.macro_nelems3d = grid.nelems3d;        // total # of 3d elements in the global mesh
    grid.macro_nelems2d = grid.nelems2d;        // total # of 2d elements in the global mesh
    grid.macro_nelems1d = grid.nelems1d;        // total # of 1d elements in the global mesh
    grid.macro_nelems2d_bed = 0;    // total # of 2d bed/surface elements across all PEs on the global mesh
    grid.macro_nTets = grid.nTets;           // total # of tetrahedron elements in global mesh
    grid.macro_nPrisms = grid.nPrisms;         // total # of prism elements in global mesh
    grid.macro_nQuads = grid.nQuads;          // total # of quadrilateral elements in global mesh
    grid.macro_nTris = grid.nTris;           // total # of triangle elements in global mesh
    grid.orig_macro_nnodes = grid.nnodes;     // total number of original nodes in the global mesh
    grid.orig_macro_nnodes_sur = grid.nnodes; // total number of original surface nodes in the global mesh
    grid.orig_macro_nnodes_bed = 0; // total number of original bed nodes in the global mesh
    grid.orig_macro_nelems1d = grid.macro_nelems1d;   // total number of original elements in the global mesh
    grid.orig_macro_nelems2d = grid.macro_nelems2d;   // total number of original elements in the global mesh
    grid.orig_macro_nelems3d = grid.macro_nelems3d;   // total number of original elements in the global mesh
    grid.my_nnodes = grid.nnodes;          // total # of residential only nodes in the global mesh
    grid.my_nnodes_sur = grid.nnodes_sur;      // total # of residential only surface nodes in the global mesh
    grid.my_nnodes_bed = grid.nnodes_bed;      // total # of residential only bed nodes in the global mesh
    grid.my_nelems3d = grid.nelems3d;        // total # of residential only 3d elements in the global mesh
    grid.my_nelems2d = grid.nelems2d;        // total # of residential only 2d elements in the global mesh
    grid.my_nelems1d = grid.nelems1d;        // total # of residential only 1d elements in the global mesh
    grid.my_nTets = grid.nTets;           // total # of residential only tetrahedron elements in global mesh
    grid.my_nPrisms = grid.nPrisms;         // total # of residential only prism elements in global mesh
    grid.my_nQuads = grid.nQuads;          // total # of residential only quadrilateral elements in global mesh
    grid.my_nTris = grid.nTris;           // total # of residential only triangle elements in global mesh

	if (grid.nelems1d > 0) {selem1d_init_alloc_array(&(grid.elem1d), grid.nelems1d);}
    if (grid.nelems2d > 0) {selem2d_init_alloc_array(&(grid.elem2d), grid.nelems2d);}
    if (grid.nelems3d > 0) {selem3d_init_alloc_array(&(grid.elem3d), grid.nelems3d);}
    
    snode_init_alloc_array(&(grid.node), grid.nnodes);


    //hard code for triangles now
    int nnodes_on_elem = 3;
    //have a temp vector of nodes
    int temp_ids[nnodes_on_elem];
    SVECT temp_node[nnodes_on_elem];
    svect_init_array(temp_node, nnodes_on_elem);

	//printf("len is %d\n", (npx*npy));

	int i, j, k=1;
	double x = xmin, xr, zr;
	double y = ymin, yr;
	int npoints = 0;
	for (i=0; i<npx; i++) {
		for (j=0; j<npy; j++) {
			x = xmin + i*dx;
			y = ymin + j*dy; // Originanlly ymax - j*dy;...not sure why
			
			// Apply rotation to xy coordinates
			xr = x*cos(theta) - y*sin(theta);
			yr = x*sin(theta) + y*cos(theta);
			zr = a0 + ax*x + ax2*x*x + ay*y + ay2*y*y + axy*x*y + ax2y*x*x*y + axy2*x*y*y + ax2y2*x*x*y*y;
			//fprintf(fp,"ND %d %20.10f %20.10f %20.10f\n", k, xr, yr, zr);
			grid.node[npoints].id = npoints;
			grid.node[npoints].original_id = npoints;
			//serial for now
			grid.node[npoints].gid = grid.node[npoints].id;
			grid.node[npoints].x = xr;
			grid.node[npoints].y = yr;
			grid.node[npoints].z = zr;
			
			k++;
			npoints++;
		}
	}

	// Define counter-clockwise elements
	// Use a cross-product to verify
	// For counter-clockwise node ordering ABC,
	// the cross product AB x BC will point up
	SVECT cross;
	SVECT side1, side2;

	int inode = 0, nd1, nd2, nd3, temp;
	k = 0; // Element counter
	for (i = 0; i < npx - 1; i++) { // Exclude last column
		for (j = 0; j < npy - 1; j++) { // Exclude last row
			inode = i * npy + j; // i, j mapping to vector index
					     // Bottom-left corner
			//
			//
			// "forward slash" square [/]
			//
			//
			if (j % 2 == 0) {
				//
				// Lower triangle '/]'
				//

				nd1 = inode;
				nd2 = inode + npy;
				nd3 = inode + npy + 1;


				
				side1.x = grid.node[nd2].x - grid.node[nd1].x;
               	side1.y = grid.node[nd2].y - grid.node[nd1].y;
               	side1.z = grid.node[nd2].z - grid.node[nd1].z;
               	side2.x = grid.node[nd3].x - grid.node[nd1].x;
        	    side2.y = grid.node[nd3].y - grid.node[nd1].y;
	            side2.z = grid.node[nd3].z - grid.node[nd1].z;

				cross = svect_cross(side1, side2); // Should be positive
				// printf("Cross=%.2f\n", cross); // ******************

	            if (svect_dotp(cross,ref_vector) <= 0) { // Vectors do not align
					temp = nd3;
					nd3 = nd2;
					nd2 = temp;

                    side1.x = grid.node[nd2].x - grid.node[nd1].x;
                    side1.y = grid.node[nd2].y - grid.node[nd1].y;
                    side1.z = grid.node[nd2].z - grid.node[nd1].z;
                    side2.x = grid.node[nd3].x - grid.node[nd1].x;
                    side2.y = grid.node[nd3].y - grid.node[nd1].y;
                    side2.z = grid.node[nd3].z - grid.node[nd1].z;

                    cross = svect_cross(side1, side2);
                    //printf("Revised cross=%.2f\n", cross); // ******************
                    assert(svect_dotp(cross,ref_vector) >= 0);
				}

				// Note usage of one-index
				// Console
				//printf("nodes: %d %d %d \t cross: %20.10f %20.10f %20.10f \n", \
				//		nd1+1,nd2+1,nd3+1,cross.x,cross.y,cross.z);
            			// File
				//allocate and save to element
				temp_node[0].x = grid.node[nd1].x;
				temp_node[0].y = grid.node[nd1].y;
				temp_node[0].z = grid.node[nd1].z;
				temp_node[1].x = grid.node[nd2].x;
				temp_node[1].y = grid.node[nd2].y;
				temp_node[1].z = grid.node[nd2].z;
				temp_node[2].x = grid.node[nd3].x;
				temp_node[2].y = grid.node[nd3].y;
				temp_node[2].z = grid.node[nd3].z;
				temp_ids[0] = nd1;
				temp_ids[1] = nd2;
				temp_ids[2] = nd3;
				selem2d_init(&grid.elem2d[k]);
        		selem2d_load(&(grid.elem2d[k]),k,k,nnodes_on_elem,temp_ids,0,temp_node,0);
				grid.elem2d[k].nodes[0] = nd1;
				grid.elem2d[k].nodes[1] = nd2;
				grid.elem2d[k].nodes[2] = nd3;
				//printf("E3T %d %d %d %d 1\n",k+1,nd1+1,nd2+1,nd3+1);
            	k++;

				//
				// Upper triangle '[/'
				//

				nd1 = inode;
				nd2 = inode + npy + 1;
				nd3 = inode + 1;
				
				side1.x = grid.node[nd2].x - grid. node[nd1].x;
                side1.y = grid.node[nd2].y - grid.node[nd1].y;
                side1.z = grid.node[nd2].z - grid.node[nd1].z;
                side2.x = grid.node[nd3].x - grid.node[nd1].x;
                side2.y = grid.node[nd3].y - grid.node[nd1].y;
                side2.z = grid.node[nd3].z - grid.node[nd1].z;

                cross = svect_cross(side1, side2); // Should be positive
                // printf("Cross=%.2f\n", cross); // ******************

                if (svect_dotp(cross,ref_vector) <= 0) { // Vectors do not align
                    temp = nd3;
                    nd3 = nd2;
                    nd2 = temp;
					side1.x = grid.node[nd2].x - grid.node[nd1].x;
					side1.y = grid.node[nd2].y - grid.node[nd1].y;
					side1.z = grid.node[nd2].z - grid.node[nd1].z;
					side2.x = grid.node[nd3].x - grid.node[nd1].x;
					side2.y = grid.node[nd3].y - grid.node[nd1].y;
					side2.z = grid.node[nd3].z - grid.node[nd1].z;

					cross = svect_cross(side1, side2);
                    //printf("Revised cross=%.2f\n", cross); // ******************
                    assert(svect_dotp(cross,ref_vector) >= 0);
                }

                // Note usage of one-index
                // Console
                //printf("nodes: %d %d %d \t cross: %20.10f %20.10f %20.10f \n", \
                //nd1+1,nd2+1,nd3+1,cross.x,cross.y,cross.z);
                // save to grid
                //allocate and save to element
				temp_node[0].x = grid.node[nd1].x;
				temp_node[0].y = grid.node[nd1].y;
				temp_node[0].z = grid.node[nd1].z;
				temp_node[1].x = grid.node[nd2].x;
				temp_node[1].y = grid.node[nd2].y;
				temp_node[1].z = grid.node[nd2].z;
				temp_node[2].x = grid.node[nd3].x;
				temp_node[2].y = grid.node[nd3].y;
				temp_node[2].z = grid.node[nd3].z;
				temp_ids[0] = nd1;
				temp_ids[1] = nd2;
				temp_ids[2] = nd3;
				selem2d_init(&grid.elem2d[k]);
        		selem2d_load(&(grid.elem2d[k]),k,k,nnodes_on_elem,temp_ids,0,temp_node,0);
                grid.elem2d[k].nodes[0] = nd1;
				grid.elem2d[k].nodes[1] = nd2;
				grid.elem2d[k].nodes[2] = nd3;
                //printf("E3T %d %d %d %d 1\n",k+1,nd1+1,nd2+1,nd3+1);
                k++;

			}
			//
			//
			// "backward slash" square [\]
			//
			//
			else {
				//
				// Lower triangle '[\'
				//
				
				nd1 = inode;
				nd2 = inode + npy;
				nd3 = inode + 1;
				
				side1.x = grid.node[nd2].x - grid.node[nd1].x;
                side1.y = grid.node[nd2].y - grid.node[nd1].y;
                side1.z = grid.node[nd2].z - grid.node[nd1].z;
                side2.x = grid.node[nd3].x - grid.node[nd1].x;
                side2.y = grid.node[nd3].y - grid.node[nd1].y;
                side2.z = grid.node[nd3].z - grid.node[nd1].z;

                cross = svect_cross(side1, side2); // Should be positive
                // printf("Cross=%.2f\n", cross); // ******************

                if (svect_dotp(cross,ref_vector) <= 0) { // Vectors do not align
                	temp = nd3;
                    nd3 = nd2;
                    nd2 = temp;

                    side1.x = grid.node[nd2].x - grid.node[nd1].x;
                    side1.y = grid.node[nd2].y - grid.node[nd1].y;
                    side1.z = grid.node[nd2].z - grid.node[nd1].z;
                    side2.x = grid.node[nd3].x - grid.node[nd1].x;
                    side2.y = grid.node[nd3].y - grid.node[nd1].y;
                    side2.z = grid.node[nd3].z - grid.node[nd1].z;

                    cross = svect_cross(side1, side2);
                    // printf("Revised cross=%.2f\n", cross); // ******************
                    assert(svect_dotp(cross,ref_vector) >= 0);
                }

                // Note usage of one-index
                // Console
                //printf("nodes: %d %d %d \t cross: %20.10f %20.10f %20.10f \n", \
                //nd1+1,nd2+1,nd3+1,cross.x,cross.y,cross.z);
                // File
                //save to grid
                //allocate and save to element
				temp_node[0].x = grid.node[nd1].x;
				temp_node[0].y = grid.node[nd1].y;
				temp_node[0].z = grid.node[nd1].z;
				temp_node[1].x = grid.node[nd2].x;
				temp_node[1].y = grid.node[nd2].y;
				temp_node[1].z = grid.node[nd2].z;
				temp_node[2].x = grid.node[nd3].x;
				temp_node[2].y = grid.node[nd3].y;
				temp_node[2].z = grid.node[nd3].z;

				temp_ids[0] = nd1;
				temp_ids[1] = nd2;
				temp_ids[2] = nd3;
				selem2d_init(&grid.elem2d[k]);
        		selem2d_load(&(grid.elem2d[k]),k,k,nnodes_on_elem,temp_ids,0,temp_node,0);
                grid.elem2d[k].nodes[0] = nd1;
				grid.elem2d[k].nodes[1] = nd2;
				grid.elem2d[k].nodes[2] = nd3;

                //printf("E3T %d %d %d %d 1\n",k+1,nd1+1,nd2+1,nd3+1);
                k++;

				//
				// Upper triangle '\]'
				//
				
				nd1 = inode + npy;
				nd2 = inode + npy + 1;
				nd3 = inode + 1;
				
				side1.x = grid.node[nd2].x - grid.node[nd1].x;
                side1.y = grid.node[nd2].y - grid.node[nd1].y;
                side1.z = grid.node[nd2].z - grid.node[nd1].z;
                side2.x = grid.node[nd3].x - grid.node[nd1].x;
                side2.y = grid.node[nd3].y - grid.node[nd1].y;
                side2.z = grid.node[nd3].z - grid.node[nd1].z;

                cross = svect_cross(side1, side2); // Should be positive
                // printf("Cross=%.2f\n", cross); // ******************

                if (svect_dotp(cross,ref_vector) <= 0) { // Vectors do not align
                    temp = nd3;
                    nd3 = nd2;
                    nd2 = temp;

					side1.x = grid.node[nd2].x - grid.node[nd1].x;
                    side1.y = grid.node[nd2].y - grid.node[nd1].y;
                    side1.z = grid.node[nd2].z - grid.node[nd1].z;
                	side2.x = grid.node[nd3].x - grid.node[nd1].x;
        	        side2.y = grid.node[nd3].y - grid.node[nd1].y;
	                side2.z = grid.node[nd3].z - grid.node[nd1].z;

					cross = svect_cross(side1, side2);
					// printf("Revised cross=%.2f\n", cross); // ******************
					assert(svect_dotp(cross,ref_vector) >= 0);
                }

                // Note usage of one-index
                // Console
				//printf("nodes: %d %d %d \t cross: %20.10f %20.10f %20.10f \n", \
				//		nd1+1,nd2+1,nd3+1,cross.x,cross.y,cross.z);
                // save to frid object
                //allocate and save to element
				temp_node[0].x = grid.node[nd1].x;
				temp_node[0].y = grid.node[nd1].y;
				temp_node[0].z = grid.node[nd1].z;
				temp_node[1].x = grid.node[nd2].x;
				temp_node[1].y = grid.node[nd2].y;
				temp_node[1].z = grid.node[nd2].z;
				temp_node[2].x = grid.node[nd3].x;
				temp_node[2].y = grid.node[nd3].y;
				temp_node[2].z = grid.node[nd3].z;
				temp_ids[0] = nd1;
				temp_ids[1] = nd2;
				temp_ids[2] = nd3;
				selem2d_init(&grid.elem2d[k]);
        		selem2d_load(&(grid.elem2d[k]),k,k,nnodes_on_elem,temp_ids,0,temp_node,0);
                grid.elem2d[k].nodes[0] = nd1;
				grid.elem2d[k].nodes[1] = nd2;
				grid.elem2d[k].nodes[2] = nd3;
				//printf("E3T %d %d %d %d 1\n",k+1,nd1+1,nd2+1,nd3+1);
				k++;
			}
		}
	}
	
	return grid;
}
