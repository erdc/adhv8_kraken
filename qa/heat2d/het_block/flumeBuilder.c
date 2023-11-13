#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

typedef struct {
    double x, y, z;               /* coordinates of the vector */
} SVECT;

SVECT svect_cross(SVECT v1, SVECT v2) {
    SVECT result;
    result.x = v1.y*v2.z - v1.z*v2.y;
    result.y = v1.z*v2.x - v1.x*v2.z;
    result.z = v1.x*v2.y - v1.y*v2.x;
    return result;
}

double svect_dotp(SVECT v1, SVECT v2) {
    return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

int main() {

    SVECT ref_vector;
    ref_vector.x = 0 - 0;
    ref_vector.y = 0 - 0;
    ref_vector.z = -1 - 0;

    double xmin = 0;
    double xmax = 50;
    double ymin = -50;
    double ymax = 0;
    double zinit = 1.0; /*!SD added User Defined initial depth of water surface (used in hotstart file)*/
    double cinit_in = 1.0; /*!SD added User Defined initial concentration near left boundary (used in hotstart file)*/
    double cinit_out = 0.2; /*!SD added User Defined initial concentration near right boundary (used in hotstart file)*/
    double c_in = 0.0; /*!SD added User Defined no-flow side walls boundary (used in hotstart file)*/
    double cDir_base = 0.0; /*!SD added Homogeneous Dirichlet boundary data*/
    double cDir_top = 1.0; /*!SD added Homogeneous Dirichlet boundary data*/
    double diff_mat1 = 1.0e-2; /*diffusion in material 1 (center)*/
    double diff_mat2 = 1.0;   /*diffusion in material 2 */
    /*!SD: change the values below for different offline runs*/
    int npx = 101; //51
    int npy = 101; //11
    double DT = 0.5; // For computation
    double DT_save = 2.0; // For writing solution data
    double T_final = 100.0; // Final time point for simulation

   // assert(npy % 2 != 0);
    double zb = -10.0;

    /* change the angle below according to orientation of domain */
    double theta = 0.0;//45 * 3.141592653589793 / 180.; //0.0

    /*!SD: velocity components in the rotated axes for a constant horizontal vel U = 1*/
    double xvel = 0.0;//0.070710678;
    double yvel = 0.0;//0.070710678;


    FILE *fp;
    fp=fopen("heat.3dm", "w");  // !SD changed from 2dm to 3dm

    fprintf(fp,"MESH2D\n");

    double dx = (xmax - xmin)/(double)(npx-1);
    double dy = (ymax - ymin)/(double)(npy-1);

    SVECT node[npx*npy];

    int i, j, k=1;
    double x = xmin, xr, zr;
    double y = ymin, yr;
    double H0 = 0.0; //1.e-9;  // !SD set this to zero to match ROM flume domain
    int npoints = 0;
    for (i=0; i<npx; i++) {
        for (j=0; j<npy; j++) {
            x = xmin + i*dx;
            y = ymax - j*dy; //ymin + j*dy;

            // rotate 45 degrees
            xr = x*cos(theta) - y*sin(theta);
            yr = x*sin(theta) + y*cos(theta);
            zr = -H0 * xr * xr;
            fprintf(fp,"ND %d %20.10f %20.10f %20.10f\n",k,xr,yr,zr);

            node[npoints].x = xr;
            node[npoints].y = yr;
            node[npoints].z = -H0 * xr * xr;

            k++;
            npoints++;
        }
    }

    SVECT cross;
    SVECT side1, side2;
    int igrp = 0, inode=0, nd1, nd2, nd3;
    k=0;

    for (i=0; i<npx-1; i++) {
        for (igrp=0; igrp<npy-1; igrp++) {
            inode = i*npy + igrp;

            // make sure node numbering is counter clock-wise
            nd1 = inode; nd2 = inode+1; nd3 = inode+npy+1;
            side1.x =  node[nd2].x -  node[nd1].x; side1.y =  node[nd2].y -  node[nd1].y; side1.z =  node[nd2].z -  node[nd1].z;
            side2.x =  node[nd3].x -  node[nd1].x; side2.y =  node[nd3].y -  node[nd1].y; side2.z =  node[nd3].z -  node[nd1].z;
            cross = svect_cross(side1,side2);

            if (svect_dotp(cross,ref_vector) > 0) {
                //if (cross.z > 0) {
                nd1 = inode; nd2 = inode+npy+1; nd3 = inode+1;
                side1.x =  node[nd2].x -  node[nd1].x; side1.y =  node[nd2].y -  node[nd1].y; side1.z =  node[nd2].z -  node[nd1].z;
                side2.x =  node[nd3].x -  node[nd1].x; side2.y =  node[nd3].y -  node[nd1].y; side2.z =  node[nd3].z -  node[nd1].z;
                cross = svect_cross(side1,side2);
                assert(svect_dotp(cross,ref_vector) < 0);
                //assert(cross.z < 0);
            }
            printf("nodes: %d %d %d \t cross: %20.10f %20.10f %20.10f \n",nd1+1,nd2+1,nd3+1,cross.x,cross.y,cross.z);
	    
	    if (i>= npx/4 && i <= 3*npx/4 && igrp >= npy/4 && igrp <=3*npy/4) {
            // original output
            	fprintf(fp,"E3T %d %d %d %d 1\n",k+1,nd1+1,nd2+1,nd3+1); }
	    else {
		fprintf(fp,"E3T %d %d %d %d 2\n",k+1,nd1+1,nd2+1,nd3+1); }
            k++;


            // make sure node numbering is counter clock-wise
            nd1 = inode; nd2 = inode+npy+1; nd3 = inode+npy;
            side1.x =  node[nd2].x -  node[nd1].x; side1.y =  node[nd2].y -  node[nd1].y; side1.z =  node[nd2].z -  node[nd1].z;
            side2.x =  node[nd3].x -  node[nd1].x; side2.y =  node[nd3].y -  node[nd1].y; side2.z =  node[nd3].z -  node[nd1].z;
            cross = svect_cross(side1,side2);
            if (svect_dotp(cross,ref_vector) > 0) {
                //if (cross.z > 0) {
                nd1 = inode; nd2 = inode+npy; nd3 = inode+npy+1;
                side1.x =  node[nd2].x -  node[nd1].x; side1.y =  node[nd2].y -  node[nd1].y; side1.z =  node[nd2].z -  node[nd1].z;
                side2.x =  node[nd3].x -  node[nd1].x; side2.y =  node[nd3].y -  node[nd1].y; side2.z =  node[nd3].z -  node[nd1].z;
                cross = svect_cross(side1,side2);
                //assert(cross.z < 0);
                assert(svect_dotp(cross,ref_vector) < 0);
            }

            printf("nodes: %d %d %d \t cross: %20.10f %20.10f %20.10f \n",nd1+1,nd2+1,nd3+1,cross.x,cross.y,cross.z);
            if (i>= npx/4 && i <= 3*npx/4 && igrp >= npy/4 && igrp <=3*npy/4) {
            // original output
            	fprintf(fp,"E3T %d %d %d %d 1\n",k+1,nd1+1,nd2+1,nd3+1); }
	    else {
		fprintf(fp,"E3T %d %d %d %d 2\n",k+1,nd1+1,nd2+1,nd3+1); }
            k++;

        }
    }


    int nelemsTotal = k;
    fclose(fp);

    fp=fopen("heat_adh.3dm", "w");  // !SD changed from 2dm to 3dm
    fprintf(fp,"MESH2D\n");

    SVECT cross_new;
    SVECT side1_new, side2_new;

    SVECT node_new[npx*npy];

    k=0;
    igrp = 0, inode=0;
    for (i=0; i<npx-1; i++) {
        for (igrp=0; igrp<npy-1; igrp++) {
            inode = i*npy + igrp;

            // make sure node numbering is counter clock-wise
            nd1 = inode; nd2 = inode+1; nd3 = inode+npy+1;
            side1_new.x =  node_new[nd2].x -  node_new[nd1].x; side1_new.y =  node_new[nd2].y -  node_new[nd1].y; side1_new.z =  node_new[nd2].z -  node_new[nd1].z;
            side2_new.x =  node_new[nd3].x -  node_new[nd1].x; side2_new.y =  node_new[nd3].y -  node_new[nd1].y; side2_new.z =  node_new[nd3].z -  node_new[nd1].z;
            cross_new = svect_cross(side1_new,side2_new);

            if (svect_dotp(cross_new,ref_vector) > 0) {
                //if (cross.z > 0) {
                nd1 = inode; nd2 = inode+npy+1; nd3 = inode+1;
                side1_new.x =  node_new[nd2].x -  node_new[nd1].x; side1_new.y =  node_new[nd2].y -  node_new[nd1].y; side1_new.z =  node_new[nd2].z -  node_new[nd1].z;
                side2_new.x =  node_new[nd3].x -  node_new[nd1].x; side2_new.y =  node_new[nd3].y -  node_new[nd1].y; side2_new.z =  node_new[nd3].z -  node_new[nd1].z;
                cross_new = svect_cross(side1_new,side2_new);
                assert(svect_dotp(cross_new,ref_vector) < 0);
                //assert(cross.z < 0);
            }
            //printf("nodes: %d %d %d \t cross: %20.10f %20.10f %20.10f \n",nd1+1,nd2+1,nd3+1,cross_new.x,cross_new.y,cross_new.z);

            // original output
            if (i>= npx/4 && i <= 3*npx/4 && igrp >= npy/4 && igrp <=3*npy/4) {
            // original output
            	fprintf(fp,"E3T %d %d %d %d 1\n",k+1,nd1+1,nd2+1,nd3+1); }
	    else {
		fprintf(fp,"E3T %d %d %d %d 2\n",k+1,nd1+1,nd2+1,nd3+1); }
            k++;


            // make sure node numbering is counter clock-wise
            nd1 = inode; nd2 = inode+npy+1; nd3 = inode+npy;
            side1_new.x =  node_new[nd2].x -  node_new[nd1].x; side1_new.y =  node_new[nd2].y -  node_new[nd1].y; side1_new.z =  node_new[nd2].z -  node_new[nd1].z;
            side2_new.x =  node_new[nd3].x -  node_new[nd1].x; side2_new.y =  node_new[nd3].y -  node_new[nd1].y; side2_new.z =  node_new[nd3].z -  node_new[nd1].z;
            cross_new = svect_cross(side1_new,side2_new);
            if (svect_dotp(cross_new,ref_vector) > 0) {
                //if (cross.z > 0) {
                nd1 = inode; nd2 = inode+npy; nd3 = inode+npy+1;
                side1_new.x =  node_new[nd2].x -  node_new[nd1].x; side1_new.y =  node_new[nd2].y -  node_new[nd1].y; side1_new.z =  node_new[nd2].z -  node_new[nd1].z;
                side2_new.x =  node_new[nd3].x -  node_new[nd1].x; side2_new.y =  node_new[nd3].y -  node_new[nd1].y; side2_new.z =  node_new[nd3].z -  node_new[nd1].z;
                cross_new = svect_cross(side1_new,side2_new);
                //assert(cross.z < 0);
                assert(svect_dotp(cross_new,ref_vector) < 0);
            }

            //printf("nodes: %d %d %d \t cross: %20.10f %20.10f %20.10f \n",nd1+1,nd2+1,nd3+1,cross_new.x,cross_new.y,cross_new.z);
            if (i>= npx/4 && i <= 3*npx/4 && igrp >= npy/4 && igrp <=3*npy/4) {
            // original output
            	fprintf(fp,"E3T %d %d %d %d 1\n",k+1,nd1+1,nd2+1,nd3+1); }
	    else {
		fprintf(fp,"E3T %d %d %d %d 2\n",k+1,nd1+1,nd2+1,nd3+1); }
            k++;
        }
    }

    k=1;
    x = xmin;
    y = ymin;
    H0 = 0.0; //1.e-9;  // !SD set this to zero to match ROM flume domain
    npoints = 0;
    for (i=0; i<npx; i++) {
        for (j=0; j<npy; j++) {
            x = xmin + i*dx;
            y = ymax - j*dy; //ymin + j*dy;

            // rotate 45 degrees
            xr = x*cos(theta) - y*sin(theta);
            yr = x*sin(theta) + y*cos(theta);
            zr = -H0 * xr * xr;
            fprintf(fp,"ND %d %20.10f %20.10f %20.10f\n",k,xr,yr,zr);

            node_new[npoints].x = xr;
            node_new[npoints].y = yr;
            node_new[npoints].z = -H0 * xr * xr;

            k++;
            npoints++;
        }
    }

    fclose(fp);
    
 
    fp=fopen("heat.bc", "w");
    fprintf(fp,"# operational parameters\n");
    fprintf(fp,"OP SW2\n");
    fprintf(fp,"OP TRN 1\n");
    fprintf(fp,"OP BLK 1\n");
    fprintf(fp,"OP INC 40\n");
    fprintf(fp,"OP PRE 1\n");
    fprintf(fp,"\n");
    fprintf(fp,"OP TEM 1\n");
    fprintf(fp,"\n");
    fprintf(fp,"OP TPG 0\n");
    fprintf(fp,"\n");
    fprintf(fp,"// Hydro is being skipped\n");
    fprintf(fp,"// NOTERM HYDRO 'X Vel' 'Y Vel' 'Z Vel' 'Displacement'\n");
    fprintf(fp,"NOTERM HYDRO %11.9f %11.9f 0.0 1.0\n",xvel,yvel);
    fprintf(fp,"\n");

    fprintf(fp,"# solver parameters\n");
    fprintf(fp,"IP NTL 1e-9\n");
    fprintf(fp,"IP ITL 1e-9\n");
    fprintf(fp,"IP NIT 10\n");
    fprintf(fp,"IP MIT 100\n");
    fprintf(fp,"\n");
    fprintf(fp,"CN CON 1 1 0 0 0\n");
    fprintf(fp,"SOUT RESID\n");
    fprintf(fp,"\n");

    fprintf(fp,"# general material properties\n");
    fprintf(fp,"MP ML 1 0\n");
    fprintf(fp,"MP TRT 1 1 1\n");
    fprintf(fp,"MP MUC 1.0\n");
    fprintf(fp,"MP MU 0.000000\n");
    fprintf(fp,"MP RHO 1000\n");
    fprintf(fp,"MP G 9.8\n");
    fprintf(fp,"\n");

    fprintf(fp,"# Material 1 properties\n");
    fprintf(fp,"MP DF 1 1 %11.9f\n",diff_mat1);  // Turbulent diffusion = diffusion coefficient in transport equation
    fprintf(fp,"MP EVS 1 0.0 0.0 0.0\n"); //Eddy viscosity = diffusion coefficient in sw2 equation
    fprintf(fp,"# 2d element material string\n");
    fprintf(fp,"MTS 1 1\n");
    fprintf(fp,"MP SRT 1 100\n");

    fprintf(fp,"# Material 2 properties\n");
    fprintf(fp,"MP DF 2 1 %11.9f\n",diff_mat2); // Turbulent diffusion = diffusion coefficient in transport equation
    fprintf(fp,"MP EVS 2 0.0 0.0 0.0\n"); //Eddy viscosity = diffusion coefficient in sw2 equation
    fprintf(fp,"# 2d element material string\n");
    fprintf(fp,"MTS 2 7\n");
    fprintf(fp,"MP SRT 2 100\n");

  
    fprintf(fp,"# friction coefficients\n");
    fprintf(fp,"FR MNG 1 0.0\n");
    fprintf(fp,"\n");

    fprintf(fp,"#TIME-SERIES --------------------\n\n");
    fprintf(fp,"# Tail water\n");
    fprintf(fp,"SERIES BC 1 2 0 0 0 0\n");
    fprintf(fp,"%-20.2f %-20.2f\n",0.0, 1.0);
    fprintf(fp,"%-20.2f %-20.2f\n",T_final, 1.0);
    fprintf(fp,"\n");

    fprintf(fp,"# Normal Natural Velocity\n");
    fprintf(fp,"SERIES BC 2 2 0 0 0 0\n");
    fprintf(fp,"%-20.2f %-20.9f\n",0.0, xvel); //0.070710678); //xvel);
    fprintf(fp,"%-20.2f %-20.9f\n",T_final, xvel); //0.070710678); // xvel);
    fprintf(fp,"\n");

    fprintf(fp,"# Output time-series \n");
    fprintf(fp,"SERIES AWRITE 3 1 0 0 0 0\n");
    fprintf(fp,"%-20.2f %-20.2f %-20.2f 0\n",0.0,T_final,DT_save);
    fprintf(fp,"\n");

    fprintf(fp,"# Time-step time-series\n");
    fprintf(fp,"SERIES DT 4 2 0 0 0 0\n");
    fprintf(fp,"%-20.2f %-20.2f\n",0.0,DT);
    fprintf(fp,"%-20.2f %-20.2f\n",T_final,DT);
    fprintf(fp,"\n");

    fprintf(fp,"# Transport Neumann time-series\n");
    fprintf(fp,"SERIES BC 5 2 0 0\n");
    fprintf(fp,"%-20.2f %-20.2f\n",0.0, c_in);
    fprintf(fp,"%-20.2f %-20.2f\n",T_final, c_in);
    fprintf(fp,"\n");
    fprintf(fp,"SERIES BC 6 2 0 0\n");
    fprintf(fp,"%-20.2f %-20.2f\n",0.0, c_in);
    fprintf(fp,"%-20.2f %-20.2f\n",T_final, c_in);
    fprintf(fp,"\n");

    fprintf(fp,"# Transport Dirichlet time-series\n");
    fprintf(fp,"SERIES BC 20 2 0 0\n");
    fprintf(fp,"%-20.2f %-20.2f\n",0.0, cDir_base);
    fprintf(fp,"%-20.2f %-20.2f\n",T_final, cDir_base);
    fprintf(fp,"\n");
    fprintf(fp,"SERIES BC 21 2 0 0\n");
    fprintf(fp,"%-20.2f %-20.2f\n",0.0, cDir_top);
    fprintf(fp,"%-20.2f %-20.2f\n",T_final, cDir_top);
    fprintf(fp,"\n");
    //fprintf(fp,"SERIES BC 10 2 0 0\n");
    //fprintf(fp,"%-20.2f %-20.2f\n",0.0, cDir);
    //fprintf(fp,"%-20.2f %-20.2f\n",T_final, cDir);
    //fprintf(fp,"\n");


    fprintf(fp,"# Rainfall time-series\n");
    fprintf(fp,"SERIES BC 8 2 0 0\n");
    fprintf(fp,"%-20.2f %-20.2f\n",0.0, 1.0);
    fprintf(fp,"%-20.2f %-20.2f\n",T_final, 1.0);
    fprintf(fp,"\n");
    
    fprintf(fp,"TC T0 %-20.2f 0\n",0.0);
    fprintf(fp,"TC TF %-20.2f 0\n",T_final);
    fprintf(fp,"\n");

    fprintf(fp,"#STRINGS ------------------------\n");
    fprintf(fp,"! Set no-flow to West\n");
    fprintf(fp,"NB VEL 2 2 \n"); //!sd changed middle value from 3 to 2
    fprintf(fp,"NB TRN 2 1 5\n");
    fprintf(fp,"! Set Dirichlet transport boundary\n");
    //fprintf(fp,"DB TRN 11 1 8\n");
    fprintf(fp,"DB TRN 12 1 20\n");
    fprintf(fp,"DB TRN 13 1 21\n");
    //fprintf(fp,"! Set tailwater to East\n");
    //fprintf(fp,"NB OTW 3 1\n");
    fprintf(fp,"! Set no-flow to East\n");
    fprintf(fp,"NB TRN 3 1 6\n");
    //fprintf(fp,"! Set global rainfall/evaporation\n");
    //fprintf(fp,"NB SOURCE 1 8\n");
    //fprintf(fp,"NB SOURCE 7 8\n");


    // Node String for Dirichlet boundary
    fprintf(fp,"\n! Dirichlet boundary nodes\n");
    k=1; 
    x = xmin;
    y = ymax;
    for (i=0; i<npx; i++) {
        for (j=0; j<npy; j++) {
            x = xmin + i*dx;
            y = ymax - j*dy; //ymin + j*dy;

            // rotate 45 degrees
            xr = x*cos(theta) - y*sin(theta);
            yr = x*sin(theta) + y*cos(theta);
	    if ( i == 0 || i == npx-1 ) {
            	fprintf(fp,"NDS %d %d \n",k,11); }
	    if ( j == 0 ) {
            	fprintf(fp,"NDS %d %d \n",k,12); }
            if ( j == npy-1 ) {
            	fprintf(fp,"NDS %d %d \n",k,13); }

            k++;
            }
    }
    assert(k == npx*npy + 1);


    fprintf(fp,"\n! West boundary edges\n"); // Originally West boundary edges (before rotation)
    k=0;
    for (i=1; i<npy; i++) {
        //printf("EGS %d %d %d \n",1 + npx * (i - 1), 1 + npx * i, 2);
        fprintf(fp,"EGS %d %d %d \n",i, i + 1, 2);
        k++;
    }
    assert(k==npy-1);

    /*!SD deleted the part below for angled flume domain
    fprintf(fp,"\n! South boundary edges\n");  // Originally South boundary edges (before rotation)
    //fprintf(fp,"NB VEL 3 3\n");
    k=0;
    for (i=1; i<npx; i++) {
        //f2.write('EGS %d %d %d \n' % (i, i+1, 3))
        fprintf(fp,"EGS %d %d %d \n",1 + npy * (i - 1), 1 + npy * i, 3);
        k++;
    }
    assert(k==npx-1);
    */

    fprintf(fp,"\n-------\n");
    fprintf(fp,"\n! East boundary edges\n"); // Originally East boundary edges (before rotation)
    //fprintf(fp,"NB VEL 4 3\n");
    k=0;
    for (i=1; i<npy; i++) {
        //f2.write('EGS %d %d %d \n' % (np_x + np_x * (i - 1), np_x + np_x * i, 4));
        fprintf(fp,"EGS %d %d %d \n",1 + npy * (npx - 1) + (i - 1), 1 + npy * (npx - 1) + i, 3);
        //!SD changed the last output above from 4 to 3
        k++;
    }
    assert(k==npy-1);

    /*!SD deleted the part below for angled flume domain
    fprintf(fp,"\n-------\n");
    //fprintf(fp,"\n! North boundary edges\n"); // Originally North boundary edges (before rotation)
    //fprintf(fp,"NB VEL 5 3\n");
    k=0;
    for (i=1; i<npx; i++) {
        //f2.write('EGS %d %d %d \n' % (1 + np_x * (np_y - 1) + (i - 1), 1 + np_x * (np_y - 1)+i, 5))
        fprintf(fp,"EGS %d %d %d \n",npy + npy * (i - 1), npy + npy * i, 5);
        k++;
    }
    assert(k==npx-1);
    */
    fprintf(fp,"\n");
    fprintf(fp,"! XDMF Output\n");
    fprintf(fp,"PC XDF 1\n");
    fprintf(fp,"\nEND\n");
    fclose(fp);


    fp=fopen("heat.hot", "w"); 
    // ----------water depth initialization (REQUIRED)
    fprintf(fp,"DATASET\n");
    fprintf(fp,"OBJTYPE \"mesh2d\"\n");
    fprintf(fp,"BEGSCL\n");
    fprintf(fp,"ND %d\n",npoints);
    fprintf(fp,"NC %d\n",nelemsTotal);
    fprintf(fp,"NAME IOH\n");
    fprintf(fp,"TS 0 0\n");
    for (i=0; i<npx; i++) {
        for (j=0; j<npy; j++) {
            x = xmin + i*dx;
            y = ymax - j*dy; //ymin + j*dy;

            // rotate 45 degrees
            xr = x*cos(theta) - y*sin(theta);
            yr = x*sin(theta) + y*cos(theta);
            zr = -H0 * xr * xr;
            fprintf(fp,"%10.1f\n",zinit); //initialize the water depth
        }
    }
    fprintf(fp,"ENDDS\n");

    // ----------transport initialization (OPTIONAL)
    fprintf(fp,"DATASET\n");
    fprintf(fp,"OBJTYPE \"mesh2d\"\n");
    fprintf(fp,"BEGSCL\n");
    fprintf(fp,"ND %d\n",npoints);
    fprintf(fp,"NC %d\n",nelemsTotal);
    fprintf(fp,"NAME ICON 1\n");
    fprintf(fp,"TS 0 0\n");
    for (i=0; i<npx; i++) {
        for (j=0; j<npy; j++) {
            x = xmin + i*dx;
            y = ymax - j*dy; //ymin + j*dy;

            // rotate 45 degrees
            xr = x*cos(theta) - y*sin(theta);
            yr = x*sin(theta) + y*cos(theta);
            zr = -H0 * xr * xr;
	    //initialize the constituent concentration
	    double c_profile = cinit_in + (cinit_in - cinit_out)/(xmin-xmax) *(x - xmin);
	    fprintf(fp,"%10.1f\n",c_profile);
//	    if (x < xmax/2) {
//            	fprintf(fp,"%10.1f\n",cinit_in); } //fprintf(fp,"%-20.10f\n",-zr);
//	    else {
//            	fprintf(fp,"%10.1f\n",cinit_out); } //fprintf(fp,"%-20.10f\n",-zr);
        }
    }
    fprintf(fp,"ENDDS\n");
    fclose(fp);


    return 0;
}
