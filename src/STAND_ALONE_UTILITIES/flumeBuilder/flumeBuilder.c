#include "extrusion.h"  // includes AdH global header

char* concat(const char *s1, const char *s2) {
    char *result = malloc(strlen(s1) + strlen(s2) + 1); // +1 for the null-terminator
    // in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

int main(int argc, char *argv[]) {
 
#ifdef _MESSG
    debug_initialize(MPI_COMM_WORLD);
#else
    debug_initialize();
#endif

    if (strcmp(argv[1], "-help") == 0) {
        printf("usage: [root_file_name] [xmin] [xmax] [ymin] [ymax] [npx] [npy] [theta] [dz] [a] [b] [c] [flag_3D]\n");
        exit(0);
    }
    
    SVECT ref_vector;
    ref_vector.x = 0 - 0;
    ref_vector.y = 0 - 0;
    ref_vector.z = -1 - 0;
    
    char *file_name = argv[1];
    double xmin = atof(argv[2]);
    double xmax = atof(argv[3]);
    double ymin = atof(argv[4]);
    double ymax = atof(argv[5]);
    
    int npx = atoi(argv[6]);
    int npy = atof(argv[7]);
    printf("npx: %d\n",npx);
    assert(npx > 1);
    assert(npy > 1);
    if (npx > 2) assert(npx % 2 != 0);
    if (npy > 2) assert(npy % 2 != 0);
    
    double theta = atof(argv[8]);
    assert(theta > -1e-6);
    theta *= 3.141592653589793 / 180.;
    
    double dz = atof(argv[9]);
    assert(dz > 0);
    
    double za = atof(argv[10]);
    double zb = atof(argv[11]);
    double zc = atof(argv[12]);
    int flag_3D = atoi(argv[13]);
    
    char *file_name_2dm =  concat(file_name,".2dm");
    char *file_name_hot =  concat(file_name,".hot");
    char *file_name_bc =   concat(file_name,".bc");
    char *file_name_bin =  concat(file_name,".bin");
    char *file_name_node = concat(file_name,".node");
    
    FILE *fp;
    fp=fopen(file_name_2dm, "w");
    
    fprintf(fp,"MESH2D\n");
    
    double dx = (xmax - xmin)/(double)(npx-1);
    double dy = (ymax - ymin)/(double)(npy-1);
    
    SVECT node[npx*npy];
    
    int i, j, k=1;
    double x = xmin, xr, zr;
    double y = ymin, yr;
    int npoints = 0;
    for (i=0; i<npx; i++) {
        for (j=0; j<npy; j++) {
            x = xmin + i*dx;
            y = ymax - j*dy; //ymin + j*dy;
            
            // totate 45 degrees
            xr = x*cos(theta) - y*sin(theta);
            yr = x*sin(theta) + y*cos(theta);
            zr = -(za + zb * x + zc * x * x);
            fprintf(fp,"ND %d %20.10f %20.10f %20.10f\n",k,xr,yr,zr);
            
            node[npoints].x = xr;
            node[npoints].y = yr;
            node[npoints].z = -(za + zb * x + zc * x * x);
            
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
            
            // original output
            fprintf(fp,"E3T %d %d %d %d 1\n",k+1,nd1+1,nd2+1,nd3+1);
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
            fprintf(fp,"E3T %d %d %d %d 1\n",k+1,nd1+1,nd2+1,nd3+1);
            k++;
            
            //            if ((i+1) % 2) { // odd
            //                if ((igrp+1) % 2) { // odd
            //                    fprintf(fp,"E3T %d %d %d %d 1\n",k+1,node+1,node+npy+1,node+npy+1+1);
            //                    k++;
            //                    fprintf(fp,"E3T %d %d %d %d 1\n",k+1,node+1,node+npy+1+1,node+1+1);
            //                    k++;
            //                } else {
            //                    fprintf(fp,"E3T %d %d %d %d 1\n",k+1,node+1,node+npy+1,node+1+1);
            //                    k++;
            //                    fprintf(fp,"E3T %d %d %d %d 1\n",k+1,node+2,node+npy+1,node+npy+1+1);
            //                    k++;
            //                }
            //            } else {
            //                if ((igrp+1) % 2) { // odd
            //                    fprintf(fp,"E3T %d %d %d %d 1\n",k+1,node+1,node+npy+1,node+1+1);
            //                    k++;
            //                    fprintf(fp,"E3T %d %d %d %d 1\n",k+1,node+2,node+npy+1,node+npy+1+1);
            //                    k++;
            //                } else {
            //                    fprintf(fp,"E3T %d %d %d %d 1\n",k+1,node+1,node+npy+1,node+npy+1+1);
            //                    k++;
            //                    fprintf(fp,"E3T %d %d %d %d 1\n",k+1,node+1,node+npy+1+1,node+1+1);
            //                    k++;
            //                }
            //            }
        }
    }
    int nelemsTotal = k;
    
    fclose(fp);
    fp=fopen(file_name_bc, "w");
    fprintf(fp,"# operational parameters\n");
    fprintf(fp,"OP SW2\n");
    fprintf(fp,"OP TRN 0\n");
    fprintf(fp,"OP BLK 1\n");
    fprintf(fp,"OP INC 40\n");
    fprintf(fp,"OP PRE 1\n");
    fprintf(fp,"\n");
    fprintf(fp,"# solver parameters\n");
    fprintf(fp,"IP NTL 1e-6\n");
    fprintf(fp,"IP ITL 1e-6\n");
    fprintf(fp,"IP NIT 10\n");
    fprintf(fp,"IP MIT 100\n");
    fprintf(fp,"\n");
    fprintf(fp,"# 2d element material string\n");
    fprintf(fp,"MTS 1 1\n");
    fprintf(fp,"\n");
    fprintf(fp,"# material properties\n");
    fprintf(fp,"MP ML 1 0\n");
    fprintf(fp,"MP SRT 1 100\n");
    fprintf(fp,"MP EVS 1 0.0 0.0 0.0\n");
    fprintf(fp,"MP MUC 1.0\n");
    fprintf(fp,"MP MU 1e-6\n");
    fprintf(fp,"MP RHO 1000\n");
    fprintf(fp,"MP G 9.8\n");
    fprintf(fp,"\n");
    fprintf(fp,"# friction coefficients\n");
    fprintf(fp,"FR MNG 1 0.0\n");
    fprintf(fp,"\n");
    fprintf(fp,"#TIME-SERIES --------------------\n\n");
    fprintf(fp,"# Output time-series \n");
    fprintf(fp,"SERIES AWRITE 1 1 0 0\n");
    fprintf(fp,"%-20.10f %-20.10f %-20.10f 0\n",0.,1000.,100.);
    fprintf(fp,"\n");
    fprintf(fp,"# Time-step time-series\n");
    fprintf(fp,"SERIES DT 2 2 0 0\n");
    fprintf(fp,"%-20.10f %-20.10f\n",0.,100.);
    fprintf(fp,"%-20.10f %-20.10f\n",1000.,100.);
    fprintf(fp,"\n");
    fprintf(fp,"# No-flow time-series\n");
    fprintf(fp,"SERIES BC 3 2 0 0\n");
    fprintf(fp,"%-20.10f 0.0\n",0.);
    fprintf(fp,"%-20.10f 0.0\n",1000.);
    fprintf(fp,"\n");
    fprintf(fp,"TC T0 %-20.10f\n",0.);
    fprintf(fp,"TC TF %-20.10f\n",1000.);
    fprintf(fp,"\n");
    fprintf(fp,"#STRINGS ------------------------\n");
    fprintf(fp,"! West\n");
    fprintf(fp,"NB VEL 2 3\n");
    k=0;
    for (i=1; i<npy; i++) {
        //printf("EGS %d %d %d \n",1 + npx * (i - 1), 1 + npx * i, 2);
        fprintf(fp,"EGS %d %d %d \n",i, i + 1, 2);
        k++;
    }
    assert(k==npy-1);
    
    fprintf(fp,"! South\n");
    fprintf(fp,"NB VEL 3 3\n");
    k=0;
    for (i=1; i<npx; i++) {
        //f2.write('EGS %d %d %d \n' % (i, i+1, 3))
        fprintf(fp,"EGS %d %d %d \n",1 + npy * (i - 1), 1 + npy * i, 3);
        k++;
    }
    assert(k==npx-1);
    
    fprintf(fp,"! East\n");
    fprintf(fp,"NB VEL 4 3\n");
    k=0;
    for (i=1; i<npy; i++) {
        //f2.write('EGS %d %d %d \n' % (np_x + np_x * (i - 1), np_x + np_x * i, 4));
        fprintf(fp,"EGS %d %d %d \n",1 + npy * (npx - 1) + (i - 1), 1 + npy * (npx - 1) + i, 4);
        k++;
    }
    assert(k==npy-1);
    
    fprintf(fp,"! North\n");
    fprintf(fp,"NB VEL 5 3\n");
    k=0;
    for (i=1; i<npx; i++) {
        //f2.write('EGS %d %d %d \n' % (1 + np_x * (np_y - 1) + (i - 1), 1 + np_x * (np_y - 1)+i, 5))
        fprintf(fp,"EGS %d %d %d \n",npy + npy * (i - 1), npy + npy * i, 5);
        k++;
    }
    assert(k==npx-1);
    fprintf(fp,"\nEND\n");
    fclose(fp);
    printf("bc file written\n");
    
    
    fp=fopen(file_name_hot, "w");
    fprintf(fp,"DATASET\n");
    fprintf(fp,"OBJTYPE \"mesh2d\"\n");
    fprintf(fp,"BEGSCL\n");
    fprintf(fp,"ND %d\n",npoints);
    fprintf(fp,"NC %d\n",nelemsTotal);
    fprintf(fp,"NAME ioh\n");
    fprintf(fp,"TS 0 0\n");
    for (i=0; i<npx; i++) {
        for (j=0; j<npy; j++) {
            x = xmin + i*dx;
            y = ymax - j*dy; //ymin + j*dy;
            
            // totate 45 degrees
            xr = x*cos(theta) - y*sin(theta);
            yr = x*sin(theta) + y*cos(theta);
            zr = -(za + zb * x + zc * x * x);
            fprintf(fp,"%-20.10f\n",-zr);
        }
    }
    fprintf(fp,"ENDDS\n");
    fclose(fp);
    printf("hotstart written\n");
    
    
    if (flag_3D) {
        fp=fopen(file_name_bin, "w");
        fprintf(fp,"BIN 1 1 0. %f \n",dz);
        fprintf(fp,"BIN 1 2 -1000 %f \n",dz);
        fclose(fp);
        printf("extrusion bin file written\n");
        
        fp=fopen(file_name_node, "w");
        for (i=0; i< npoints; i++) {
            fprintf(fp,"ND %d 1 \n",i+1);
        }
        fclose(fp);
        printf("extrusion node file written\n");
        
        
        free(file_name_2dm);
        free(file_name_bc);
        free(file_name_node);
        free(file_name_hot);
        free(file_name_bin);
        
        
        // NOW CALL 3D EXTRUSION
        printf("extruding files now ... \n");
        int flag_bins = 1;
        int flag_min_layer = 0;
        int minimum_layers = 1;
        int mesh_type = TETRAHEDRON; //MIXED_ELEMENT_MESH;
        extrudeAdH(flag_bins,flag_min_layer,minimum_layers,file_name,mesh_type);

    }
    
    return 0;
}


