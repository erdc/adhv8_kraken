#include "global_header.h"

void create_velocity_file_test(SGRID *grid, double t1, double t2, int nt, char *file_base, SVECT (*get_vel)(SVECT p, double t));
void create_water_quality_file_test(SGRID *grid, double t1, double t2, int nt, char *file_base, SVECT (*get_wq)(SVECT p, double t));
void create_nd_dpl_file_test(SGRID *grid, double t1, double t2, int nt, char *file_base);
void error_report(char *file, int line, int right_elem, int ielem) {
    printf("ERROR :: file: %s :: line: %d :: correct element: %d :: element_found: %d\n",file,line,right_elem,ielem);
    exit(-1);
}

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

int pt_test_searchEngine() {
    int i, ielem, flag, call_ID=0, node_id;
    SVECT p1, p2, pi, pi0;
    double weights[4];
    pi0.x = -99999.0;
    pi0.y = -99999.0;
    
    RUN_VERBOSE = false;
 
     // allocate hydro grid
    SGRID *grid = (SGRID *)malloc(sizeof(SGRID));

    double t0 = 0.0;
    double tf = 100.0;
    double nt = 50.;
    char *filebase = "teset";
    int test_flag = 1;

    if (test_flag == 1) {
        
        // ----------------------------------------
        // ----------------------------------------
        //printf("**** FULL TEST ON ****\n");
        // ----------------------------------------
        
        // build 2D test grid
        sgrid_create2d(grid,0,100,0,100,21,21,0,200,100,0,0);
        
        // ----------------------------------------
        // ----------------------------------------
        printf("-- TESTING GEOMETRY SEARCH ENGINE\n");
        test_geometry();
        
        // ----------------------------------------
        // ----------------------------------------
        printf("-- TESTING 2D VELOCITY PROJECTION ");
        SVECT v1,v2,v1_proj,v2_proj,nd1,nd2;
        v1.x =  1.0; v1.y = 1.0;     // along the segment (give back original velocity)
        v2.x = -1.0; v2.y = 0.0;     // in the direction opposite the segment normal (give back original velocity)
        nd1.x = 0.0; nd1.y = 0.0;
        nd2.x = 5.0; nd2.y = 5.0;
        svect_init(&v1_proj); svect_init(&v2_proj);
        project_velocity_2D(v1,v2,nd1,nd2,&v1_proj,&v2_proj);
        if(fabs(v1_proj.x-v1.x) > 1e-10 || fabs(v1_proj.y-v1.y) > 1e-10)  {
            printf("TEST ERROR || FILE: %s || LINE: %d \n",__FILE__,__LINE__);
            printf("p1: %f %f p2: %f %f || v1: %f %f v1_proj: %f %f || v2: %f %f v2_proj: %f %f\n",
                   nd1.x,nd1.y,nd2.x,nd2.y,v1.x,v1.y,v1_proj.x,v1_proj.y,v2.x,v2.y,v2_proj.x,v2_proj.y);
            exit(-1);
        }
        if(fabs(v2_proj.x-v2.x) > 1e-10 || fabs(v2_proj.y-v2.y) > 1e-10)  {
                  printf("TEST ERROR || FILE: %s || LINE: %d \n",__FILE__,__LINE__);
                  printf("p1: %f %f p2: %f %f || v1: %f %f v1_proj: %f %f || v2: %f %f v2_proj: %f %f\n",
                         nd1.x,nd1.y,nd2.x,nd2.y,v1.x,v1.y,v1_proj.x,v1_proj.y,v2.x,v2.y,v2_proj.x,v2_proj.y);
                  exit(-1);
              }

        v1.x = 1.0; v1.y = -2.5;   // this will be projected
        v2.x = 1.0; v2.y = -1.0;   // this will be projected, but orhongal to edge so gives v = 0,0
        nd1.x = 0.0; nd1.y = 0.0;
        nd2.x = 5.0; nd2.y = 5.0;
        svect_init(&v1_proj); svect_init(&v2_proj);
        project_velocity_2D(v1,v2,nd1,nd2,&v1_proj,&v2_proj);
        if(fabs(v1_proj.x+0.75) > 1e-10 || fabs(v1_proj.y+0.75) > 1e-10)  {
                  printf("TEST ERROR || FILE: %s || LINE: %d \n",__FILE__,__LINE__);
                  printf("p1: %f %f p2: %f %f || v1: %f %f v1_proj: %f %f || v2: %f %f v2_proj: %f %f\n",
                         nd1.x,nd1.y,nd2.x,nd2.y,v1.x,v1.y,v1_proj.x,v1_proj.y,v2.x,v2.y,v2_proj.x,v2_proj.y);
                  exit(-1);
              }
        if(fabs(v2_proj.x) > 1e-10 || fabs(v2_proj.y) > 1e-10)   {
                  printf("TEST ERROR || FILE: %s || LINE: %d \n",__FILE__,__LINE__);
                  printf("p1: %f %f p2: %f %f || v1: %f %f v1_proj: %f %f || v2: %f %f v2_proj: %f %f\n",
                         nd1.x,nd1.y,nd2.x,nd2.y,v1.x,v1.y,v1_proj.x,v1_proj.y,v2.x,v2.y,v2_proj.x,v2_proj.y);
                  exit(-1);
              }
        //printf("p1: %f %f p2: %f %f || v1: %f %f v1_proj: %f %f || v2: %f %f v2_proj: %f %f\n",
        //       nd1.x,nd1.y,nd2.x,nd2.y,v1.x,v1.y,v1_proj.x,v1_proj.y,v2.x,v2.y,v2_proj.x,v2_proj.y);
        printf("|| PASSED\n");
        
        // ----------------------------------------
        // ----------------------------------------
        printf("-- TESTING 3D VELOCITY PROJECTION ");
        SVECT v3,v3_proj,nd3;
        v1.x =  1.0; v1.y = 1.0; v1.z = 0.0;    // along the face (give back original velocity)
        v2.x = -1.0; v2.y = 0.0; v2.z = 0.0;    // in the direction opposite the face normal (give back original velocity)
        v3.x =  0.0; v3.y = 0.0; v3.z = 1.0;    // in the direction of face normal, but perpendicular (gives back 0)
        nd1.x = 0.0; nd1.y = 0.0; nd1.z = 0.0;
        nd2.x = 5.0; nd2.y = 0.0; nd2.z = 0.0;
        nd3.x = 5.0; nd3.y = 5.0; nd3.z = 0.0;
        svect_init(&v1_proj); svect_init(&v2_proj); svect_init(&v3_proj);
        project_velocity_3D(v1,v2,v3,nd1,nd2,nd3,&v1_proj,&v2_proj,&v3_proj);
        if(fabs(v1_proj.x-v1.x) > 1e-10 || fabs(v1_proj.y-v1.y) > 1e-10 || fabs(v1_proj.z-v1.z) > 1e-10)  {
            printf("TEST ERROR || FILE: %s || LINE: %d \n",__FILE__,__LINE__);
            printf("p1: %f %f %f p2: %f %f %f || p3: %f %f %f || v1: %f %f %f v1_proj: %f %f %f || v2: %f %f %f v2_proj: %f %f %f || v3: %f %f %f v3_proj: %f %f %f\n",
                   nd1.x,nd1.y,nd1.z,nd2.x,nd2.y,nd2.z,nd3.x,nd3.y,nd3.z,
                   v1.x,v1.y,v1.z,v1_proj.x,v1_proj.y,v1_proj.z,
                   v2.x,v2.y,v2.z,v2_proj.x,v2_proj.y,v2_proj.z,
                   v3.x,v3.y,v3.z,v3_proj.x,v3_proj.y,v3_proj.z);
            exit(-1);
        }
        if(fabs(v2_proj.x-v2.x) > 1e-10 || fabs(v2_proj.y-v2.y) > 1e-10 || fabs(v2_proj.z-v2.z) > 1e-10)  {
            printf("TEST ERROR || FILE: %s || LINE: %d \n",__FILE__,__LINE__);
            printf("p1: %f %f %f p2: %f %f %f || p3: %f %f %f || v1: %f %f %f v1_proj: %f %f %f || v2: %f %f %f v2_proj: %f %f %f || v3: %f %f %f v3_proj: %f %f %f\n",
                   nd1.x,nd1.y,nd1.z,nd2.x,nd2.y,nd2.z,nd3.x,nd3.y,nd3.z,
                   v1.x,v1.y,v1.z,v1_proj.x,v1_proj.y,v1_proj.z,
                   v2.x,v2.y,v2.z,v2_proj.x,v2_proj.y,v2_proj.z,
                   v3.x,v3.y,v3.z,v3_proj.x,v3_proj.y,v3_proj.z);
            exit(-1);
        }
        
        if(fabs(v3_proj.x) > 1e-10 || fabs(v3_proj.y) > 1e-10 || fabs(v3_proj.z) > 1e-10)  {
            printf("TEST ERROR || FILE: %s || LINE: %d \n",__FILE__,__LINE__);
            printf("p1: %f %f %f p2: %f %f %f || p3: %f %f %f || v1: %f %f %f v1_proj: %f %f %f || v2: %f %f %f v2_proj: %f %f %f || v3: %f %f %f v3_proj: %f %f %f\n",
                   nd1.x,nd1.y,nd1.z,nd2.x,nd2.y,nd2.z,nd3.x,nd3.y,nd3.z,
                   v1.x,v1.y,v1.z,v1_proj.x,v1_proj.y,v1_proj.z,
                   v2.x,v2.y,v2.z,v2_proj.x,v2_proj.y,v2_proj.z,
                   v3.x,v3.y,v3.z,v3_proj.x,v3_proj.y,v3_proj.z);
            exit(-1);
        }
        v1.x =  1.0; v1.y = 1.0; v1.z =  1.0;    // Rides an edge, gives same velocity back
        v2.x =  0.0; v2.y = 0.0; v2.z = -1.0;    // inflow, give velocity back
        v3.x = -1.0; v3.y = 0.0; v3.z =  0.0;    // project!
        nd1.x = 0.0; nd1.y = 0.0; nd1.z = 0.0;
        nd2.x = 1.0; nd2.y = 0.0; nd2.z = 1.0;
        nd3.x = 1.0; nd3.y = 1.0; nd3.z = 1.0;
        svect_init(&v1_proj); svect_init(&v2_proj); svect_init(&v3_proj);
        project_velocity_3D(v1,v2,v3,nd1,nd2,nd3,&v1_proj,&v2_proj,&v3_proj);
        if(fabs(v1_proj.x-v1.x) > 1e-10 || fabs(v1_proj.y-v1.y) > 1e-10 || fabs(v1_proj.z-v1.z) > 1e-10)  {
            printf("TEST ERROR || FILE: %s || LINE: %d \n",__FILE__,__LINE__);
            printf("p1: %f %f %f p2: %f %f %f || p3: %f %f %f || v1: %f %f %f v1_proj: %f %f %f || v2: %f %f %f v2_proj: %f %f %f || v3: %f %f %f v3_proj: %f %f %f\n",
                   nd1.x,nd1.y,nd1.z,nd2.x,nd2.y,nd2.z,nd3.x,nd3.y,nd3.z,
                   v1.x,v1.y,v1.z,v1_proj.x,v1_proj.y,v1_proj.z,
                   v2.x,v2.y,v2.z,v2_proj.x,v2_proj.y,v2_proj.z,
                   v3.x,v3.y,v3.z,v3_proj.x,v3_proj.y,v3_proj.z);
            exit(-1);
        }
        if(fabs(v2_proj.x-v2.x) > 1e-10 || fabs(v2_proj.y-v2.y) > 1e-10 || fabs(v2_proj.z-v2.z) > 1e-10)  {
            printf("TEST ERROR || FILE: %s || LINE: %d \n",__FILE__,__LINE__);
            printf("p1: %f %f %f p2: %f %f %f || p3: %f %f %f || v1: %f %f %f v1_proj: %f %f %f || v2: %f %f %f v2_proj: %f %f %f || v3: %f %f %f v3_proj: %f %f %f\n",
                   nd1.x,nd1.y,nd1.z,nd2.x,nd2.y,nd2.z,nd3.x,nd3.y,nd3.z,
                   v1.x,v1.y,v1.z,v1_proj.x,v1_proj.y,v1_proj.z,
                   v2.x,v2.y,v2.z,v2_proj.x,v2_proj.y,v2_proj.z,
                   v3.x,v3.y,v3.z,v3_proj.x,v3_proj.y,v3_proj.z);
            exit(-1);
        }
        
        if(fabs(v3_proj.x+0.5) > 1e-10 || fabs(v3_proj.y) > 1e-10 || fabs(v3_proj.z+0.5) > 1e-10)  {
            printf("TEST ERROR || FILE: %s || LINE: %d \n",__FILE__,__LINE__);
            printf("p1: %f %f %f p2: %f %f %f || p3: %f %f %f || v1: %f %f %f v1_proj: %f %f %f || v2: %f %f %f v2_proj: %f %f %f || v3: %f %f %f v3_proj: %f %f %f\n",
                   nd1.x,nd1.y,nd1.z,nd2.x,nd2.y,nd2.z,nd3.x,nd3.y,nd3.z,
                   v1.x,v1.y,v1.z,v1_proj.x,v1_proj.y,v1_proj.z,
                   v2.x,v2.y,v2.z,v2_proj.x,v2_proj.y,v2_proj.z,
                   v3.x,v3.y,v3.z,v3_proj.x,v3_proj.y,v3_proj.z);
            exit(-1);
        }
        printf("|| PASSED\n");
        
        // ----------------------------------------
        // ----------------------------------------
        printf("-- TESTING 2D SEARCH ENGINE\n");
        
        double dx = 5.0, dy = 5.0, x0 = 50.0, y0 = 50.0;
        double xc = x0 - dx/4.0, x2 = x0 + dx/4.0;
        double yc = y0 - dy/4.0, y2 = y0 + dy/4.0;
        double zc = 0.0;

        //----------------------------------------------------------------------------------------------
        printf("---- TESTING 2D ENGINE NODE START\n");
        
        p1.x = x0;  p1.y = y0; p1.z = -100.0;
        printf("-------- TEST :: particle trajectory starts on node and ends on same element center :: ");
        p2.x = xc;  p2.y = yc; p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 381) error_report(__FILE__,__LINE__,381,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and ends on same element node :: ");
        p2.x = x0-dx;  p2.y = y0; p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 381) error_report(__FILE__,__LINE__,381,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and ends on same element edge :: ");
        p2.x = xc;  p2.y = y0; p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 381) error_report(__FILE__,__LINE__,381,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and ends on next element center :: ");
        p2.x = (1./3.) * (grid->node[ grid->elem2d[418].nodes[0]].x +
                          grid->node[ grid->elem2d[418].nodes[1]].x +
                          grid->node[ grid->elem2d[418].nodes[2]].x);
        p2.y = (1./3.) * (grid->node[ grid->elem2d[418].nodes[0]].y +
                          grid->node[ grid->elem2d[418].nodes[1]].y +
                          grid->node[ grid->elem2d[418].nodes[2]].y);
        p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 418) error_report(__FILE__,__LINE__,418,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and ends on next element node :: ");
        p2.x = x0+dx;  p2.y = y0;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 418 && ielem != 421) error_report(__FILE__,__LINE__,418,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and ends on next element edge :: ");
        p2.x = x0-dx/2.0; p2.y = y0-dy, p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 380 && ielem != 383) error_report(__FILE__,__LINE__,380,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and rides grid to top-left grid corner :: ");
        p2.x = 0;   p2.y = 100;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 0 && ielem != 1) error_report(__FILE__,__LINE__,1,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and rides grid to bottom-left grid corner :: ");
        p2.x = 0;   p2.y = 0;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 38 && ielem != 39) error_report(__FILE__,__LINE__,1,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and rides grid to top-right grid corner :: ");
        p2.x = 100;   p2.y = 100;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 760 && ielem != 761) error_report(__FILE__,__LINE__,1,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and rides grid to bottom-right grid corner :: ");
        p2.x = 100;   p2.y = 0;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 798 && ielem != 799) error_report(__FILE__,__LINE__,1,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and leaves grid @ node :: ");
        p2.x = 110.0;   p2.y = 110.0;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 760 && ielem != 761) error_report(__FILE__,__LINE__,1,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and leaves grid @ edge :: ");
        p2.x = 99.85;   p2.y = -.5;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 798 && ielem != 799) error_report(__FILE__,__LINE__,1,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        
        //----------------------------------------------------------------------------------------------
        printf("---- TESTING 2D ENGINE ELEMENT CENTER START\n");
        p1.x = xc; p1.y = yc; p1.z = -100.0;
        printf("-------- TEST :: particle trajectory starts on center and ends on same element center :: ");
        p2.x = xc;  p2.y = yc; p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 381) error_report(__FILE__,__LINE__,381,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on center and ends on same element node :: ");
        p2.x = x0-dx;  p2.y = y0; p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 381) error_report(__FILE__,__LINE__,381,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on center and ends on same element edge :: ");
        p2.x = xc;  p2.y = y0; p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 381) error_report(__FILE__,__LINE__,381,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on center and ends on next element center :: ");
        p2.x = (1./3.) * (grid->node[ grid->elem2d[418].nodes[0]].x +
                          grid->node[ grid->elem2d[418].nodes[1]].x +
                          grid->node[ grid->elem2d[418].nodes[2]].x);
        p2.y = (1./3.) * (grid->node[ grid->elem2d[418].nodes[0]].y +
                          grid->node[ grid->elem2d[418].nodes[1]].y +
                          grid->node[ grid->elem2d[418].nodes[2]].y);
        p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 418) error_report(__FILE__,__LINE__,418,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST ::  particle trajectory starts on center and ends on next element node :: ");
        p2.x = x0+dx;  p2.y = y0;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 418 && ielem != 421) error_report(__FILE__,__LINE__,418,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on center and ends on next element edge :: ");
        p2.x = x0-dx/2.0; p2.y = y0-dy, p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 380 && ielem != 383) error_report(__FILE__,__LINE__,380,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on center and ends on corner of grid :: ");
        p2.x = 0;   p2.y = 0;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 38 && ielem != 39) error_report(__FILE__,__LINE__,39,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on center and rides edges all the way to end of grid :: ");
        p2.x = 0;   p2.y = 100;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 0 && ielem !=1) error_report(__FILE__,__LINE__,1,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on center and leaves grid @ node :: ");
        p2.x = 110.0;   p2.y = 110.0;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 760 && ielem != 761) error_report(__FILE__,__LINE__,1,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on center and leaves grid @ edge :: ");
        p2.x = 99.85;   p2.y = -.5;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 798 && ielem != 799) error_report(__FILE__,__LINE__,1,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        
        //----------------------------------------------------------------------------------------------
        printf("---- TESTING 2D ENGINE ELEMENT EDGE START\n");
        p1.x = xc;  p1.y = y0;   p1.z = -100.0;
        printf("-------- TEST :: particle trajectory starts on edge and ends on same element center :: ");
        p2.x = xc;  p2.y = yc; p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 381) error_report(__FILE__,__LINE__,381,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on edge and ends on same element node :: ");
        p2.x = x0-dx;  p2.y = y0; p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 381) error_report(__FILE__,__LINE__,381,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on edge and ends on same element edge :: ");
        p2.x = xc;  p2.y = y0; p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 381) error_report(__FILE__,__LINE__,381,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on edge and ends on next element center :: ");
        p2.x = (1./3.) * (grid->node[ grid->elem2d[418].nodes[0]].x +
                          grid->node[ grid->elem2d[418].nodes[1]].x +
                          grid->node[ grid->elem2d[418].nodes[2]].x);
        p2.y = (1./3.) * (grid->node[ grid->elem2d[418].nodes[0]].y +
                          grid->node[ grid->elem2d[418].nodes[1]].y +
                          grid->node[ grid->elem2d[418].nodes[2]].y);
        p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 418) error_report(__FILE__,__LINE__,418,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST ::  particle trajectory starts on edge and ends on next element node :: ");
        p2.x = x0+dx;  p2.y = y0;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 418 && ielem != 421) error_report(__FILE__,__LINE__,418,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on edge and ends on next element edge :: ");
        p2.x = x0-dx/2.0; p2.y = y0-dy, p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 380 && ielem != 383) error_report(__FILE__,__LINE__,380,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on edge and ends on corner of grid :: ");
        p2.x = 0;   p2.y = 0;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 38 && ielem != 39) error_report(__FILE__,__LINE__,39,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on edge and rides edges all the way to end of grid :: ");
        p2.x = 0;   p2.y = 100;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 0 && ielem !=1) error_report(__FILE__,__LINE__,1,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on edge and leaves grid @ node :: ");
        p2.x = 110.0;   p2.y = 110.0;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 760 && ielem != 761) error_report(__FILE__,__LINE__,1,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on edge and leaves grid @ edge :: ");
        p2.x = 99.85;   p2.y = -.5;  p2.z = -100.0;
        ielem = 381; flag = elementSearch2D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 798 && ielem != 799) error_report(__FILE__,__LINE__,1,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        
        // ----------------------------------------
        printf("-- TESTING 3D SEARCH ENGINE\n");
        
        // build 2D test grid
        sgrid_free(grid);
        grid = (SGRID *)malloc(sizeof(SGRID));
        sgrid_read_adh(&grid,"test_3d");
        
        // pick a starting element
        int ielem0 = 11446;
        int ielem_adj = 11355;
        
        xc =   (grid->node[grid->elem3d[ielem0].nodes[0]].x +
                grid->node[grid->elem3d[ielem0].nodes[1]].x +
                grid->node[grid->elem3d[ielem0].nodes[2]].x +
                grid->node[grid->elem3d[ielem0].nodes[3]].x) / 4.0;
        yc =   (grid->node[grid->elem3d[ielem0].nodes[0]].y +
                grid->node[grid->elem3d[ielem0].nodes[1]].y +
                grid->node[grid->elem3d[ielem0].nodes[2]].y +
                grid->node[grid->elem3d[ielem0].nodes[3]].y) / 4.0;
        zc =   (grid->node[grid->elem3d[ielem0].nodes[0]].z +
                grid->node[grid->elem3d[ielem0].nodes[1]].z +
                grid->node[grid->elem3d[ielem0].nodes[2]].z +
                grid->node[grid->elem3d[ielem0].nodes[3]].z) / 4.0;
        
        assert(grid->elem3d[ielem0].nodes[1] == 2425);
        
        //----------------------------------------------------------------------------------------------
        printf("---- TESTING 3D ENGINE NODE START\n");
        p1.x = grid->node[grid->elem3d[ielem0].nodes[1]].x;
        p1.y = grid->node[grid->elem3d[ielem0].nodes[1]].y;
        p1.z = grid->node[grid->elem3d[ielem0].nodes[1]].z;
        printf("-------- TEST :: particle trajectory starts on node and ends on same element center :: ");
        p2.x = xc; p2.y = yc; p2.z = zc;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem0) error_report(__FILE__,__LINE__,ielem0,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and ends on same element node :: ");
        p2.x = grid->node[grid->elem3d[ielem0].nodes[2]].x;
        p2.y = grid->node[grid->elem3d[ielem0].nodes[2]].y;
        p2.z = grid->node[grid->elem3d[ielem0].nodes[2]].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != grid->elem3d[ielem0].nodes[2] &&
           grid->elem3d[ielem].nodes[1] != grid->elem3d[ielem0].nodes[2] &&
           grid->elem3d[ielem].nodes[2] != grid->elem3d[ielem0].nodes[2] &&
           grid->elem3d[ielem].nodes[3] != grid->elem3d[ielem0].nodes[2]) error_report(__FILE__,__LINE__,UNSET_INT,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and ends on same element edge :: ");
        p2.x = grid->node[grid->elem3d[ielem0].nodes[1]].x;
        p2.y = grid->node[grid->elem3d[ielem0].nodes[1]].y;
        p2.z = 0.5 * (grid->node[grid->elem3d[ielem0].nodes[1]].z + grid->node[grid->elem3d[ielem0].nodes[3]].z);
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem0) error_report(__FILE__,__LINE__,ielem0,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and ends on same element face :: ");
        p2.x = (1.0/3.0) * (grid->node[grid->elem3d[ielem0].nodes[0]].x + grid->node[grid->elem3d[ielem0].nodes[2]].x + grid->node[grid->elem3d[ielem0].nodes[3]].x);
        p2.y = (1.0/3.0) * (grid->node[grid->elem3d[ielem0].nodes[0]].y + grid->node[grid->elem3d[ielem0].nodes[2]].y + grid->node[grid->elem3d[ielem0].nodes[3]].y);
        p2.z = (1.0/3.0) * (grid->node[grid->elem3d[ielem0].nodes[0]].z + grid->node[grid->elem3d[ielem0].nodes[2]].z + grid->node[grid->elem3d[ielem0].nodes[3]].z);
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem0) error_report(__FILE__,__LINE__,ielem0,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and ends on next element center:: ");
        p2.x = (grid->node[grid->elem3d[ielem_adj].nodes[0]].x +
                grid->node[grid->elem3d[ielem_adj].nodes[1]].x +
                grid->node[grid->elem3d[ielem_adj].nodes[2]].x +
                grid->node[grid->elem3d[ielem_adj].nodes[3]].x) / 4.0;
        p2.y = (grid->node[grid->elem3d[ielem_adj].nodes[0]].y +
                grid->node[grid->elem3d[ielem_adj].nodes[1]].y +
                grid->node[grid->elem3d[ielem_adj].nodes[2]].y +
                grid->node[grid->elem3d[ielem_adj].nodes[3]].y) / 4.0;
        p2.z = (grid->node[grid->elem3d[ielem_adj].nodes[0]].z +
                grid->node[grid->elem3d[ielem_adj].nodes[1]].z +
                grid->node[grid->elem3d[ielem_adj].nodes[2]].z +
                grid->node[grid->elem3d[ielem_adj].nodes[3]].z) / 4.0;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem_adj) error_report(__FILE__,__LINE__,ielem_adj,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and ends on next element node :: ");
        p2.x = grid->node[grid->elem3d[ielem_adj].nodes[2]].x;
        p2.y = grid->node[grid->elem3d[ielem_adj].nodes[2]].y;
        p2.z = grid->node[grid->elem3d[ielem_adj].nodes[2]].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != grid->elem3d[ielem_adj].nodes[2] &&
           grid->elem3d[ielem].nodes[1] != grid->elem3d[ielem_adj].nodes[2] &&
           grid->elem3d[ielem].nodes[2] != grid->elem3d[ielem_adj].nodes[2] &&
           grid->elem3d[ielem].nodes[3] != grid->elem3d[ielem_adj].nodes[2]) error_report(__FILE__,__LINE__,UNSET_INT,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and ends on next element edge :: ");
        p2.x = (grid->node[grid->elem3d[ielem_adj].nodes[2]].x + grid->node[grid->elem3d[ielem0].nodes[1]].x)/2.0;
        p2.y = (grid->node[grid->elem3d[ielem_adj].nodes[2]].y + grid->node[grid->elem3d[ielem0].nodes[1]].y)/2.0;
        p2.z = grid->node[grid->elem3d[ielem0].nodes[1]].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 11354) error_report(__FILE__,__LINE__,ielem_adj,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and ends on next element face :: ");
        p2.x = (1.0/3.0) * (grid->node[grid->elem3d[ielem_adj].nodes[0]].x + grid->node[grid->elem3d[ielem_adj].nodes[2]].x + grid->node[grid->elem3d[ielem_adj].nodes[3]].x);
        p2.y = (1.0/3.0) * (grid->node[grid->elem3d[ielem_adj].nodes[0]].y + grid->node[grid->elem3d[ielem_adj].nodes[2]].y + grid->node[grid->elem3d[ielem_adj].nodes[3]].y);
        p2.z = (1.0/3.0) * (grid->node[grid->elem3d[ielem_adj].nodes[0]].z + grid->node[grid->elem3d[ielem_adj].nodes[2]].z + grid->node[grid->elem3d[ielem_adj].nodes[3]].z);
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem_adj) error_report(__FILE__,__LINE__,ielem_adj,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and ends on corner of grid :: ");
        p2.x = 0.0; p2.y = 0.0; p2.z = -100.0;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 1168 && ielem != 1169) error_report(__FILE__,__LINE__,1168,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and rides edges all the way to end of grid :: ");
        p2.x = 0.0; p2.y = 50.0; p2.z = -50.0;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 553) error_report(__FILE__,__LINE__,553,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and ends on external boundary face :: ");
        p2.x = 0.0; p2.y = 0.1; p2.z = -99.9;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 1169) error_report(__FILE__,__LINE__,1169,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        
        printf("-------- TEST :: particle trajectory starts on node and rides grid to top-left grid corner :: ");
        node_id = 10;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and rides grid to bottom-left grid corner :: ");
        node_id = 230;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and rides grid to top-right grid corner :: ");
        node_id = 4630;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and rides grid to bottom-right grid corner :: ");
        node_id = 4850;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and leaves grid @ node :: ");
        node_id = 4630;
        p2.x = grid->node[node_id].x + 5.0;   p2.y = grid->node[node_id].y + 5.0;  p2.z = grid->node[node_id].z - 5.0;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and leaves grid @ edge :: ");
        node_id = 10;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z - .5;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and leaves grid @ face :: ");
        node_id = 10;
        p2.x = grid->node[node_id].x-1;   p2.y = grid->node[node_id].y-1;  p2.z = grid->node[node_id].z + .5;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);

        //----------------------------------------------------------------------------------------------
        printf("---- TESTING 3D ENGINE ELEMENT CENTER START\n");
        p1.x = xc; p1.y = yc; p1.z = zc;
        printf("-------- TEST :: particle trajectory starts on element center and ends on same element insides :: ");
        p2.x = xc+1e-3; p2.y = yc+1e-3; p2.z = zc+1e-3;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem0) error_report(__FILE__,__LINE__,ielem0,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element center and ends on same element node :: ");
        p2.x = grid->node[grid->elem3d[ielem0].nodes[1]].x;
        p2.y = grid->node[grid->elem3d[ielem0].nodes[1]].y;
        p2.z = grid->node[grid->elem3d[ielem0].nodes[1]].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem0) error_report(__FILE__,__LINE__,ielem0,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element center and ends on same element edge :: ");
        p2.x = grid->node[grid->elem3d[ielem0].nodes[1]].x;
        p2.y = grid->node[grid->elem3d[ielem0].nodes[1]].y;
        p2.z = 0.5 * (grid->node[grid->elem3d[ielem0].nodes[1]].z + grid->node[grid->elem3d[ielem0].nodes[3]].z);
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem0) error_report(__FILE__,__LINE__,ielem0,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element center and ends on same element face :: ");
        p2.x = (1.0/3.0) * (grid->node[grid->elem3d[ielem0].nodes[0]].x + grid->node[grid->elem3d[ielem0].nodes[2]].x + grid->node[grid->elem3d[ielem0].nodes[3]].x);
        p2.y = (1.0/3.0) * (grid->node[grid->elem3d[ielem0].nodes[0]].y + grid->node[grid->elem3d[ielem0].nodes[2]].y + grid->node[grid->elem3d[ielem0].nodes[3]].y);
        p2.z = (1.0/3.0) * (grid->node[grid->elem3d[ielem0].nodes[0]].z + grid->node[grid->elem3d[ielem0].nodes[2]].z + grid->node[grid->elem3d[ielem0].nodes[3]].z);
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem0) error_report(__FILE__,__LINE__,ielem0,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element center and ends on next element center :: ");
        p2.x = (grid->node[grid->elem3d[ielem_adj].nodes[0]].x +
                grid->node[grid->elem3d[ielem_adj].nodes[1]].x +
                grid->node[grid->elem3d[ielem_adj].nodes[2]].x +
                grid->node[grid->elem3d[ielem_adj].nodes[3]].x) / 4.0;
        p2.y = (grid->node[grid->elem3d[ielem_adj].nodes[0]].y +
                grid->node[grid->elem3d[ielem_adj].nodes[1]].y +
                grid->node[grid->elem3d[ielem_adj].nodes[2]].y +
                grid->node[grid->elem3d[ielem_adj].nodes[3]].y) / 4.0;
        p2.z = (grid->node[grid->elem3d[ielem_adj].nodes[0]].z +
                grid->node[grid->elem3d[ielem_adj].nodes[1]].z +
                grid->node[grid->elem3d[ielem_adj].nodes[2]].z +
                grid->node[grid->elem3d[ielem_adj].nodes[3]].z) / 4.0;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem_adj) error_report(__FILE__,__LINE__,ielem_adj,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element center and ends on next element node :: ");
         p2.x = grid->node[grid->elem3d[ielem_adj].nodes[2]].x;
         p2.y = grid->node[grid->elem3d[ielem_adj].nodes[2]].y;
         p2.z = grid->node[grid->elem3d[ielem_adj].nodes[2]].z;
         //printf("node1: %d node2: %d\n",grid->elem3d[ielem0].nodes[1],grid->elem3d[ielem_adj].nodes[2]);
         ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
         if(ielem != ielem_adj) error_report(__FILE__,__LINE__,ielem_adj,ielem);
         printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element center and ends on next element edge :: ");
        p2.x = (grid->node[grid->elem3d[ielem_adj].nodes[2]].x + grid->node[grid->elem3d[ielem0].nodes[1]].x)/2.0;
        p2.y = (grid->node[grid->elem3d[ielem_adj].nodes[2]].y + grid->node[grid->elem3d[ielem0].nodes[1]].y)/2.0;
        p2.z = grid->node[grid->elem3d[ielem0].nodes[1]].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem_adj) error_report(__FILE__,__LINE__,ielem_adj,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element center and ends on next element face :: ");
        p2.x = (1.0/3.0) * (grid->node[grid->elem3d[ielem_adj].nodes[0]].x + grid->node[grid->elem3d[ielem_adj].nodes[2]].x + grid->node[grid->elem3d[ielem_adj].nodes[3]].x);
        p2.y = (1.0/3.0) * (grid->node[grid->elem3d[ielem_adj].nodes[0]].y + grid->node[grid->elem3d[ielem_adj].nodes[2]].y + grid->node[grid->elem3d[ielem_adj].nodes[3]].y);
        p2.z = (1.0/3.0) * (grid->node[grid->elem3d[ielem_adj].nodes[0]].z + grid->node[grid->elem3d[ielem_adj].nodes[2]].z + grid->node[grid->elem3d[ielem_adj].nodes[3]].z);
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem_adj) error_report(__FILE__,__LINE__,ielem_adj,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element center and ends on corner of grid :: ");
        p2.x = 0.0; p2.y = 0.0; p2.z = -100.0;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 1168 && ielem != 1169) error_report(__FILE__,__LINE__,1168,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element center and rides edges all the way to end of grid :: ");
        p2.x = 0.0; p2.y = 50.0; p2.z = -50.0;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 646) error_report(__FILE__,__LINE__,646,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element center and ends on external boundary face :: ");
        p2.x = 0.0; p2.y = 0.1; p2.z = -99.9;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 1169) error_report(__FILE__,__LINE__,1169,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        
        printf("-------- TEST :: particle trajectory starts on element center and rides grid to top-left grid corner :: ");
        node_id = 10;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element center and rides grid to bottom-left grid corner :: ");
        node_id = 230;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element center and rides grid to top-right grid corner :: ");
        node_id = 4630;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element center and rides grid to bottom-right grid corner :: ");
        node_id = 4850;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and leaves grid @ node :: ");
        node_id = 4630;
        p2.x = grid->node[node_id].x + 5.0;   p2.y = grid->node[node_id].y + 5.0;  p2.z = grid->node[node_id].z - 5.0;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and leaves grid @ edge :: ");
        node_id = 10;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z - .5;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and leaves grid @ face :: ");
        node_id = 10;
        p2.x = grid->node[node_id].x-1;   p2.y = grid->node[node_id].y-1;  p2.z = grid->node[node_id].z + .5;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);


        //----------------------------------------------------------------------------------------------
        printf("---- TESTING 3D ENGINE ELEMENT EDGE START\n");
        p1.x = grid->node[grid->elem3d[ielem0].nodes[1]].x;
        p1.y = grid->node[grid->elem3d[ielem0].nodes[1]].y;
        p1.z = 0.5 * (grid->node[grid->elem3d[ielem0].nodes[1]].z + grid->node[grid->elem3d[ielem0].nodes[3]].z);
        printf("-------- TEST :: particle trajectory starts on element edge and ends on same element insides :: ");
        p2.x = xc+1e-3; p2.y = yc+1e-3; p2.z = zc+1e-3;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem0) error_report(__FILE__,__LINE__,ielem0,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element edge and ends on same element node :: ");
        p2.x = grid->node[grid->elem3d[ielem0].nodes[1]].x;
        p2.y = grid->node[grid->elem3d[ielem0].nodes[1]].y;
        p2.z = grid->node[grid->elem3d[ielem0].nodes[1]].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem0) error_report(__FILE__,__LINE__,ielem0,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element edge and ends on same element edge :: ");
        p2.x = grid->node[grid->elem3d[ielem0].nodes[1]].x;
        p2.y = grid->node[grid->elem3d[ielem0].nodes[1]].y;
        p2.z = 0.5 * (grid->node[grid->elem3d[ielem0].nodes[1]].z + grid->node[grid->elem3d[ielem0].nodes[3]].z);
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem0) error_report(__FILE__,__LINE__,ielem0,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element edge and ends on same element face :: ");
        p2.x = (1.0/3.0) * (grid->node[grid->elem3d[ielem0].nodes[0]].x + grid->node[grid->elem3d[ielem0].nodes[2]].x + grid->node[grid->elem3d[ielem0].nodes[3]].x);
        p2.y = (1.0/3.0) * (grid->node[grid->elem3d[ielem0].nodes[0]].y + grid->node[grid->elem3d[ielem0].nodes[2]].y + grid->node[grid->elem3d[ielem0].nodes[3]].y);
        p2.z = (1.0/3.0) * (grid->node[grid->elem3d[ielem0].nodes[0]].z + grid->node[grid->elem3d[ielem0].nodes[2]].z + grid->node[grid->elem3d[ielem0].nodes[3]].z);
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem0) error_report(__FILE__,__LINE__,ielem0,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element edge and ends on next element center :: ");
        p2.x = (grid->node[grid->elem3d[ielem_adj].nodes[0]].x +
                grid->node[grid->elem3d[ielem_adj].nodes[1]].x +
                grid->node[grid->elem3d[ielem_adj].nodes[2]].x +
                grid->node[grid->elem3d[ielem_adj].nodes[3]].x) / 4.0;
        p2.y = (grid->node[grid->elem3d[ielem_adj].nodes[0]].y +
                grid->node[grid->elem3d[ielem_adj].nodes[1]].y +
                grid->node[grid->elem3d[ielem_adj].nodes[2]].y +
                grid->node[grid->elem3d[ielem_adj].nodes[3]].y) / 4.0;
        p2.z = (grid->node[grid->elem3d[ielem_adj].nodes[0]].z +
                grid->node[grid->elem3d[ielem_adj].nodes[1]].z +
                grid->node[grid->elem3d[ielem_adj].nodes[2]].z +
                grid->node[grid->elem3d[ielem_adj].nodes[3]].z) / 4.0;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem_adj) error_report(__FILE__,__LINE__,ielem_adj,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element edge and ends on next element node :: ");
        p2.x = grid->node[grid->elem3d[ielem_adj].nodes[2]].x;
        p2.y = grid->node[grid->elem3d[ielem_adj].nodes[2]].y;
        p2.z = grid->node[grid->elem3d[ielem_adj].nodes[2]].z;
        //printf("end node %d\n",grid->elem3d[ielem_adj].nodes[2]);
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if (grid->elem3d[ielem].nodes[0]!=2183 && grid->elem3d[ielem].nodes[1]!=2183 && grid->elem3d[ielem].nodes[2]!=2183 && grid->elem3d[ielem].nodes[3]!=2183) error_report(__FILE__,__LINE__,11354,ielem); // can be a number of elements
        //if(ielem != 11354) error_report(__FILE__,__LINE__,11354,ielem); // can be a number of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element edge and ends on next element edge :: ");
        p2.x = (grid->node[grid->elem3d[ielem_adj].nodes[2]].x + grid->node[grid->elem3d[ielem0].nodes[1]].x)/2.0;
        p2.y = (grid->node[grid->elem3d[ielem_adj].nodes[2]].y + grid->node[grid->elem3d[ielem0].nodes[1]].y)/2.0;
        p2.z = grid->node[grid->elem3d[ielem0].nodes[1]].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 11354) error_report(__FILE__,__LINE__,11354,ielem); // can be a number of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element edge and ends on next element face :: ");
        p2.x = (1.0/3.0) * (grid->node[grid->elem3d[ielem_adj].nodes[0]].x + grid->node[grid->elem3d[ielem_adj].nodes[2]].x + grid->node[grid->elem3d[ielem_adj].nodes[3]].x);
        p2.y = (1.0/3.0) * (grid->node[grid->elem3d[ielem_adj].nodes[0]].y + grid->node[grid->elem3d[ielem_adj].nodes[2]].y + grid->node[grid->elem3d[ielem_adj].nodes[3]].y);
        p2.z = (1.0/3.0) * (grid->node[grid->elem3d[ielem_adj].nodes[0]].z + grid->node[grid->elem3d[ielem_adj].nodes[2]].z + grid->node[grid->elem3d[ielem_adj].nodes[3]].z);
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem_adj) error_report(__FILE__,__LINE__,ielem_adj,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element edge and ends on corner of grid :: ");
        p2.x = 0.0; p2.y = 0.0; p2.z = -100.0;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 1168 && ielem != 1169) error_report(__FILE__,__LINE__,1168,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element edge and rides edges all the way to end of grid :: ");
        p2.x = 0.0; p2.y = 50.0; p2.z = -50.0; //
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != 115 &&
           grid->elem3d[ielem].nodes[1] != 115 &&
           grid->elem3d[ielem].nodes[2] != 115 &&
           grid->elem3d[ielem].nodes[3] != 115) error_report(__FILE__,__LINE__,UNSET_INT,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element edge and ends on external boundary face :: ");
        p2.x = 0.0; p2.y = 0.1; p2.z = -99.9;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 1169) error_report(__FILE__,__LINE__,1169,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        
        printf("-------- TEST :: particle trajectory starts on edge and rides grid to top-left grid corner :: ");
        node_id = 10;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on edge and rides grid to bottom-left grid corner :: ");
        node_id = 230;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on edge and rides grid to top-right grid corner :: ");
        node_id = 4630;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on edge and rides grid to bottom-right grid corner :: ");
        node_id = 4850;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and leaves grid @ node :: ");
        node_id = 4630;
        p2.x = grid->node[node_id].x + 5.0;   p2.y = grid->node[node_id].y + 5.0;  p2.z = grid->node[node_id].z - 5.0;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and leaves grid @ edge :: ");
        node_id = 10;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z - .5;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and leaves grid @ face :: ");
        node_id = 10;
        p2.x = grid->node[node_id].x-1;   p2.y = grid->node[node_id].y-1;  p2.z = grid->node[node_id].z + .5;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        
        //----------------------------------------------------------------------------------------------
        printf("---- TESTING 3D ENGINE ELEMENT FACE START\n");
        p1.x = (1.0/3.0) * (grid->node[grid->elem3d[ielem0].nodes[0]].x + grid->node[grid->elem3d[ielem0].nodes[2]].x + grid->node[grid->elem3d[ielem0].nodes[3]].x);
        p1.y = (1.0/3.0) * (grid->node[grid->elem3d[ielem0].nodes[0]].y + grid->node[grid->elem3d[ielem0].nodes[2]].y + grid->node[grid->elem3d[ielem0].nodes[3]].y);
        p1.z = (1.0/3.0) * (grid->node[grid->elem3d[ielem0].nodes[0]].z + grid->node[grid->elem3d[ielem0].nodes[2]].z + grid->node[grid->elem3d[ielem0].nodes[3]].z);
        printf("-------- TEST :: particle trajectory starts on element face and ends on same element insides :: ");
        p2.x = xc+1e-3; p2.y = yc+1e-3; p2.z = zc+1e-3;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem0) error_report(__FILE__,__LINE__,ielem0,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element face and ends on same element node :: ");
        p2.x = grid->node[grid->elem3d[ielem0].nodes[1]].x;
        p2.y = grid->node[grid->elem3d[ielem0].nodes[1]].y;
        p2.z = grid->node[grid->elem3d[ielem0].nodes[1]].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem0) error_report(__FILE__,__LINE__,ielem0,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element face and ends on same element edge :: ");
        p2.x = grid->node[grid->elem3d[ielem0].nodes[1]].x;
        p2.y = grid->node[grid->elem3d[ielem0].nodes[1]].y;
        p2.z = 0.5 * (grid->node[grid->elem3d[ielem0].nodes[1]].z + grid->node[grid->elem3d[ielem0].nodes[3]].z);
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem0) error_report(__FILE__,__LINE__,ielem0,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element face and ends on same element face :: ");
        p2.x = (1.0/3.0) * (grid->node[grid->elem3d[ielem0].nodes[0]].x + grid->node[grid->elem3d[ielem0].nodes[2]].x + grid->node[grid->elem3d[ielem0].nodes[3]].x);
        p2.y = (1.0/3.0) * (grid->node[grid->elem3d[ielem0].nodes[0]].y + grid->node[grid->elem3d[ielem0].nodes[2]].y + grid->node[grid->elem3d[ielem0].nodes[3]].y);
        p2.z = (1.0/3.0) * (grid->node[grid->elem3d[ielem0].nodes[0]].z + grid->node[grid->elem3d[ielem0].nodes[2]].z + grid->node[grid->elem3d[ielem0].nodes[3]].z);
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem0) error_report(__FILE__,__LINE__,ielem0,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element face and ends on next element center :: ");
        p2.x = (grid->node[grid->elem3d[ielem_adj].nodes[0]].x +
                grid->node[grid->elem3d[ielem_adj].nodes[1]].x +
                grid->node[grid->elem3d[ielem_adj].nodes[2]].x +
                grid->node[grid->elem3d[ielem_adj].nodes[3]].x) / 4.0;
        p2.y = (grid->node[grid->elem3d[ielem_adj].nodes[0]].y +
                grid->node[grid->elem3d[ielem_adj].nodes[1]].y +
                grid->node[grid->elem3d[ielem_adj].nodes[2]].y +
                grid->node[grid->elem3d[ielem_adj].nodes[3]].y) / 4.0;
        p2.z = (grid->node[grid->elem3d[ielem_adj].nodes[0]].z +
                grid->node[grid->elem3d[ielem_adj].nodes[1]].z +
                grid->node[grid->elem3d[ielem_adj].nodes[2]].z +
                grid->node[grid->elem3d[ielem_adj].nodes[3]].z) / 4.0;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem_adj) error_report(__FILE__,__LINE__,ielem_adj,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element face and ends on next element node :: ");
        p2.x = grid->node[grid->elem3d[ielem_adj].nodes[2]].x;
        p2.y = grid->node[grid->elem3d[ielem_adj].nodes[2]].y;
        p2.z = grid->node[grid->elem3d[ielem_adj].nodes[2]].z;
        //printf("node1: %d node2: %d\n",grid->elem3d[ielem0].nodes[1],grid->elem3d[ielem_adj].nodes[2]);
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != grid->elem3d[ielem_adj].nodes[2] &&
           grid->elem3d[ielem].nodes[1] != grid->elem3d[ielem_adj].nodes[2] &&
           grid->elem3d[ielem].nodes[2] != grid->elem3d[ielem_adj].nodes[2] &&
           grid->elem3d[ielem].nodes[3] != grid->elem3d[ielem_adj].nodes[2]) error_report(__FILE__,__LINE__,UNSET_INT,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element face and ends on next element edge :: ");
        p2.x = (grid->node[grid->elem3d[ielem_adj].nodes[2]].x + grid->node[grid->elem3d[ielem0].nodes[1]].x)/2.0;
        p2.y = (grid->node[grid->elem3d[ielem_adj].nodes[2]].y + grid->node[grid->elem3d[ielem0].nodes[1]].y)/2.0;
        p2.z = grid->node[grid->elem3d[ielem0].nodes[1]].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 11355) error_report(__FILE__,__LINE__,11355,ielem); // can be a number of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element face and ends on next element face :: ");
        p2.x = (1.0/3.0) * (grid->node[grid->elem3d[ielem_adj].nodes[0]].x + grid->node[grid->elem3d[ielem_adj].nodes[2]].x + grid->node[grid->elem3d[ielem_adj].nodes[3]].x);
        p2.y = (1.0/3.0) * (grid->node[grid->elem3d[ielem_adj].nodes[0]].y + grid->node[grid->elem3d[ielem_adj].nodes[2]].y + grid->node[grid->elem3d[ielem_adj].nodes[3]].y);
        p2.z = (1.0/3.0) * (grid->node[grid->elem3d[ielem_adj].nodes[0]].z + grid->node[grid->elem3d[ielem_adj].nodes[2]].z + grid->node[grid->elem3d[ielem_adj].nodes[3]].z);
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != ielem_adj) error_report(__FILE__,__LINE__,ielem_adj,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element face and ends on corner of grid :: ");
        p2.x = 0.0; p2.y = 0.0; p2.z = -100.0;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 1168 && ielem != 1169) error_report(__FILE__,__LINE__,1168,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element face and rides edges all the way to end of grid :: ");
        p2.x = 0.0; p2.y = 50.0; p2.z = -50.0; //
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != 115 &&
           grid->elem3d[ielem].nodes[1] != 115 &&
           grid->elem3d[ielem].nodes[2] != 115 &&
           grid->elem3d[ielem].nodes[3] != 115) error_report(__FILE__,__LINE__,UNSET_INT,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on element face and ends on external boundary face :: ");
        p2.x = 0.0; p2.y = 0.1; p2.z = -99.9;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(ielem != 1169) error_report(__FILE__,__LINE__,1169,ielem);
        printf("IELEM: %d || PASSED\n",ielem);
        
        printf("-------- TEST :: particle trajectory starts on face and rides grid to top-left grid corner :: ");
        node_id = 10;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on face and rides grid to bottom-left grid corner :: ");
        node_id = 230;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on face and rides grid to top-right grid corner :: ");
        node_id = 4630;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on face and rides grid to bottom-right grid corner :: ");
        node_id = 4850;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and leaves grid @ node :: ");
        node_id = 4630;
        p2.x = grid->node[node_id].x + 5.0;   p2.y = grid->node[node_id].y + 5.0;  p2.z = grid->node[node_id].z - 5.0;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and leaves grid @ edge :: ");
        node_id = 10;
        p2.x = grid->node[node_id].x;   p2.y = grid->node[node_id].y;  p2.z = grid->node[node_id].z - .5;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        printf("-------- TEST :: particle trajectory starts on node and leaves grid @ face :: ");
        node_id = 10;
        p2.x = grid->node[node_id].x-1;   p2.y = grid->node[node_id].y-1;  p2.z = grid->node[node_id].z + .5;
        ielem = ielem0; flag = elementSearch3D(grid,p1,p2,&ielem,true,&pi);
        if(grid->elem3d[ielem].nodes[0] != node_id &&
           grid->elem3d[ielem].nodes[1] != node_id &&
           grid->elem3d[ielem].nodes[2] != node_id &&
           grid->elem3d[ielem].nodes[3] != node_id) error_report(__FILE__,__LINE__,node_id,ielem); // can be a # of elements
        printf("IELEM: %d || PASSED\n",ielem);
        
        
        // TEST SHAKE AND BAKE
        //    int ielem = 38, lnode = 1;
        //    SVECT pp = shakeAndBake(grid,ielem,lnode);
        //    printf("node: %20.10e %20.10e || p: %20.10e %20.10e\n",
        //           grid->node[ grid->elem2d[ielem].nodes[lnode] ].x,
        //           grid->node[ grid->elem2d[ielem].nodes[lnode] ].y,
        //           pp.x,pp.y);
        //    exit(-1);
        
        //create_velocity_file_test(grid,t0,tf,50,file_base,getV_2d_uniform_oscillating);
    }
    
    return 1;
}

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

double run_verification(SGRID *grid, char *file_base, int np, int *node_ID, int *elem_ID, double dt, double tf,
                     SVECT (*get_v)(SVECT p,  double t),
                     SVECT (*get_p)(SVECT p0, double t),
                     SVECT (*get_wq)(SVECT p0, double t),
                     int vel_spool, int behavior_flag, int spool,
                     double error_tolerance, int ptype) {
    int ip, snaps = ceil((tf)/(double)vel_spool);
    double t0 = 0.;
    char vel_file[100],pt_file[100],error_file[100];

    
    if (RUN_VERBOSE) printf("------ creating test velocity file\n");
    create_velocity_file_test(grid,t0,tf,snaps,file_base,get_v);
    
    if (get_wq != NULL) {
        if (RUN_VERBOSE) printf("------ creating test water quality file\n");
        create_water_quality_file_test(grid,t0,tf,snaps,file_base,get_wq);
    }
    
    if (grid->vertical_dpl_flag == ON) {
        if (RUN_VERBOSE) printf("------ creating test nodal displacement file\n");
        create_nd_dpl_file_test(grid,t0,tf,snaps,file_base);
    }
    
    if (RUN_VERBOSE) printf("------ allocated and initialize test model\n");
    SMODEL *mod = (SMODEL *) malloc(sizeof(SMODEL));
    strcpy(mod->file_base,file_base);
    mod->get_analytic_p = get_p;
    mod->get_analytic_v = get_v;
    mod->get_analytic_wq = get_wq;
    
    if (RUN_VERBOSE) printf("------ initialize test model\n");
    smodel_init_noRead(mod,grid,file_base,t0,tf,dt,np,node_ID,elem_ID,spool,ptype);
    for (ip=0; ip<mod->np; ip++) mod->p[ip].behavoiral_vel_flag = behavior_flag;
    mod->screen_print_timestep = false;
    
    if (RUN_VERBOSE) printf("------ running the test\n");
    smodel_run_lockstep(mod);
    double max_t_abs_error_mag = mod->max_t_abs_error_mag;
    
    if (RUN_VERBOSE) printf("------ finalizing the test model\n");
    smodel_finalize(mod);
    
    // now test error
    if (max_t_abs_error_mag > error_tolerance + 1e-6) {
        printf("MAX TIME ABS ERROR: %f || ERROR TOLERANCE %f || FAILED\n",max_t_abs_error_mag,error_tolerance);
        exit(-1);
    } else {
        printf("MAX TIME ABS ERROR: %f || ERROR TOLERANCE %f || PASSED\n",max_t_abs_error_mag,error_tolerance);
    }
    
    return max_t_abs_error_mag;
}

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

int pt_verification() {
    
    RUN_VERBOSE = false;
    
    int i,j,k;
    double abs_error_tmax = 0.;
    
    printf("\n PARTICLE TRAJECTORY VERIFICATION\n");
    char file_base[100], vel_file[100],pt_file[100],error_file[100],pt_file_hdf5[100];
    
    // -------------------------------------------------------------
    // build 2D test grid ------------------------------------------
    //printf("-- CREATING 2D VERIFICATION GRID\n");
    SGRID *grid;
    //sgrid_create2d(grid,0,100,0,100,21,21,0,200,100,0,0);
    
    printf("---- 2D TEST :: PARTICLE DIAGONAL OSCILLATION ON GRID || "); fflush(stdout);
    {
        //RUN_VERBOSE = true;
        grid = (SGRID *)malloc(sizeof(SGRID)); sgrid_read_adh(&grid,"test_2d");
        strcpy(file_base, "test_2d_oscillation");
        int node_ID[5], elem_ID[5];
        node_ID[0] = 220; elem_ID[0] = 381;
        node_ID[1] = 200; elem_ID[1] = 380;
        node_ID[2] = 180; elem_ID[2] = 342;
        node_ID[3] = 160; elem_ID[3] = 304;
        node_ID[4] = 140; elem_ID[4] = 266;
        abs_error_tmax = run_verification(grid,file_base,5/*np*/,node_ID,elem_ID,1/*dt*/,100./*tf*/,
                         get_vel_2d_uniform_oscillating,
                         get_analytic_p_2d_uniform_oscillating,
                         NULL,2,-1,1,0.19,GENERAL);
    }
    
    printf("---- 2D TEST :: 225 PARTICLE ROTATION || "); fflush(stdout);
    {
        //RUN_VERBOSE = true;
        grid = (SGRID *)malloc(sizeof(SGRID)); sgrid_read_adh(&grid,"test_2d");
        strcpy(file_base, "test_2d_rotation_169");
        int node_ID[169], elem_ID[169], k=0;
        for (i=0; i<grid->nnodes; i++) {
            if (grid->node[i].x > 15 && grid->node[i].x < 85 && grid->node[i].y > 15 && grid->node[i].y < 85 ) {
                node_ID[k] = i;
                for (j=0; j<grid->nelems2d; j++) {
                    if (grid->elem2d[j].nodes[0] == i || grid->elem2d[j].nodes[1] == i || grid->elem2d[j].nodes[2] == i) {
                        elem_ID[k] = j;
                        break;
                    }
                }
                k++;
            }
        }
        abs_error_tmax = run_verification(grid,file_base,169/*np*/,node_ID,elem_ID,1/*dt*/,100./*tf*/,
                         get_vel_2d_rotation,
                         get_analytic_p_2d_rotation,
                         NULL,2,-1,1,0.1754029,GENERAL);
    }

    printf("---- 2D TEST :: CLOSED BOUNDARY TEST || "); fflush(stdout);
    {
        // NOTE :: not sure how to test this for error, it's more of a qualitative test
        //RUN_VERBOSE = true;
        grid = (SGRID *)malloc(sizeof(SGRID)); sgrid_read_adh(&grid,"test_2d");
        strcpy(file_base, "test_2d_closed_boundary");
        int node_ID[4], elem_ID[4];
        node_ID[0] = 430; elem_ID[0] = 779; // (100,50)
        node_ID[1] = 420; elem_ID[1] = 761; // (100,100)
        node_ID[2] = 0;   elem_ID[2] = 0;   // (0,0)
        node_ID[3] = 220; elem_ID[3] = 418; // (50,50)
        abs_error_tmax = run_verification(grid,file_base,4/*np*/,node_ID,elem_ID,1/*dt*/,300./*tf*/,
                         get_vel_2d_closed_boundary,
                         get_analytic_p_2d_closed_boundary,
                         NULL,1,-1,1,1e-6,GENERAL);
    }
    
//    printf("---- 2D TEST :: PARTICLE DIAGONAL OSCILLATION OFF GRID *OPEN* BOUNDARY || ");
//    {
//        //RUN_VERBOSE = true;
//        grid = (SGRID *)malloc(sizeof(SGRID)); sgrid_read_adh(&grid,"test_2d");
//        strcpy(file_base, "test_2d_oscillation_off_grid");
//        int node_ID = 220, elem_ID = 379;
//        run_verification(grid,file_base,1/*np*/,&node_ID,&elem_ID,1/*dt*/,100./*tf*/,
//                         get_vel_2d_uniform_oscillating_off_grid,
//                         get_analytic_p_2d_uniform_oscillating_off_grid,
//                         &abs_error_tmax,0.104945);
//    }
    

    
    
    // -------------------------------------------------------------
    // build 3D test grid ------------------------------------------
    //printf("-- READING 3D VERIFICATION GRID\n");
//    free(grid);
//    grid = (SGRID *)malloc(sizeof(SGRID));
//    sgrid_read_adh(&grid,"test_3d");
    
    printf("---- 3D TEST :: PARTICLE SPIRAL (TRANSIENT VELOCITY) || ");
    {
        //RUN_VERBOSE = true;
        grid = (SGRID *)malloc(sizeof(SGRID)); sgrid_read_adh(&grid,"test_3d");
        strcpy(file_base, "test_3d_rotation_507");
        int node_ID[507], elem_ID[507], k=0;
        for (i=0; i<grid->nnodes; i++) {
            if (grid->node[i].x > 15 && grid->node[i].x < 85 && grid->node[i].y > 15 && grid->node[i].y < 85 && grid->node[i].z > -65 && grid->node[i].z < -35) {
                node_ID[k] = i;
                for (j=0; j<grid->nelems3d; j++) {
                    if (grid->elem3d[j].nodes[0] == i || grid->elem3d[j].nodes[1] == i || grid->elem3d[j].nodes[2] == i || grid->elem3d[j].nodes[3] == i) {
                        elem_ID[k] = j;
                        break;
                    }
                }
                k++;
                //printf("k: %d\n",k); fflush(stdout);
            }
        }
        assert(k == 507);
        //exit(-1);
        abs_error_tmax = run_verification(grid,file_base,507/*np*/,node_ID,elem_ID,1/*dt*/,100./*tf*/,
                         get_vel_3d_rotation,
                         get_analytic_p_3d_rotation,
                         NULL,2,-1,1,0.1754469,GENERAL);
    }
    
//    printf("---- 3D TEST :: PARTICLE OFF GRID NODE (LINEAR) ** OPEN BOUNDARY ** || ");
//    {
//        RUN_VERBOSE = true;
//        grid = (SGRID *)malloc(sizeof(SGRID)); sgrid_read_adh(&grid,"test_3d");
//        strcpy(file_base, "test_3d_off_grid_node");
//        int node_ID = 3580, elem_ID = 17352;
//        run_verification(grid,file_base,1/*np*/,&node_ID,&elem_ID,1/*dt*/,100./*tf*/,
//                         get_vel_3d_off_grid_node,
//                         get_analytic_p_3d_off_grid_node,
//                         &abs_error_tmax,1e-10);
//    }
    printf("---- 3D TEST :: PARTICLE OFF GRID EDGE (LINEAR) || ");
    {
        //RUN_VERBOSE = true;
        grid = (SGRID *)malloc(sizeof(SGRID)); sgrid_read_adh(&grid,"test_3d");
        strcpy(file_base, "test_3d_off_grid_edge");
        int node_ID = 3580, elem_ID = 17352;
        abs_error_tmax = run_verification(grid,file_base,1/*np*/,&node_ID,&elem_ID,1/*dt*/,100./*tf*/,
                         get_vel_3d_off_grid_edge,
                         get_analytic_p_3d_off_grid_edge,
                         NULL,2,-1,1,10000,GENERAL);
    }
//    printf("---- 3D TEST :: PARTICLE OFF GRID FACE (OPEN BOUNDARY - LINEAR) || ");
//    {
          //RUN_VERBOSE = true;
//        strcpy(file_base, "test_3d_off_grid_face");
//        int node_ID = 3580, elem_ID = 17352;
//        run_verification(grid,file_base,1/*np*/,&node_ID,&elem_ID,1/*dt*/,100./*tf*/,
//                         get_vel_3d_off_grid_face,
//                         get_analytic_p_3d_off_grid_face,
//                         &abs_error_tmax,1e-10);
//    }
    
    //-----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------
    // OYSTER TEST
    printf("---- 3D TEST :: OYSTER FALL VELOCITY || ");
    {
        //RUN_VERBOSE = true;
        grid = (SGRID *)malloc(sizeof(SGRID)); sgrid_read_adh(&grid,"test_3d");
        strcpy(file_base, "test_3d_oyster_fall");
        int node_ID = 2420, elem_ID = 11430;
        abs_error_tmax = run_verification(grid,file_base,1/*np*/,&node_ID,&elem_ID,100/*dt*/,100000./*tf*/,
                                          get_vel_3d_oyster_fall,
                                          get_analytic_p_3d_oyster_fall,
                                          NULL,1000,0,10,0.022487,OYSTER);
    }
    printf("---- 3D TEST :: OYSTER SWIM VELOCITY || ");
    {
        //RUN_VERBOSE = true;
        grid = (SGRID *)malloc(sizeof(SGRID)); sgrid_read_adh(&grid,"test_3d");
        strcpy(file_base, "test_3d_oyster_swim");
        int node_ID = 2430, elem_ID = 11459;
        abs_error_tmax = run_verification(grid,file_base,1/*np*/,&node_ID,&elem_ID,100/*dt*/,100000./*tf*/,
                                          get_vel_3d_oyster_swim,
                                          get_analytic_p_3d_oyster_swim,
                                          NULL,1000,1,10,0.006255,OYSTER);
    }
    printf("---- 3D TEST :: OYSTER TOTAL VELOCITY || ");
    {
        //RUN_VERBOSE = true;
        grid = (SGRID *)malloc(sizeof(SGRID)); sgrid_read_adh(&grid,"test_3d");
        strcpy(file_base, "test_3d_oyster_tvel");
        int node_ID = 2430, elem_ID = 11459;
        abs_error_tmax = run_verification(grid,file_base,1/*np*/,&node_ID,&elem_ID,100/*dt*/,100000./*tf*/,
                                          get_vel_3d_oyster_tvel,
                                          get_analytic_p_3d_oyster_tvel,
                                          NULL,1000,2,10,0.004094,OYSTER);
    }
    printf("---- 3D TEST :: OYSTER OXYGEN TEST || ");
    { // oyster starts on bed & swims/falls until it gets into a oxygen rich depth, then it dies (inactive) and only falls
        //RUN_VERBOSE = true;
        grid = (SGRID *)malloc(sizeof(SGRID)); sgrid_read_adh(&grid,"test_3d");
        strcpy(file_base, "test_3d_oyster_oxygen");
        int node_ID = 2430, elem_ID = 11459;
        abs_error_tmax = run_verification(grid,file_base,1/*np*/,&node_ID,&elem_ID,100/*dt*/,680000./*tf*/,
                                          get_vel_3d_oyster_tvel,
                                          soyster_get_analytic_oxygen,
                                          get_wq_3d_oyster_oxygen,
                                          4000,2,10,0.02,OYSTER);
    }
    printf("---- 3D TEST :: OYSTER SALNITY TEST || ");
    { // oyster starts on bed & swims/falls until it gets into a oxygen rich depth, then it dies (inactive) and only falls
        //RUN_VERBOSE = true;
        grid = (SGRID *)malloc(sizeof(SGRID)); sgrid_read_adh(&grid,"test_3d");
        strcpy(file_base, "test_3d_oyster_salinity");
        int node_ID = 2430, elem_ID = 11459;
        abs_error_tmax = run_verification(grid,file_base,1/*np*/,&node_ID,&elem_ID,100/*dt*/,680000./*tf*/,
                                          get_vel_3d_oyster_tvel,
                                          soyster_get_analytic_oxygen,
                                          get_wq_3d_oyster_salinity,
                                          4000,2,10,0.02,OYSTER);
    }
    printf("---- 3D TEST :: OYSTER SUNLIGHT TEST || ");
    { // oyster starts on bed & swims/falls until it gets into a oxygen rich depth, then it dies (inactive) and only falls
        //RUN_VERBOSE = true;
        grid = (SGRID *)malloc(sizeof(SGRID)); sgrid_read_adh(&grid,"test_3d");
        strcpy(file_base, "test_3d_oyster_sunlight");
        int node_ID = 2430, elem_ID = 11459;
        abs_error_tmax = run_verification(grid,file_base,1/*np*/,&node_ID,&elem_ID,100/*dt*/,680000./*tf*/,
                                          get_vel_3d_oyster_tvel,
                                          soyster_get_analytic_oxygen,
                                          get_wq_3d_oyster_sunlight,
                                          4000,2,10,0.02,OYSTER);
    }
    
    //-----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------
    printf("---- 3D TEST :: SURFACE DISPLACEMENT || ");
    { // particle rides along a tidal displaced surface
        //RUN_VERBOSE = true;
        grid = (SGRID *)malloc(sizeof(SGRID)); sgrid_read_adh(&grid,"test_3d");
        grid->vertical_dpl_flag = ON;
        strcpy(file_base, "test_3d_surface_dpl");
        int node_ID = 2420, elem_ID = 12540;
        abs_error_tmax = run_verification(grid,file_base,1/*np*/,&node_ID,&elem_ID,1/*dt*/,600./*tf*/,
                                          get_vel_3d_surface_dpl,
                                          get_analytic_p_3d_surface_dpl,
                                          NULL,2,-1,1,0.047519,GENERAL);
        grid->vertical_dpl_flag = OFF;
    }

    return 1;
}

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

void create_velocity_file_test(SGRID *grid, double t1, double t2, int nt, char *file_base, SVECT (*get_vel)(SVECT p, double t)) {
    int i;
    double t = t1,dt;
    
    SVECT v;
    FILE *fp_vel;
    char vel_file[100];
    strcpy(vel_file,file_base);
    strcat(vel_file, "_vel.dat");
    fp_vel = fopen(vel_file,"w");
    if (RUN_VERBOSE) printf("-- CREATING VELOCITY FILE: %s\n",vel_file);
    fprintf(fp_vel,"dum\ndum\ndum\ndum\ndum\ndum\ndum\n");
    fprintf(fp_vel, "TIMEUNITS SECONDS\n");
    
    dt = (t2 - t1) / ((double) nt);
    while(1) {
        //printf("dt: %f t: %f\n",dt,t);
        if (t > t2) break;
        fprintf(fp_vel,"TS 0 %f\n",t);
        
        for (i=0; i<grid->nnodes; i++) {
            v = get_vel(grid->node[i],t);
            fprintf(fp_vel,"%20.10e %20.10e %20.10e\n",v.x,v.y,v.z);
        }
        t += dt;
    }
    fprintf(fp_vel,"ENDDS\n");
    fclose(fp_vel);
}

void create_water_quality_file_test(SGRID *grid, double t1, double t2, int nt, char *file_base, SVECT (*get_wq)(SVECT p, double t)) {
    int i;
    double t = t1,dt;
    
    SVECT v; // stores x - oxygen, y - salinity, z - sunlight
    FILE *fp;
    char file[100];
    strcpy(file,file_base);
    strcat(file, "_wq.dat");
    fp = fopen(file,"w");
    if (RUN_VERBOSE) printf("-- CREATING WATER QUALITY FILE: %s\n",file);
    fprintf(fp,"dum\ndum\ndum\ndum\ndum\ndum\ndum\n");
    
    dt = (t2 - t1) / ((double) nt);
    while(1) {
        //printf("dt: %f t: %f\n",dt,t);
        if (t > t2) break;
        fprintf(fp,"TS 0 %f\n",t);
        
        for (i=0; i<grid->nnodes; i++) {
            v = get_wq(grid->node[i],t);
            fprintf(fp,"%20.10e %20.10e %20.10e\n",v.x,v.y,v.z);
        }
        t += dt;
    }
    fprintf(fp,"ENDDS\n");
    fclose(fp);
}

 
void create_nd_dpl_file_test(SGRID *grid, double t1, double t2, int nt, char *file_base) {
    int i;
    double t = t1,dt;
    
    double dpl,dpl_surface;
    
    FILE *fp;
    char file[100];
    strcpy(file,file_base);
    strcat(file, "_dpl.dat");
    fp = fopen(file,"w");
    if (RUN_VERBOSE) printf("-- CREATING GRID DISPLACEMENT FILE: %s\n",file);
    fprintf(fp,"dum\ndum\ndum\ndum\ndum\ndum\ndum\n");
    
    dt = (t2 - t1) / ((double) nt);
    while(1) {
        //printf("dt: %f t: %f\n",dt,t);
        if (t > t2) break;
        fprintf(fp,"TS 0 %f\n",t);
        
        dpl_surface = 3*sin(1.57079633*t/150.);
        double max_dpl = -1e6;
        for (i=0; i<grid->nnodes; i++) {
            dpl = (grid->node[i].z + 100.0) * 0.01 * dpl_surface; // AdH accordian style displacement of 1m/s
            fprintf(fp,"%20.10e\n",dpl);
            //printf("node: %d \t dpl: %f\n",i,dpl);
            if (fabs(dpl) > max_dpl) max_dpl = fabs(dpl);
        }
        //printf("create_nd_dpl_file_test :: time: %f \t max_abs(dpl): %f :: dpl_surface: %f\n",t,max_dpl,dpl_surface);
        t += dt;
    }
    fprintf(fp,"ENDDS\n");
    fclose(fp);
}

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

SVECT get_errors_test(SVECT p, SVECT p0, double t, SVECT (*get_analytic_p)(SVECT p0, double t), SVECT *abs_error, SVECT *rel_error) {
    SVECT pa = get_analytic_p(p0,t);
    
    (*abs_error).x = fabs(p.x - pa.x);
    (*abs_error).y = fabs(p.y - pa.y);
    (*abs_error).z = fabs(p.z - pa.z);
    
    if (fabs(pa.x) > 1e-6) {
        (*rel_error).x = 100*fabs(p.x - pa.x)/fabs(pa.x);
    } else {
        (*rel_error).x = (*abs_error).x;
    }
    if (fabs(pa.y) > 1e-6) {
        (*rel_error).y = 100*fabs(p.y - pa.y)/fabs(pa.y);
    } else {
        (*rel_error).y = (*abs_error).y;
    }
    if (fabs(pa.z) > 1e-6) {
        (*rel_error).z = 100*fabs(p.z - pa.z)/fabs(pa.z);
    } else {
        (*rel_error).z = (*abs_error).z;
    }
    
    if (RUN_VERBOSE)  printf("abs_error: {%f,%f,%f}  || rel_error: {%f%%,%f%%,%f%%}\n",abs_error->x,abs_error->y,abs_error->z,rel_error->x,rel_error->y,rel_error->z);
    
    return pa;
}

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
// 2D :: spatially constant grid velocity with vy=vy, linearly oscillating time
SVECT get_vel_2d_uniform_oscillating(SVECT p, double t) {
    SVECT v;
    v.x = 2 * sin(t/10.0);
    v.y = 2 * sin(t/10.0);
    v.z = 0.0;
    return v;
}
SVECT get_vel_2d_uniform_oscillating_off_grid(SVECT p, double t) {
    SVECT v;
    v.x = 3 * sin(t/10.0);
    v.y = 3 * sin(t/10.0);
    v.z = 0.0;
    return v;
}
SVECT get_analytic_p_2d_uniform_oscillating(SVECT p0, double t) {
    SVECT p;
    p.x = p0.x + 20 * (-cos(t/10.) + 1);
    p.y = p0.y + 20 * (-cos(t/10.) + 1);
    p.z = p0.z;
    return p;
}
SVECT get_analytic_p_2d_uniform_oscillating_off_grid(SVECT p0, double t) {
    SVECT p;
    if (t<23) {
        p.x = p0.x + 30 * (-cos(t/10.) + 1);
        p.y = p0.y + 30 * (-cos(t/10.) + 1);
        p.z = p0.z;
    } else {
        p.x = 100.0; p.y = 100.0; p.z = -100.0;
    }
    return p;
}
SVECT get_vel_2d_rotation(SVECT p, double t) {
    SVECT v;
    double r = pow((p.x-50.0) * (p.x-50.0) + (p.y-50.) * (p.y-50.),0.5);
    double theta = atan2 ((p.y-50.),(p.x-50.0));
    double w = PI / 50.0; // radian/second -- time for one orbit, 2 * pi / w = 100s
    // note :: w = d(theta)/dt
    v.x =  -w * r * sin(theta);
    v.y =   w * r * cos(theta);
    v.z = 0.0;
    return v;
}
SVECT get_analytic_p_2d_rotation(SVECT p0, double t) {
    double r = pow((p0.x-50.0) * (p0.x-50.0) + (p0.y-50.) * (p0.y-50.),0.5);
    double w = PI / 50.0; // radian/second -- time for one orbit, 2 * pi / w = 100s
    SVECT p;
    p.x =  (p0.x - 50.) * cos(w*t) - (p0.y - 50.) * sin(w*t) + 50.0;
    p.y =  (p0.x - 50.) * sin(w*t) + (p0.y - 50.) * cos(w*t) + 50.0;
    p.z = p0.z;
    return p;
}
SVECT get_vel_3d_rotation(SVECT p, double t) {
    SVECT v;
    double r = pow((p.x-50.0) * (p.x-50.0) + (p.y-50.) * (p.y-50.),0.5);
    double theta = atan2 ((p.y-50.),(p.x-50.0));
    double w = PI / 50.0; // radian/second -- time for one orbit, 2 * pi / w = 100s
    v.x =  -w * r * sin(theta);
    v.y =   w * r * cos(theta);
    v.z =   sin(2 * PI * t / 100. - PI/2.0);
    return v;
}
SVECT get_analytic_p_3d_rotation(SVECT p0, double t) {
    double r = pow((p0.x-50.0) * (p0.x-50.0) + (p0.y-50.) * (p0.y-50.),0.5);
    double w = PI / 50.0; // radian/second -- time for one orbit, 2 * pi / w = 100s
    SVECT p;
    p.x =  (p0.x - 50.) * cos(w*t) - (p0.y - 50.) * sin(w*t) + 50.0;
    p.y =  (p0.x - 50.) * sin(w*t) + (p0.y - 50.) * cos(w*t) + 50.0;
    p.z = p0.z + (1.0/w) * ( cos(PI/2.0) - cos(0.02*t*PI - PI/2.0) );
    return p;
}
SVECT get_vel_3d_off_grid_node(SVECT p, double t) {
    SVECT v;
    v.x =  1.0;
    v.y =  1.0;
    v.z =  0;
    return v;
}
SVECT get_analytic_p_3d_off_grid_node(SVECT p0, double t) {
    SVECT p;
    if (t < 25) {
        p.x = p0.x + t;
        p.y = p0.y + t;
        p.z = p0.z;
    } else {
        p.x = p0.x + 25;
        p.y = p0.y + 25.;
        p.z = p0.z;
    }
    return p;
}
SVECT get_vel_3d_off_grid_edge(SVECT p, double t) {
    SVECT v;
    if (t<50) {
        v.x =  1.0;
        v.y =  1.0;
        v.z =  0.25;
    } else {
        double r = pow((p.x-50.0) * (p.x-50.0) + (p.y-50.) * (p.y-50.),0.5);
        double theta = atan2 ((p.y-50.),(p.x-50.0));
        double w = PI / 50.0; // radian/second -- time for one orbit, 2 * pi / w = 100s
        v.x =  -w * r * sin(theta);
        v.y =   w * r * cos(theta);
        v.z =   sin(2 * PI * t / 100. - PI/2.0);
    }
    return v;
}
SVECT get_analytic_p_3d_off_grid_edge(SVECT p0, double t) {
    SVECT p;
    if (t < 50) {
        p.x = p0.x + t;
        p.y = p0.y + t;
        p.z = p0.z + 0.25 * t;
        if (p.x > 100) p.x = 100;
        if (p.x < 0)   p.x = 0;
        if (p.y > 100) p.y = 100;
        if (p.y < 0)   p.y = 0;
        if (p.z > 100) p.z = 100;
        if (p.z < 0)   p.z = 0;
    } else {
        p.x = p0.x + 25;
        p.y = p0.y + 25.;
        p.z = p0.z + 0.25 * 25;
    }
    return p;
}
SVECT get_vel_3d_off_grid_face(SVECT p, double t) {
    SVECT v;
    v.x =  1.0;
    v.y =  1.001;
    v.z =  0.25;
    return v;
}
SVECT get_analytic_p_3d_off_grid_face(SVECT p0, double t) {
    SVECT p;
    if (t > 25) t = 25;
    p.x = p0.x + t;
    p.y = p0.y + 1.001 * t;
    p.z = p0.z + 0.25 * t;
    return p;
}
SVECT get_vel_2d_closed_boundary(SVECT p, double t) {
    SVECT v;
    if (t >-1e-12 && t < 50) {
        v.x = 2;
        v.y = 2;
    } else if (t>=50 && t<100) {
        v.x = -2;
        v.y = 2;
    } else if (t>=100 && t<150) {
        v.x = -2;
        v.y = -2;
    } else if (t>=150 && t<200) {
        v.x = 2;
        v.y = -2;
    } else {
        v.x = 2;
        v.y = 2;
    }
    v.z = 0.0;
    return v;
}
SVECT get_analytic_p_2d_closed_boundary(SVECT p0, double t) {
    SVECT p,p49,p99,p149,p199,v1,v2,v3,v4,v5;
    p = p0;
    
    v1.x = 2;  v1.y = 2;
    v2.x = -2; v2.y = 2;
    v3.x = -2; v3.y = -2;
    v4.x = 2;  v4.y = -2;
    v5.x = 2;  v5.y = 2;
    
    double time = 0., dt = 1;
    while (1) {
        if (fabs(time - t) < 1e-10) return  p;
        if (time >-1e-12 && time < 50) {
            p.x += v1.x * dt;
            p.y += v1.y * dt;
            if (p.x > 100) p.x = 100.;
            if (p.x < 0)   p.x = 0.;
            if (p.y > 100) p.y = 100.;
            if (p.y < 0)   p.y = 0.;
        } else if (time>=50 && time<100) {
            p.x += v2.x * dt;
            p.y += v2.y * dt;
            if (p.x > 100) p.x = 100.;
            if (p.x < 0)   p.x = 0.;
            if (p.y > 100) p.y = 100.;
            if (p.y < 0)   p.y = 0.;
        } else if (time>=100 && time<150) {
            p.x += v3.x * dt;
            p.y += v3.y * dt;
            if (p.x > 100) p.x = 100.;
            if (p.x < 0)   p.x = 0.;
            if (p.y > 100) p.y = 100.;
            if (p.y < 0)   p.y = 0.;
        } else if (time>=150 && time<200) {
            p.x += v4.x * dt;
            p.y += v4.y * dt;
            if (p.x > 100) p.x = 100.;
            if (p.x < 0)   p.x = 0.;
            if (p.y > 100) p.y = 100.;
            if (p.y < 0)   p.y = 0.;
        } else {
            p.x += v5.x * dt;
            p.y += v5.y * dt;
            if (p.x > 100) p.x = 100.;
            if (p.x < 0)   p.x = 0.;
            if (p.y > 100) p.y = 100.;
            if (p.y < 0)   p.y = 0.;
        }
        time += dt;
    }
    
}

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
// 3D :: spatially constant grid velocity with vy=vy, linearly oscillating time
SVECT get_vel_test3d_1(SVECT p, double t) {
    SVECT v;
    v.x = 2 * sin(t/10.0);
    v.y = 2 * sin(t/10.0);
    v.z = 2 * sin(t/10.0);
    return v;
}
SVECT get_analytic_p_test3d_1(SVECT p0, double t) {
    SVECT p;
    p.x = p0.x + 20 * (-cos(t/10.) + 1);
    p.y = p0.y + 20 * (-cos(t/10.) + 1);
    p.z = p0.z + 20 * (-cos(t/10.) + 1);
    return p;
}

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
// 3D OYSTER BEHAVOIR

SVECT get_vel_3d_oyster_fall(SVECT p, double t) {
    SVECT v;
    v.x = 0.;
    v.y = 0.;
    v.z = 0.;
    return v;
}
SVECT get_analytic_p_3d_oyster_fall(SVECT p0, double t) {
    return soyster_get_analytic_fall_position(p0,t);
}

SVECT get_vel_3d_oyster_swim(SVECT p, double t) {
    SVECT v;
    v.x = 0.;
    v.y = 0.;
    v.z = 0.;
    return v;
}
SVECT get_analytic_p_3d_oyster_swim(SVECT p0, double t) {
    return soyster_get_analytic_swim_position(p0,t);
}

SVECT get_vel_3d_oyster_tvel(SVECT p, double t) {
    SVECT v;
    v.x = 0.;
    v.y = 0.;
    v.z = 0.;
    return v;
}
SVECT get_analytic_p_3d_oyster_tvel(SVECT p0, double t) {
    return soyster_get_analytic_position(p0,t);
}

SVECT get_wq_3d_oyster_oxygen(SVECT p, double t) {
    SVECT v;
    if (p.z > -70) {
        v.x = 100;
    } else {
        v.x = 49.999999;
    }
    v.y = -1;
    v.z = -1;
    return v;
}

SVECT get_wq_3d_oyster_salinity(SVECT p, double t) {
    SVECT v;
    if (p.z > -70) {
        v.y = 0;
    } else {
        v.y = 30.00001;
    }
    v.x = -1;
    v.z = -1;
    return v;
}

SVECT get_vel_stillwater(SVECT p, double t) {SVECT v; svect_init(&v); return v;}
SVECT get_wq_3d_oyster_sunlight(SVECT p, double t) {
    SVECT v;
    if (p.z > -70) {
        v.z = 100;
    } else {
        v.z = 24.99999;
    }
    v.x = -1;
    v.y = -1;
    return v;
}

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
// SURFACE DPL
SVECT get_vel_3d_surface_dpl(SVECT p, double t) {
    SVECT v;
    v.x = 0.03535;
    v.y = 0.03535;
    //v.z = grid->node[i].z + 100.0) * 0.01 * 1.57079633 * 3 * cos(1.57079633*t);
    v.z = 100.;
    return v;
}

SVECT get_analytic_p_3d_surface_dpl(SVECT p0, double t) {
    SVECT p;
    double dpl_surface = 3*sin(1.57079633*t/150.);
    double dpl = (p0.z + 100.0) * 0.01 * dpl_surface; // AdH accordian style displacement of 1m/s
    p.x = 50.0 + 0.03535 * t;
    p.y = 50.0 + 0.03535 * t;
    p.z = dpl;
    return p;
}

