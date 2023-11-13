#include "global_header.h"
#define EPSILON  0.0000001
#define MODULUS(p) (sqrt(p.x*p.x + p.y*p.y + p.z*p.z))
#define TWOPI 6.283185307179586476925287
#define RTOD 57.2957795

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void spolygon_alloc_init(SPOLYGON *poly, int ndim, int np, SVECT *p, double zmax) {
    
    poly->np = np;
    poly->ne = -1;
    poly->ndim = ndim;
    poly->p = (SVECT *)malloc(sizeof(SVECT) * np);
    int i;
    for (i=0; i<poly->np; i++) {
        poly->p[i].x = p[i].x;
        poly->p[i].y = p[i].y;
        poly->p[i].z = p[i].z;
    }
    
    poly->xmax = -1e-6;
    poly->xmin = +1e+6;
    poly->ymax = -1e-6;
    poly->ymin = +1e+6;
    poly->zmax = -1e-6;
    poly->zmin = +1e+6;
    
    // build edge/face connectivity
    if (ndim == 2) {
        
        poly->zmin = -10000000; // for bed to some depth extruded polygons
        poly->zmax = zmax;      // for bed to some depth extruded polygons
        for (i=0; i<poly->np; i++) {
            if (poly->p[i].x > poly->xmax) poly->xmax = poly->p[i].x;
            if (poly->p[i].x < poly->xmin) poly->xmin = poly->p[i].x;
            if (poly->p[i].y > poly->ymax) poly->ymax = poly->p[i].y;
            if (poly->p[i].y < poly->ymin) poly->ymin = poly->p[i].y;
        }
        
        poly->edges = (int **)malloc(sizeof(int *) * np);
        for (i=0; i<np; i++) {
            poly->edges[i] = (int *)malloc(sizeof(int) * 2);
            poly->edges[i][0] = i;
            poly->edges[i][1] = (i+1)%poly->np;
        }
    } else if (ndim == 3) {

    }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void spolygon_geoConform_alloc_init(SPOLYGON *poly, int ndim, int ne, int *elem) {
    
    poly->np = -1;
    poly->ne = ne;
    poly->ndim = ndim;
    poly->elem = (int *)malloc(sizeof(int) * ne);
    int i;
    for (i=0; i<poly->ne; i++) poly->elem[i] = elem[i];
    
    poly->xmax = -1e-6;
    poly->xmin = +1e+6;
    poly->ymax = -1e-6;
    poly->ymin = +1e+6;
    poly->zmax = -1e-6;
    poly->zmin = +1e+6;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void spolygon_free(SPOLYGON *poly) {
    assert(poly != NULL);
    int i;
    if (poly->np > 0) {
        free(poly->p);
        if (poly->ndim == 2) {
            for (i=0; i<poly->np; i++) {
                free(poly->edges[i]);
            }
            free(poly->edges);
        }
    } else {
        assert(poly->ne > 0);
        free(poly->elem);
    }
    free(poly);
    poly=NULL;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double Angle2D(double x1, double y1, double x2, double y2) {
   double theta1 = atan2(y1,x1);
   double theta2 = atan2(y2,x2);
   double dtheta = theta2 - theta1;
   while (dtheta > PI)  dtheta -= TWOPI;
   while (dtheta < -PI) dtheta += TWOPI;
   return(dtheta);
}

int spolygon_containsPoint_2D(SPOLYGON *poly, SVECT p) {
    int i;
    double angle=0;
    SVECT p1,p2;
    for (i=0;i<poly->np;i++) {
        // if point is on polygon vertex, include it
        if (fabs(p.x - poly->p[i].x) + fabs(p.y - poly->p[i].y) < 1e-6) return TRUE;
        
        p1.x = poly->p[i].x - p.x;
        p1.y = poly->p[i].y - p.y;
        p2.x = poly->p[(i+1)%poly->np].x - p.x;
        p2.y = poly->p[(i+1)%poly->np].y - p.y;
        angle += Angle2D(p1.x,p1.y,p2.x,p2.y);
    }
    if (fabs(angle) < PI)
        return(FALSE);
    else
        return(TRUE);
}

// NOTE: reef stops a fixed depth
int spolygon_containsPoint_2D_Extrude(SPOLYGON *poly, SVECT p) {
    int i;
    double angle=0;
    SVECT p1,p2;
    
    if (p.z >= poly->zmax) return(FALSE); // point above polygon extrusion
    
    for (i=0;i<poly->np;i++) {
        // if point is on polygon vertex, include it
        if (fabs(p.x - poly->p[i].x) + fabs(p.y - poly->p[i].y) < 1e-6) return TRUE;
        
        p1.x = poly->p[i].x - p.x;
        p1.y = poly->p[i].y - p.y;
        p2.x = poly->p[(i+1)%poly->np].x - p.x;
        p2.y = poly->p[(i+1)%poly->np].y - p.y;
        angle += Angle2D(p1.x,p1.y,p2.x,p2.y);
    }
    if (fabs(angle) < PI)
        return(FALSE);
    else
        return(TRUE);
}

// mesh conforming! --> just hand particle's element and check
int spolygon_containsPoint(SPOLYGON *poly, int ielem) {
    int i;
    assert(poly->ne > 0);
    assert(ielem > -1);
    for (i=0; i<poly->ne;i++) {
        if (poly->elem[i] == ielem) return(TRUE);
    }
    return(FALSE);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void spolygon_read_file(SPOLYGON *poly, char *pfile) {
    
    int i,npes=1,myid=0;
    double x,y,z;
    SVECT p[1000];
    int ielem[1000];
    char dum[100];
    
#ifdef _MPI
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif
    
    // open model parameter file
    if (myid == 0) {fflush(stdout); printf("OPENING & READING MODEL INPUT FILE: %s\n",pfile);}
    FILE *fp = fopen(pfile,"r");
    
    
    // read parameters
    size_t buffer_size = 80;
    char *buffer = malloc(buffer_size * sizeof(char));
    char delim[] = " ";  // line delimiter
    char *ptr = NULL;    // to read the card on line
    
    // read to count either polygon vertices or the # of elements that comprise it
    int np=0, ne=0, ndim=-1;
    double zmax = -10000.;
    while(getline(&buffer, &buffer_size, fp) != -1) {
        ptr = strtok(buffer, delim);
        //printf("%s ", ptr);
        if (strcmp(ptr,"NDIM") == 0) {
            ptr = strtok(NULL, delim);
            ndim = atoi(ptr);
        }
        if (strcmp(ptr,"PT") == 0) {
            sscanf(buffer,"%s %lf %lf %lf ",dum,&p[np].x,&p[np].y,&p[np].z);
            np++;
        }
        if (strcmp(ptr,"ELEM") == 0) {
            sscanf(buffer,"%s %d ",dum,&ielem[ne]);
            ne++;
        }
        if (strcmp(ptr,"ZMAX") == 0) {
            ptr = strtok(NULL, delim);
            zmax = atoi(ptr);
        }
    }
    rewind(fp);
    
    assert(ndim==2 || ndim==3);
    assert(np>0 || ne>0);
    assert(np<1000 && ne<1000);
    
    // allocate polygon
    if (np > 0) {
        spolygon_alloc_init(poly,ndim,np,p,zmax);
    } else {
        spolygon_geoConform_alloc_init(poly,ndim,ne,ielem);
    }
    
    fclose(fp);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void spolygon_test() {
    SVECT pt;              // test point
    SPOLYGON *poly = NULL; // test polygon
    SVECT p[100];          // test polygon vertices
    
    printf("-- TESTING POINT IN POLYGON\n");

    // ++++++++++++++++++++++++++++++++++++++++++++
    // create a 2D polygon with 6 control points
    // ++++++++++++++++++++++++++++++++++++++++++++
    p[0].x =  0.0;    p[0].y = 0.0;     p[0].z = 0.0;
    p[1].x =  1.0;    p[1].y = 1.0;     p[1].z = 0.0;
    p[2].x = -1.0;    p[2].y = 3.0;     p[2].z = 0.0;
    p[3].x = -5.0;    p[3].y = 4.0;     p[3].z = 0.0;
    p[4].x = -8.0;    p[4].y = -1.0;    p[4].z = 0.0;
    p[5].x = -3.0;    p[5].y = -2.0;    p[5].z = 0.0;
    //p[6].x = p[0].x;  p[6].y = p[0].y;  p[6].z = p[0].z;
    poly = (SPOLYGON *)malloc(sizeof(SPOLYGON));
    spolygon_alloc_init(poly, 2, 6, p, -50.0);

    // check that a point on polygon vertex is inside
    pt.x = p[0].x;
    pt.y = p[0].y;
    pt.z = p[0].z;
    printf("test1\n"); assert(spolygon_containsPoint_2D(poly,pt));

    // check that a point on polygon edge is inside
    pt.x = 0.5 * (poly->p[0].x + poly->p[1].x);
    pt.y = 0.5 * (poly->p[0].y + poly->p[1].y);
    pt.z = p[0].z;
    printf("test2\n"); assert(spolygon_containsPoint_2D(poly,pt));

    // check that a point inside is inside
    pt.x = -3.0;
    pt.y = -1.0;
    pt.z = p[0].z;
    printf("test3\n"); assert(spolygon_containsPoint_2D(poly,pt));

    // check that a point outside is ouside
    pt.x = 100;
    pt.y = 100;
    pt.z = p[0].z;
    printf("test4\n"); assert(spolygon_containsPoint_2D(poly,pt)==FALSE);
    spolygon_free(poly);
    
    // try a more complicated 2D polygon
    p[0].x =  0.0;    p[0].y = 0.0;    p[0].z = 0.0;
    p[1].x =  0.0;    p[1].y = 6.0;    p[1].z = 0.0;
    p[2].x = -9.0;    p[2].y = 6.0;    p[2].z = 0.0;
    p[3].x = -9.0;    p[3].y = 1.0;    p[3].z = 0.0;
    p[4].x = -7.0;    p[4].y = 1.0;    p[4].z = 0.0;
    p[5].x = -7.0;    p[5].y = 4.0;    p[5].z = 0.0;
    p[6].x = -5.0;    p[6].y = 4.0;    p[6].z = 0.0;
    p[7].x = -5.0;    p[7].y = 2.0;    p[7].z = 0.0;
    p[8].x = -4.0;    p[8].y = 2.0;    p[8].z = 0.0;
    p[9].x = -4.0;    p[9].y = 4.0;    p[9].z = 0.0;
    p[10].x = -2.0;   p[10].y = 4.0;   p[10].z = 0.0;
    p[11].x = -2.0;   p[11].y = 0.0;   p[11].z = 0.0;
    poly = (SPOLYGON *)malloc(sizeof(SPOLYGON));
    spolygon_alloc_init(poly, 2, 12, p, -50.0);

    // check that a point on polygon vertex is inside
    pt.x = p[0].x;
    pt.y = p[0].y;
    pt.z = p[0].z;
    printf("test1\n"); assert(spolygon_containsPoint_2D(poly,pt));

    // check that a point on polygon edge is inside
    pt.x = -1;
    pt.y =  1;
    pt.z = p[0].z;
    printf("test2\n"); assert(spolygon_containsPoint_2D(poly,pt));

    // check that a point inside is inside
    pt.x = -3.0;
    pt.y =  2.0;
    pt.z = p[0].z;
    printf("test3\n"); assert(spolygon_containsPoint_2D(poly,pt)==FALSE);

    // check that a point outside is ouside
    pt.x = -4.5;
    pt.y = 2.5;
    pt.z = p[0].z;
    printf("test4\n"); assert(spolygon_containsPoint_2D(poly,pt));
    spolygon_free(poly);
    
    printf("|| PASSED\n");
    //exit(-1);
}
