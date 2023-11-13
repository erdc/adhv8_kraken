#ifndef H_POLYGON_
#define H_POLYGON_


typedef struct {
    int ndim; // polygon dimension
    int np;   // number of polygon vertices
    int ne;   // number of elements that make polygon (geo mesh conforming)
    SVECT *p; // polygon vertices
    int **edges; // segment edge connectivity for 2D polygons
    int **faces; // triangular face connectivity for 3D polygons
    int *elem;  // geo elements that comprise polygon
    double xmax,xmin,ymax,ymin,zmax,zmin; // polygon bounds
} SPOLYGON;

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// struct methods +++++++++++++++++++++++++++++++++++++++++

void spolygon_alloc_init(SPOLYGON *poly, int ndim, int np, SVECT *p, double zmax);
void spolygon_geoConform_alloc_init(SPOLYGON *poly, int ndim, int ne, int *elem);
int spolygon_containsPoint_2D(SPOLYGON *poly, SVECT p1);
int spolygon_containsPoint_2D_Extrude(SPOLYGON *poly, SVECT p);
void spolygon_read_file(SPOLYGON *poly, char *pfile);
void spolygon_test();
#endif
