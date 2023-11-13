#ifndef H_STESTCASE_
#define H_STESTCASE_

typedef struct {

    // 2d hydro
    int nb;
    int winds_and_waves;
    int wet_dry;
    int dam2d;

    // 3d hydro
    int on;             // turn on test case checking
    int coriolis;
    int xcoriolis;
    int ycoriolis;
    int ptest;
    int discharge;
    int slosh;
    int winds_Huang;
    int water_source;
    int lock;
    int cslosh;
    int outflow;
    int tide3d;
    int tide2d;

    // diffusive wave equation
    int dwe_hunter;

    // 2D-3D hydro     // gkc
    int slosh2d3d;
    int slosh2d3d_model_id;
    int tide2d3d;
    int model_id_2d3d;

	// barotropic transport
	int steep_front;
	double init_conc;

    // baroclinic transport
    int salt;
    
    // user input variables
    int nodeID, node1, node2;
    int nodeID_2d, node1_2d, node2_2d;
    int nodeID_3d, node1_3d, node2_3d;
    double a;
    double L;
    double B;
    double H;
    double R;
    double xslice, yslice;
    double theta;
    double hl, hr, dam_location;


} STESTCASE;

/*********************************************************/
/* struct methods -------------------------------------- */

void stestcase_init();


#endif

