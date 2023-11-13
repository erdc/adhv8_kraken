#include "global_header.h"

static int firstcall = 1; // gkc

void stestcase_init() {

    test_case_flag.on =  OFF;
    
    test_case_flag.nb = OFF;
    test_case_flag.ptest = OFF;
    test_case_flag.coriolis = OFF;
    test_case_flag.xcoriolis = OFF;
    test_case_flag.ycoriolis = OFF;
    test_case_flag.discharge = OFF;
    test_case_flag.slosh = OFF;
    test_case_flag.winds_Huang = OFF;
    test_case_flag.salt = OFF;
    test_case_flag.water_source = OFF;
    test_case_flag.lock = OFF;
    test_case_flag.winds_and_waves = OFF;
    test_case_flag.wet_dry = OFF;
    test_case_flag.dam2d = OFF;
    test_case_flag.cslosh = OFF;
    test_case_flag.outflow = OFF;
    test_case_flag.tide3d = OFF;
    test_case_flag.tide2d = OFF;

    test_case_flag.dwe_hunter = OFF;


    test_case_flag.slosh2d3d = OFF; // gkc
    test_case_flag.slosh2d3d_model_id = UNSET_INT; // gkc
     test_case_flag.model_id_2d3d = UNSET_INT;
     test_case_flag.tide2d3d = OFF;

    // user input variables
    test_case_flag.nodeID = UNSET_INT;
    test_case_flag.nodeID_2d = UNSET_INT;
    test_case_flag.nodeID_3d = UNSET_INT;
    test_case_flag.a = 0.0;
    test_case_flag.L = 0.0;
    test_case_flag.B = 0.0;
    test_case_flag.H = 0.0;
    test_case_flag.R = 0.0;

    test_case_flag.node1 = UNSET_INT;
    test_case_flag.node1_2d = UNSET_INT;
    test_case_flag.node1_3d = UNSET_INT;
    test_case_flag.node2 = UNSET_INT;
    test_case_flag.node2_2d = UNSET_INT;
    test_case_flag.node2_3d = UNSET_INT;
    test_case_flag.xslice = 0.0;
    test_case_flag.yslice = 0.0;
    test_case_flag.theta = 0.0;

    test_case_flag.hl = 0.;
    test_case_flag.hr = 0.;
    test_case_flag.dam_location = 0.;
    firstcall = 0;
} 

