////any constants related to surface water
//typedef struct {
////nodal variable types
//// vertical grid speeds
//#define GS 0
//#define GS_OLD 1
//// bed displacement from sediment transport
//#define BED_DPL 2
//#define BED_DPL_OLD 3
//#define BED_DPL_OLDER 4
//// pressure and pressure pertubations
//#define PRS 5
//#define PRS_PLUS 6
//#define PRS_MINUS 7
////other dependent variables
//#define DENSITY 8            // density
//#define DPL_PERTURBATION 9   // displacement perturbation
//#define CONT_ERROR 10              // nodal continuity errors
//#define VERTICAL_NODE_FLUX 11
//#define HYD_VISCOSITY 12
//#define TRN 13
////must deconstruct into double fields
//#define GRAD_BED_X 14
//#define GRAD_BED_Y 15          // bed gradient
////tangent vector
//#define TAN_VEC_X 16
//#define TAN_VEC_Y 17
////SVECT2D tanvec;
//#define BED_ELEVATION 18 // node.z + bed displacement from sediment
//#define SURFACE_VEL_X 19    //SVECT2D surface_vel;
//#define SURFACE_VEL_Y 20
//#define BOTTOM_VEL_X 21
//#define BOTTOM_VEL_Y 22   //SVECT2D bottom_vel;
//// wind and wave stresses
////not just doubles, need to think about this
////SWIND winds;
////SWAVE waves;
//#define WIND_STRESS_X 23
//#define WIND_STRESS_Y 24
//#define WAVE_RADS_XX 25
//#define WAVE_RADS_XY 26
//#define WAVE_RADS_YY 27
//#define WAVE_STRESS_X 28
//#define WAVE_STRESS_Y 29
//}SSW_DVAR_CODE


//elemental variable types (doubles), none so far

//elemental variable types (integers)
#define WD_FLAG 0
