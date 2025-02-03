#ifndef H_SSW_
#define H_SSW_
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
typedef struct {

    // node and element based variables
    int nsw_nodes;
    int nsw_elems;
    SDVAR dvar;

    // uniform surface water parameters (read from model parameter file)                                   
    double drying_lower_limit;                                                 
    double drying_upper_limit;
    double wd_lower_tol;
    double wd_upper_tol;
    double wd_rate_lower_tol;
    double wd_rate_upper_tol;
    double viscosity;
    double manning_units_constant;
    double density;
    double tau_pg;
    // surface water flags (set defaults)
    int elem_rhs_realloc;

    //doesnt fit in dvar, will use dvar stuff tho (map)
    double **elem_rhs_supg_dacont; //for SW3 only
    double **elem_rhs_supg_cont; // for SW3 only
    double **elem_rhs_dacont_extra_terms; // for SW2 only, FLIPPING ORDER FOR EFFICIENCY [nsw_elems][nnode_on_elem]

    //stuff to store the DVAR column indices
    //nodal variable types
    // vertical grid speeds
    int GS; // 0
    int GS_OLD; //1
    // bed displacement from sediment transport
    int BED_DPL;// 2
    int BED_DPL_OLD;// 3
    int BED_DPL_OLDER;// 4
    // pressure and pressure pertubations
    int PRS;// 5
    int PRS_PLUS;// 6
    int PRS_MINUS;// 7
    //other dependent variables
    int DENSITY;// 8            // density
    int DPL_PERTURBATION;// 9   // displacement perturbation
    int CONT_ERROR;// 10              // nodal continuity errors
    int VERTICAL_NODE_FLUX;// 11
    int HYD_VISCOSITY;// 12
    int TRN;// 13
    //must deconstruct into double fields
    int GRAD_BED_X;// 14
    int GRAD_BED_Y;// 15          // bed gradient
    //tangent vector
    int TAN_VEC_X;// 16
    int TAN_VEC_Y;// 17
    //SVECT2D tanvec;
    int BED_ELEVATION;// 18 // node.z + bed displacement from sediment
    int SURFACE_VEL_X;// 19    //SVECT2D surface_vel;
    int SURFACE_VEL_Y;// 20
    int BOTTOM_VEL_X;// 21
    int BOTTOM_VEL_Y;// 22   //SVECT2D bottom_vel;
    // wind and wave stresses
    //not just doubles, need to think about this
    //SWIND winds;
    //SWAVE waves;
    int WIND_STRESS_X;// 23
    int WIND_STRESS_Y;// 24
    int WAVE_RADS_XX;// 25
    int WAVE_RADS_XY;// 26
    int WAVE_RADS_YY;// 27
    int WAVE_STRESS_X;// 28
    int WAVE_STRESS_Y;// 29

    //elemental variable types (doubles), none so far

    //elemental variable types (integers)
    int WD_FLAG;// 0

} SSW;

void ssw_alloc_init(SSW *sw);
#endif

