#ifndef H_SENSOR_PT_
#define H_SENSOR_PT_

//FSO=2 ==> xy-plane sensory point #2 ~=   0.0  (Front)
//    3 ==> xy-plane sensory point #3 ~= 180.0  (Back)
//    4 ==> xy-plane sensory point #4 ~=  90.0  (Left side)
//    5 ==> xy-plane sensory point #5 ~= -90.0  (Right side)
//    6 ==> vertical sensory point #6 ~=  90.0  (Above)
//    7 ==> vertical sensory point #7 ~= -90.0  (Below)

typedef struct {
    
    SVECT r;      // position
    SVECT r_old;   // last position
    SVECT dpl;
 
    double distance;

    SVECT v_old;  // Velocity at sensor point :: IndvSensoryVelocity_NPm1(1:3,FN)
    SVECT v;     // Velocity at sensor point :: IndvSensoryVelocity_NP(1:3,1:FSOLIMIT,FN)
    
    bool found;   // boolean to relay whether the sensor point is out of bounds
    
    double prs;   // IndvSensoryFieldVars_NP(1,FSO,FN2)   ! Pressure at sensory point
    double turb;  // IndvSensoryFieldVars_NP(2,FSO,FN2)   ! Turbulent kinetic energy at sensory point
    double accel; // IndvSensoryFieldVars_NP(3,FSO,FN2)   ! Acceleration magnitude at sensory point
    double ths;   // IndvSensoryFieldVars_NP(4,FSO,FN2)   ! THS magnitude at sensory point
    
    double CondVelM; // ! Speed (VelM) of immersed, surrounding flow field, e.g., wind or water
    
    double FVaoMshXYZ[2]; // Flow velocity vector vertical angle off xy-plane
    double FVaoSVXY;
    SVECT FVvoSVXYZ;
    
} SENSOR_PT;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// struct methods ++++++++++++++++++++++++++++++++++++++++++
void sensor_pt_init(SENSOR_PT *sp);
void sensor_pt_set(SENSOR_PT *pt, double SPDistMin, double SVaoMshXYZ_NPm1, int VldSVOrientXYZ_NPm1, int downscale_flag, double SPDistMean_Max, double SPDistStdDev_Max, double field_scale);

#endif
