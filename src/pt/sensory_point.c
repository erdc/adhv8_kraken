

#include "global_header.h"

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void sensor_pt_init(SENSOR_PT *sp) {
    
    sp->found = FALSE;
    
    svect_init(&sp->r);
    svect_init(&sp->r_old);
    svect_init(&sp->v);
    svect_init(&sp->v_old);
    svect_init(&sp->dpl);
    
    sp->distance = 0.0;
    
    sp->prs = 0.0;
    sp->turb = 0.0;
    sp->accel = 0.0;
    sp->ths = 0.0;
    
    sp->CondVelM = 0.0;
    sp->FVaoMshXYZ[0] = 0.0;
    sp->FVaoMshXYZ[1] = 0.0;
    
    sp->FVaoSVXY = 0.0;
    svect_init(&sp->FVvoSVXYZ);
    
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// input :: double :: SPDistMin - minimum sensor point distance
// input :: int    :: downscale_flag - flag for activiting downscaling of sensory point based on field
// input :: double :: SPDistMean_Max  - max meanof sensory ovoid point distance (m)
// input :: double :: SPDistStdDev_Max - max StdDev of sensory ovoid point distance (m)
// input :: double :: field_scale - (F)ield (V)ariable (for) (S)caling (O)void
// input :: double :: SVaoMshXYZ_NPm1 - the most up-to-date info on where the individual is pointing
// input :: int    :: VldSVOrientXYZ_NPm1 - is xy-plane orientation of the individual's axis valid?
// output :: SENSORY_POINT :: pt;

void sensor_pt_set(SENSOR_PT *pt, double SPDistMin, double SVaoMshXYZ_NPm1, int VldSVOrientXYZ_NPm1, int downscale_flag, double SPDistMean_Max, double SPDistStdDev_Max, double field_scale) {
    
    int i;
    double SPDistMean = 0.0;
    double SPDistStdDev = 0.0;
    
    if (downscale_flag == 1) {
        //Activate downscaling of sensory ovoid point distance based on field variable value
        double SPDistMean_Min           = 1.5;
        double SPDistStdDev_Min         = 1.0;
        double FVforSO_NPm1_Temp        = field_scale;
        double FVforSO_AtSPDistMean_Max = 0.0001;
        double FVforSO_AtSPDistMean_Min = 0.5;
        if (FVforSO_NPm1_Temp < FVforSO_AtSPDistMean_Max) {
            // Sensory ovoid does not get any bigger at values beyond FVforSO_AtSPDistMean_Max
            FVforSO_NPm1_Temp = FVforSO_AtSPDistMean_Max;
        }
        double SPDistMean_Delta   = SPDistMean_Max   - SPDistMean_Min;
        double SPDistStdDev_Delta = SPDistStdDev_Max - SPDistStdDev_Min;
        double FVforSO_Delta      = log10(FVforSO_AtSPDistMean_Max/1E-6) - log10(FVforSO_AtSPDistMean_Min/1E-6);
        double FractBtwDelta      = (log10(FVforSO_NPm1_Temp/1E-6)  - log10(FVforSO_AtSPDistMean_Max/1E-6)) / FVforSO_Delta;
        
        SPDistMean   = SPDistMean_Max   + FractBtwDelta * SPDistMean_Delta;
        SPDistStdDev = SPDistStdDev_Max + FractBtwDelta * SPDistStdDev_Delta;
        
    } else {
        //No downscaling of sensory ovoid point distance
        SPDistMean   = SPDistMean_Max;
        SPDistStdDev = SPDistStdDev_Max;
    }
    
    // calculate a random number with Mean=0, StdDev=1; values can be +/-
    double RRR_Normal01 = find_random_from_normal(0.0,1.0);
    SPDistStdDev = RRR_Normal01 * SPDistStdDev;
    double SPDist = SPDistMean + SPDistStdDev;  // Random SPDist with StdDev centered on the mean
    if (SPDist < SPDistMin) SPDist = SPDistMin;
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Define number and location of sensory points on the sensory ovoid.
    //     positive (+) "x" is in the direction the head of the individual is pointing
    //     negative (-) "x" is in the direction of the individual's tail
    //     positive (+) "y" is to the left of the individual
    //     negative (-) "y" is to the right of the individual
    //     positive (+) "z" is in the direction above the individual (parallel to gravity)
    //     negative (-) "z" is in the direction below the individual (direction of gravity)
    //
    //      Y  (0,1,0)
    //       |
    //       |
    //       |
    //       |
    //       |
    //      --------- X (1,0,0)
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    double SenOvoidX[7], SenOvoidY[7], SenOvoidZ[7];
    
    // First sensory point is the individual's centroid and, therefore, there is no displacement distance:
    SenOvoidX[0] = 0.0;
    SenOvoidY[0] = 0.0;
    SenOvoidZ[0] = 0.0;

    // Sensory point #2 is defined as a displacement distance of (x,y,z)=(SenOvoidX(2),SenOvoidY(2),SenOvoidZ(2))
    //   from the location of the individual as:
    SenOvoidX[1] = SPDist;
    SenOvoidY[1] = 0.0;
    SenOvoidZ[1] = 0.0;

    // Sensory point #3 is defined
    //   from the location of the individual as:
    SenOvoidX[2] = -SPDist;
    SenOvoidY[2] = 0.0;
    SenOvoidZ[2] = 0.0;

    // Sensory point #4 is defined
    //   from the location of the individual as:
    SenOvoidX[3] = 0.0;
    SenOvoidY[3] = SPDist;
    SenOvoidZ[3] = 0.0;

    // Sensory point #5 is defined
    //   from the location of the individual as:
    SenOvoidX[4] = 0.0;
    SenOvoidY[4] = -SPDist;
    SenOvoidZ[4] = 0.0;

    // Sensory point #6 is defined
    //   from the location of the individual as:
    SenOvoidX[5] = 0.0;
    SenOvoidY[5] = 0.0;
    SenOvoidZ[5] = SPDist;

    // Sensory point #7 is defined
    //   from the location of the individual as:
    SenOvoidX[6] = 0.0;
    SenOvoidY[6] = 0.0;
    SenOvoidZ[6] = -SPDist;

    for(i=0; i<7; i++) {
        pt[i].dpl.x = SenOvoidX[i]; // Distance (meters, x-component in Cartesian space) from the individual to sensory point
        pt[i].dpl.y = SenOvoidY[i]; // Distance (meters, y-component in Cartesian space) from the individual to sensory point
        pt[i].dpl.z = SenOvoidZ[i]; // Distance (meters, z-component in Cartesian space) from the individual to sensory point
    }
    
    // Adjust SP distances based on direction from head to tail of individual
    double FSPXYResultant = 0.; // distance from the particle to the SP
    double AbsFSPaoFwd = 0.;    // the absolute angle from the particle to the SP
    double FSPaoFwd = 0.;       // the SP angle relative to the direction pointing from tail to head
    double FSPaoSV = 0.;
    double FSPaoSVRAD = 0.;
    if (VldSVOrientXYZ_NPm1 == 1) {
        // Valid xy-plane orientation of the individual's axis can be determined
        for (i=1; i<5; i++) {
            FSPXYResultant = sqrt( SenOvoidX[i] * SenOvoidX[i] + SenOvoidY[i] * SenOvoidY[i] );
            AbsFSPaoFwd = fabs(atan2(fabs(SenOvoidY[i]),fabs(SenOvoidX[i])));
            if ( (SenOvoidX[i] >= 0.0) && (SenOvoidY[i] >= 0.0)) {
                // FSPaoFwd should be unchanged here
                FSPaoFwd = AbsFSPaoFwd * rad2deg; // convert radians to degrees
            } else if ((SenOvoidX[i] <= 0.0) && (SenOvoidY[i] >= 0.0)) {
                FSPaoFwd = 180.0 - AbsFSPaoFwd * rad2deg;
            } else if ((SenOvoidX[i] <= 0.0) && (SenOvoidY[i] <= 0.0)) {
                FSPaoFwd = 180.0 + AbsFSPaoFwd * rad2deg;
            } else if ((SenOvoidX[i] > 0.0) && (SenOvoidY[i] <= 0.0)) {
                FSPaoFwd = 360.0 - AbsFSPaoFwd * rad2deg;
            }
            
            FSPaoSV = FSPaoFwd + SVaoMshXYZ_NPm1; // add angle of individual
            if (FSPaoSV < 0.0) {
                FSPaoSV = FSPaoSV + 360.0;
            } else if (FSPaoSV > 360.0) {
                FSPaoSV = FSPaoSV - 360.0;
            }
            
            if (FSPaoSV > 180.0) {
                //Ensure that angles go from 0.0 to 180.0 and 0.0 to -180.0
                FSPaoSV -= 360.0;
            }
            if ((FSPaoSV >=  0.0) && (FSPaoSV <= 90.0)) {
                FSPaoSVRAD = FSPaoSV * deg2rad;
                pt[i].dpl.x = FSPXYResultant * cos(FSPaoSVRAD);
                pt[i].dpl.y = FSPXYResultant * cos(FSPaoSVRAD);
            } else if ((FSPaoSV >=  90.0) && (FSPaoSV <= 180.0)) {
                FSPaoSVRAD = (180.0 - FSPaoSV) * deg2rad;
                pt[i].dpl.x = -FSPXYResultant * cos(FSPaoSVRAD);
                pt[i].dpl.y =  FSPXYResultant * cos(FSPaoSVRAD);
            } else if ((FSPaoSV <= 0.0) && (FSPaoSV >= -90.0)) {
                FSPaoSVRAD = -FSPaoSV * deg2rad;
                pt[i].dpl.x =  FSPXYResultant * cos(FSPaoSVRAD);
                pt[i].dpl.y = -FSPXYResultant * cos(FSPaoSVRAD);
            } else if ((FSPaoSV <=  -90.0) && (FSPaoSV >= -180.0)) {
                FSPaoSVRAD = (180.0 + FSPaoSV) * deg2rad;
                pt[i].dpl.x = -FSPXYResultant * cos(FSPaoSVRAD);
                pt[i].dpl.y = -FSPXYResultant * cos(FSPaoSVRAD);
            } else {
                printf("Error: Improper FSPaoSV value\n");
                exit(-1);
            }
        }
    }
    return;
}
