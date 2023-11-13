#include "global_header.h"
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void PreferredDirection_FromVector(SENSOR_SPECIES *sp, int MoveWithAgainstVector, int DEBUG);
void PreferredDirection_FromGradient(SENSOR_SPECIES *sp, double *CondVarX,int DEBUG);
void UnstrMeshWallEffect(SENSOR_SPECIES *sp, double dt, int DEBUG);
void VolMovInOtherFrames(SENSOR_SPECIES *sp, int DEBUG);
void Decision_UsherMcClelland2001(SENSOR_SPECIES *sp,double dt,int *BehavChosen,bool DEBUG);
void Decision_Anderson2002(SENSOR_SPECIES *sp, double *Outp_Accumulator,int *BehavChosen,bool DEBUG);
void VectorRelation(int VectorType, SVECT MVComponentsMsh, bool *VldRVOrientXYZ, bool *VldMVOrientXYZ, double *RVaoMshXYZ, SVECT RVvoMshXYZ, double *MVaoRVXY, double *MVaoMshXYZ, SVECT *MVvoRVXYZ);
void PreferredDirection_FromSensoryPts(SENSOR_PT *sp, double *BVvaSP, bool DEBUG, double *SVaoSVXY, double *SVaoMshXYZ);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void sensor_species_alloc_init(SENSOR_SPECIES **ssp) {
    
    (*ssp) = (SENSOR_SPECIES *) malloc(sizeof(SENSOR_SPECIES));
    SENSOR_SPECIES *sp = (*ssp);
    
    int i;
    
    for (i=0; i<NIDATTRB; i++) {
        sp->IdAttrb[0] = 0;
        sp->IdAttrb_old[0] = 0;
    }
    
    for (i=0; i<NSENSORPT_MAX; i++) sensor_pt_init(&sp->sensor[i]);
        
    for (i=0; i<NBEHAVS_MAX; i++) {
        sp->initialize_flag[i] = -999;
        sp->initialize_flag_old[i] = -999;
        sp->Info_Accumulator[i] = 0.0;
        sp->Info_Accumulator_old[i] = 0.0;
        sp->Evid_Acclimatized[i] = 0.0;
        sp->Evid_Acclimatized_old[i] = 0.0;
        sp->Evid_Acclimatized_new[i] = 0.0;
        sp->Evid_Threshold[i] = 0.0;
        sp->Evid_Value[i] = 0.0;
        sp->Evid_Activity[i] = 0.0;
    }
    
    sp->NumDecisions = 0;
    sp->IdSpdRes = 0.0;
    sp->IdSpdRes_old = 0.0;
    
    sp->SVaoMshXYZ[0] = 0.0;
    sp->SVaoMshXYZ[1] = 0.0;
    sp->SVaoMshXYZ_old[0] = 0.0;
    sp->SVaoMshXYZ_old[1] = 0.0;
    sp->SVvoMshXYZ.x = 0.0;
    sp->SVvoMshXYZ.y = 0.0;
    sp->SVvoMshXYZ.z = 0.0;
    sp->SVvoMshXYZ_old.x = 0.0;
    sp->SVvoMshXYZ_old.y = 0.0;
    sp->SVvoMshXYZ_old.z = 0.0;
    
    sp->SVaoFVXY = 0.0;
    sp->SVaoFVXY_old = 0.0;
    sp->SVvoSVXYZ.x = 0.0;
    sp->SVvoSVXYZ.y = 0.0;
    sp->SVvoSVXYZ.z = 0.0;
    sp->SVvoFVXYZ.x = 0.0;
    sp->SVvoFVXYZ.y = 0.0;
    sp->SVvoFVXYZ.z = 0.0;
    sp->SVvoFVXYZ_old.x = 0.0;
    sp->SVvoFVXYZ_old.y = 0.0;
    sp->SVvoFVXYZ_old.z = 0.0;
    
    sp->VldSVOrientXYZ[0] = false;
    sp->VldSVOrientXYZ[1] = false;
    sp->VldSVOrientXYZ_old[0] = false;
    sp->VldSVOrientXYZ_old[1] = false;
    
    sp->SVaoSVXY = 0.0;
    sp->SVaoSVXY_old = 0.0;
    
    sp->decision_algorithm = 1;
    sp->sensing_algorithm = 1;
    
    stimuli_perception_init(&sp->stimuli);
    pt_interact_init(&sp->interact);
    wall_effect_init(&sp->wallEffect);
    behavoir_trans_usher_init(&sp->trans_usher);
    behavoir_trans_anderson_init(&sp->trans_anderson);
    sspeed_init(&sp->speed);
    sensing_orientation_ornstein_init(&sp->sorient_ornstein);
    sensing_orientation_codling_init(&sp->sorient_codling);
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  The purpose of this subroutine is to computer the following 3 variables:
//  SVaoSVXY      = the volitional movement vector angle off of (relative to) the previous volitional movement vector in the horizontal
//                     plane. This value must range from -180.0 to 180.0 and will be used in the Call function below to calculate SVaoMshXYZ_NP[0],
//                     which is the volitional movement vector angle off of (relative to) the mesh xy-plane axes.
//  sp->SVaoMshXYZ[1]= the vertical volitional movement vector angle off of (relative to) the mesh xy-plane itself; in other words, this is the
//                     vertical angle off the horizontal plane. This value must range from -90.0 to 90.0.
//  IdSpdRes      = movement velocity (resultant value).

void sensor_species_get_behavior(SENSOR_SPECIES *sp, double dt) {
    
    bool DEBUG = FALSE;
    int i;
    double kappa_MinForSubroutine = 1e-6;
    int Seed = 187827;
    
    sp->IdAttrb[0] = sp->IdAttrb_old[0];
    
    // make sure sensor point acceleration is positive (why?)
    for (i=0; i<NSENSORPT_MAX; i++) if (sp->sensor[i].accel < 2E-6) sp->sensor[i].accel = 2E-6;
    
    
    // initialize behavoirs if flagged
    if (sp->initialize_flag[0] == -999) {
        for (i=0; i<NBEHAVS_MAX; i++) sp->initialize_flag_old[i] = 0;
        sp->SVaoSVXY_old = 0.0;
        sp->SVaoMshXYZ_old[0] = 0.0;
        sp->SVaoMshXYZ_old[1] = 0.0;
        svect_init(&sp->SVvoMshXYZ_old);
        sp->IdSpdRes_old = 0.0;
        sp->Evid_Acclimatized[1]  = log10(sp->sensor[0].accel / 1E-6);
        if (sp->Evid_Acclimatized[1] < sp->stimuli.AccIM_min) sp->Evid_Acclimatized[1] = sp->stimuli.AccIM_min;
        sp->Evid_Acclimatized[3] = sp->sensor[0].r.z; // elevation at individual location at current t
        sp->Evid_Acclimatized[5] = 0.80;
        if (sp->decision_algorithm == 1) {
            for (i=0; i<NBEHAVS_MAX; i++) sp->Info_Accumulator_old[i] = 0.0;
        }
        if (sp->decision_algorithm == 2) {
            for (i=0; i<NBEHAVS_MAX; i++) sp->Info_Accumulator_old[i] = sp->trans_anderson.prob_reward;
        }
        sp->NumDecisions = 0;
    }
    
    if ((sp->sensor[3].found == FALSE) && (sp->sensor[4].found == FALSE)) {
        // Left/right sensory points out-of-bounds in 2-D laterally-averaged mesh data
        for (i=0; i<NSENSORPT_MAX; i++) {
            sp->sensor[i].v.y = 0.0;
            if (i == 3 || i == 4) {
                sp->sensor[i].v.x = 0.0;
                sp->sensor[i].v.y = 0.0;
                sp->sensor[i].v.z = 0.0;
            }
        }
    } else if ((sp->sensor[5].found == FALSE) && (sp->sensor[6].found == FALSE)) {
        // Vertical sensory points out-of-bounds in 2-D depth-averaged mesh data
        for (i=0; i<NSENSORPT_MAX; i++) {
            sp->sensor[i].v.z = 0.0;
            if (i == 5 || i == 6) {
                sp->sensor[i].v.x = 0.0;
                sp->sensor[i].v.y = 0.0;
                sp->sensor[i].v.z = 0.0;
            }
        }
    }
    
    sp->SVaoFVXY = 0.0;
    svect_init(&sp->SVvoFVXYZ);
    svect_init(&sp->SVvoSVXYZ);
    sp->VldSVOrientXYZ[0] = false;
    sp->VldSVOrientXYZ[1] = false;
    sp->VldSVOrientXYZ_old[0] = false;
    sp->VldSVOrientXYZ_old[1] = false;
    sp->VldFVOrientXYZ[0] = false;
    sp->VldFVOrientXYZ[1] = false;
    
    // Flow velocity vector xy-angle off swim vector axis (i.e., flow velocity xy-angle relative to volitional angle)
    // Flow velocity components off moving swim vector axis (i.e., velocity components as if individual is stationary)
    for (i=0; i<NSENSORPT_MAX; i++) {
        sp->sensor[i].FVaoSVXY = 0.0;
        sp->sensor[i].FVvoSVXYZ.x = 0.0;
        sp->sensor[i].FVvoSVXYZ.y = 0.0;
        sp->sensor[i].FVvoSVXYZ.z = 0.0;
        
        sp->sensor[i].FVaoMshXYZ[0] = 0.0;
        sp->sensor[i].FVaoMshXYZ[1] = 0.0;
    }
    
    for (i=0; i<NBEHAVS_MAX; i++) {
        sp->Evid_Value[i] = 0.0;
        sp->initialize_flag[i] = 0;
    }
    sp->NumDecisions++;
    
    //****************************************************************************************************
    //*******************************************   Section   ********************************************
    //****************************************************************************************************
    
    //Flow Vector Velocity Relative to Individual's Volitional Movement Vector Velocity.
    if (fabs(sp->SVvoMshXYZ_old.x) < 1e-12 && fabs(sp->SVvoMshXYZ_old.y) < 1e-12) {
        sp->VldSVOrientXYZ_old[0] = false; //No valid xy-plane orientation of the individual's axis => no individual movement from t-1 to t
    } else {
        sp->VldSVOrientXYZ_old[1] = true; //// Valid xy-plane orientation of the individual's axis can be determined for t-1 to t
    }
    
    for (i=0; i<NSENSORPT_MAX; i++) {
        if (sp->sensor[i].found == FALSE) continue;
        VectorRelation(1, sp->sensor[i].v, sp->VldSVOrientXYZ_old, sp->VldFVOrientXYZ, sp->SVaoMshXYZ_old,sp->SVvoMshXYZ_old, &sp->sensor[i].FVaoSVXY, sp->sensor[i].FVaoMshXYZ, &sp->sensor[i].FVvoSVXYZ);
    }
    
    //****************************************************************************************************
    //*******************************************   Section   ********************************************
    //****************************************************************************************************
    
    //Default behavior.
    for (i=0; i<NBEHAVS_MAX; i++) sp->Evid_Activity[i] = 0.0;
    sp->Evid_Activity[0] = 1.0;
    
    // Environment at the 7 sensory points, the individual's centroid and each cardinal direction.
    double CondAcclM = 0.0, dv = 0.0;;
    SVECT dr, Accl_i;
    for (i=0; i<NSENSORPT_MAX; i++) {
        if (sp->sensor[i].found == TRUE){
            sp->sensor[i].CondVelM = svect_mag(sp->sensor[i].v);
        } else {
            //  if 2-D mesh data, set non-existent sensory pts to same value as centroid position
            sp->sensor[i].CondVelM = sp->sensor[0].CondVelM;
        }
        CondAcclM = sp->sensor[0].accel; // Acceleration (AcclM) of flow field, e.g., wind or water, from field variable
        double dv = 0.0;
        
        dr.x = sp->sensor[0].r.x - sp->sensor[0].r_old.x;
        dr.y = sp->sensor[0].r.y - sp->sensor[0].r_old.y;
        dr.z = sp->sensor[0].r.z - sp->sensor[0].r_old.z;
        
        dv = sp->sensor[0].v.x - sp->sensor[0].v_old.x;
        Accl_i.x = dv / dt + sp->sensor[0].v.x * dv / dr.x + sp->sensor[0].v.y * dv / dr.y + sp->sensor[0].v.z * dv / dr.z;
        
        dv = sp->sensor[0].v.y - sp->sensor[0].v_old.y;
        Accl_i.y = dv / dt + sp->sensor[0].v.y * dv / dr.x + sp->sensor[0].v.y * dv / dr.y + sp->sensor[0].v.z * dv / dr.z;
        
        dv = sp->sensor[0].v.z - sp->sensor[0].v_old.z;
        Accl_i.z = dv / dt + sp->sensor[0].v.z * dv / dr.x + sp->sensor[0].v.y * dv / dr.y + sp->sensor[0].v.z * dv / dr.z;
        
        // Acceleration (AcclM) of flow field, e.g., wind or water; material derivative from individual's trajectory
        CondAcclM = svect_mag(Accl_i);
    }
    
    // WoE for behavior 2 (swim to faster water; is change in AcclM > kB{2}).
    sp->Evid_Value[1] = ( log10(CondAcclM / 1E-6) - sp->Evid_Acclimatized[1] ) / sp->Evid_Acclimatized[1];
    sp->Evid_Activity[1] = fabs(sp->Evid_Value[1]/ sp->Evid_Threshold[1]);
    
    // WoE for behavior 3 (swim against water flow vector; is change in AcclM > kB{3}).
    // Evid_Value[2]    = ( Log10(CondAcclM / 1E-6) - Evid_Acclimatized[2] ) / Evid_Acclimatized[2]
    sp->Evid_Activity[2] = fabs(sp->Evid_Value[1] / sp->Evid_Threshold[2]);
    
    
    // WoE for behavior 4 (swim to acclimatized depth; is |change| in Elevation/Pressure > kB{4}).
    if ((sp->sensor[5].found == FALSE) && (sp->sensor[6].found == FALSE)) {
        sp->Evid_Activity[3] = 0.0;// Vertical sensory points out-of-bounds in 2-D depth-averaged mesh data
    } else {
        sp->Evid_Value[3] = fabs(sp->sensor[0].r.z - sp->Evid_Acclimatized[3]);
        sp->Evid_Activity[3] = fabs(sp->Evid_Value[3] / sp->Evid_Threshold[3]);
    }
    
    
    //****************************************************************************************************
    //*******************************************   Section   ********************************************
    //****************************************************************************************************
    
    // Update individual's integrated condition.
    sp->Evid_Acclimatized_new[1] = (1.0 - sp->stimuli.Evid_AcclimRate[1]) * log10(CondAcclM / 1E-6) + sp->stimuli.Evid_AcclimRate[1] * sp->Evid_Acclimatized[1];
    sp->Evid_Acclimatized_new[3] = (1.0 - sp->stimuli.Evid_AcclimRate[3]) * sp->sensor[0].r.z + sp->stimuli.Evid_AcclimRate[3] * sp->Evid_Acclimatized[3];
    
    //****************************************************************************************************
    //*******************************************   Section   ********************************************
    //****************************************************************************************************
    
    // Behavior transition.
    int BehavChosen[2] = {-99,-99};
    double Outp_Accumulator[NBEHAVS_MAX]; // For output display purposes
    if (sp->decision_algorithm == 1) {
        // returns Info_Accumulator and BehavChosen
        Decision_UsherMcClelland2001(sp,dt,BehavChosen,DEBUG);
        for (i=0; i<NBEHAVS_MAX; i++) Outp_Accumulator[i] = sp->Info_Accumulator[i];
        // Info_Accumulator_old <=> Evid_Accumulated_old
        // Info_Accumulator   <=> Evid_Accumulated      // For storage (old) purposes
        // Outp_Accumulator   <=> Evid_Accumulated      // For output display purposes
    } else if (sp->decision_algorithm == 2) {
        // returns Info_Accumulator, Outp_Accumulator and BehavChosen
        Decision_Anderson2002(sp,Outp_Accumulator,BehavChosen,DEBUG);
        // Info_Accumulator_old <=> sp->Info_Accumulator_old
        // Info_Accumulator   <=> sp->Info_Accumulator
        // Outp_Accumulator   <=> Util_Behav
    } else {
        printf("Error: Invalid decision algorithm.n");
        exit(-1);
    }
    
    //****************************************************************************************************
    //*******************************************   Section   ********************************************
    //****************************************************************************************************
    
    //Behavior.
    //  Three forms of response:
    //  [0] use of a vector as the basis for response. For instance, behavior
    //      could be moving at an angle relative to the direction of water flow.
    //  [1] use of a gradient as the basis for response. For instance, behavior
    //      could be moving at an angle relative to the gradient in water speed.
    //  [2] use of a discrete entity (object) as the basis for response. For instance, behavior
    //      could be moving at an angle relative to the direction of a food item or predator.
    
    //      BehavChosen[0] = 2
    //      BehavChosen[1] = 0
    
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>> Behavior 1.
    //>>>>>>>>>>>>>>>
    double d_xy = 0.0, d_z = 0.0, c_xy = 0.0, c_z = 0.0;
    double kappa_xy = 0.0, kappa_z = 0.0, lambda_xy = 0.0, lambda_z = 0.0;
    if (BehavChosen[0] == 1) {
        // B1 (swim with water flow vector)
        for (i=0; i<NBEHAVS_MAX; i++) sp->initialize_flag[i] = 0;
        sp->initialize_flag[0] = 1;
        
        //  ---> Orientation preference.
        int MoveWithAgainstVector = 1; // 1 = with ; -1 against
        // output SVaoMshXYZ, SVaoSVXY
        PreferredDirection_FromVector(sp,MoveWithAgainstVector,DEBUG);
        
        if (sp->sensor[0].CondVelM < sp->stimuli.velM) {
            // slack flow
            //  ---> Orientation stochasticity (Codling et al., 2004).
            //   xy
            d_xy     = sp->sorient_codling.slack_d_xy[0];
            kappa_xy = sp->sorient_codling.slack_kappa_xy[0];  // orienting ability of the individual (Codling et al., 2004)
            //   z
            d_z      = sp->sorient_codling.slack_d_z[0];     // Dble(Coefficients(5,37));
            kappa_z  = sp->sorient_codling.slack_kappa_z[0]; // Dble(Coefficients(5,45));
            //  ---> Orientation stochasticity (Ornstein-Uhlenbeck).
            //   xy
            lambda_xy = sp->sorient_ornstein.slack_lambda_xy[0];  // Coefficients(7,33);
            c_xy      = sp->sorient_ornstein.slack_c_xy[0];       // Coefficients(7,41);
            //   z
            lambda_z  = sp->sorient_ornstein.slack_lambda_z[0];  // Coefficients(7,37);
            c_z       = sp->sorient_ornstein.slack_c_z[0];       // Coefficients(7,45);
        } else {
            //  ---> Orientation stochasticity (Codling et al., 2004).
            //   xy
            d_xy     = sp->sorient_codling.d_xy[0];     // Dble(Coefficients(5,1));
            kappa_xy = sp->sorient_codling.kappa_xy[0]; // Dble(Coefficients(5,17));
            //   z
            d_z      = sp->sorient_codling.d_z[0];     // Dble(Coefficients(5,9));
            kappa_z  = sp->sorient_codling.kappa_z[0]; // Dble(Coefficients(5,25));
            //  ---> Orientation stochasticity (Ornstein-Uhlenbeck).
            //   xy
            lambda_xy = sp->sorient_ornstein.lambda_xy[0];  // DCoefficients(7,1);
            c_xy      = sp->sorient_ornstein.c_xy[0];       // DCoefficients(7,17);
            //   z
            lambda_z  = sp->sorient_ornstein.lambda_z[0];  // Coefficients(7,9);
            c_z       = sp->sorient_ornstein.c_z[0];       // Coefficients(7,25);
        }
        
        //  ---> Speed
        double RRR_Normal01 = find_random_from_normal(0.0,1.0); // calculate a random number with Mean=0, StdDev=1; values can be +/-
        double StdDev = RRR_Normal01 * (sp->speed.cruise - sp->speed.drift) / sp->speed.speed_coeff1[0];
        sp->IdSpdRes = sp->speed.length * (sp->speed.cruise  + StdDev); //Mean centered on cruise speed
        
        //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        //>>> Behavior 2.
        //>>>>>>>>>>>>>>>
    } else if (BehavChosen[0] == 2) {
        //// B2 (move to faster water)
        for (i=0; i<NBEHAVS_MAX; i++) sp->initialize_flag[i] = 0;
        sp->initialize_flag[1] = 1;
        
        //  ---> Orientation preference.
        //       CondVarX = CondVelM
        int SPFound[NSENSORPT_MAX];
        double CondVarX[NSENSORPT_MAX];
        for (i=0; i<NSENSORPT_MAX; i++) {
            SPFound[i] = sp->sensor[i].found;
            CondVarX[i] = sp->sensor[i].ths;
        }
        // output SVaoMshXYZ,SVaoSVXY
        PreferredDirection_FromGradient(sp,CondVarX,DEBUG);
        
        if ((sp->initialize_flag_old[1] == 1) && (sp->Evid_Acclimatized_new[1] < sp->Evid_Acclimatized[1])) {
            // slack flow
            //  ---> Orientation stochasticity (Codling et al., 2004). // TimeStep-1 ==> Beh 2 ; AcclM ==> improving
            //   xy
            d_xy     = sp->sorient_codling.slack_d_xy[1];
            kappa_xy = sp->sorient_codling.slack_kappa_xy[1];
            //   z
            d_z      = sp->sorient_codling.slack_d_z[1];     // Dble(Coefficients(5,37));
            kappa_z  = sp->sorient_codling.slack_kappa_z[1]; // Dble(Coefficients(5,45));
            //  ---> Orientation stochasticity (Ornstein-Uhlenbeck).
            //   xy
            lambda_xy = sp->sorient_ornstein.slack_lambda_xy[1];  // Coefficients(7,33);
            c_xy      = sp->sorient_ornstein.slack_c_xy[1];       // Coefficients(7,41);
            //   z
            lambda_z  = sp->sorient_ornstein.slack_lambda_z[1];  // Coefficients(7,37);
            c_z       = sp->sorient_ornstein.slack_c_z[1];       // Coefficients(7,45);
            
        } else {
            //  ---> Orientation stochasticity (Codling et al., 2004). // TimeStep-1 ==> not Beh 2 <or> Beh 2 ; AcclM ==> not improving
            //   xy
            d_xy     = sp->sorient_codling.d_xy[1];     // Dble(Coefficients(5,3))
            kappa_xy = sp->sorient_codling.kappa_xy[1]; // Dble(Coefficients(5,19))
            //   z
            d_z      = sp->sorient_codling.d_z[1];     // Dble(Coefficients(5,11))
            kappa_z  = sp->sorient_codling.kappa_z[1]; // Dble(Coefficients(5,27))
            //  ---> Orientation stochasticity (Ornstein-Uhlenbeck).
            //   xy
            lambda_xy = sp->sorient_ornstein.lambda_xy[1];  // Coefficients(7,3)
            c_xy      = sp->sorient_ornstein.c_xy[1];       // Coefficients(7,19)
            //   z
            lambda_z  = sp->sorient_ornstein.lambda_z[1];  // Coefficients(7,11)
            c_z       = sp->sorient_ornstein.c_z[1];       // Coefficients(7,27)
        }
        
        //  ---> Speed
        double RRR_Normal01 = find_random_from_normal(0.0,1.0); // calculate a random number with Mean=0, StdDev=1; values can be +/-
        double StdDev = RRR_Normal01 * (sp->speed.cruise - sp->speed.drift);
        sp->IdSpdRes = sp->speed.length * (sp->speed.cruise  + StdDev); //Mean centered on cruise speed
        
        //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        //>>> Behavior 3.
        //>>>>>>>>>>>>>>>
    } else if (BehavChosen[0] == 3) {
        // B3 (swim against water flow vector)
        for (i=0; i<NBEHAVS_MAX; i++) sp->initialize_flag[i] = 0;
        sp->initialize_flag[2] = 1;
        
        //  ---> Orientation preference.
        int MoveWithAgainstVector = -1; // 1 = with ; -1 against
        // output SVaoMshXYZ, SVaoSVXY
        PreferredDirection_FromVector(sp,MoveWithAgainstVector,DEBUG);
        
        if ((sp->initialize_flag_old[2] == 1) && (sp->Evid_Acclimatized_new[2] < sp->Evid_Acclimatized[2])) {
            // slack flow
            //  ---> Orientation stochasticity (Codling et al., 2004).
            //   xy
            d_xy     = sp->sorient_codling.slack_d_xy[2];
            kappa_xy = sp->sorient_codling.slack_kappa_xy[2];
            //   z
            d_z      = sp->sorient_codling.slack_d_z[2];     // Dble(Coefficients(5,37));
            kappa_z  = sp->sorient_codling.slack_kappa_z[2]; // Dble(Coefficients(5,45));
            //  ---> Orientation stochasticity (Ornstein-Uhlenbeck).
            //   xy
            lambda_xy = sp->sorient_ornstein.slack_lambda_xy[2];  // Coefficients(7,33);
            c_xy      = sp->sorient_ornstein.slack_c_xy[2];       // Coefficients(7,41);
            //   z
            lambda_z  = sp->sorient_ornstein.slack_lambda_z[2];  // Coefficients(7,37);
            c_z       = sp->sorient_ornstein.slack_c_z[2];       // Coefficients(7,45);
            
            //  ---> Speed
            sp->IdSpdRes = sp->IdSpdRes_old * sp->speed.value_reward[2]; //Coefficients(4,11) TimeStep-1 ==> Beh 3 ; AcclM ==> improving
            
        } else {
            //  ---> Orientation stochasticity (Codling et al., 2004).
            //   xy
            d_xy     = sp->sorient_codling.d_xy[2];     // Dble(Coefficients(5,3))
            kappa_xy = sp->sorient_codling.kappa_xy[2]; // Dble(Coefficients(5,19))
            //   z
            d_z      = sp->sorient_codling.d_z[2];     // Dble(Coefficients(5,11))
            kappa_z  = sp->sorient_codling.kappa_z[2]; // Dble(Coefficients(5,27))
            //  ---> Orientation stochasticity (Ornstein-Uhlenbeck).
            //   xy
            lambda_xy = sp->sorient_ornstein.lambda_xy[2];  // Coefficients(7,3)
            c_xy      = sp->sorient_ornstein.c_xy[2];       // Coefficients(7,19)
            //   z
            lambda_z  = sp->sorient_ornstein.lambda_z[2];  // Coefficients(7,11)
            c_z       = sp->sorient_ornstein.c_z[2];       // Coefficients(7,27)
            
            //  ---> Speed
            sp->IdSpdRes = sp->sensor[0].CondVelM  * sp->speed.speed_coeff1[2]; //Coefficients(4,3) TimeStep-1 ==> not Beh 3 <or> Beh 3 ; AcclM ==> not improving
        }
        
    } else {
        
        printf("Error: BehavChosen[0] = %d\n",BehavChosen[0]);
        exit(-1);
        
    }
    
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>> Sensing Orientation.
    //>>>>>>>>>>>>>>>
    double RRR = 0.;
    if (sp->sensing_algorithm == 1) {
        
        double mu = 0.; // Mean turning angle
        double Angle_vonMises = 0.; // Angle is between -pi and pi
        
        //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        //>>> Orientation (Codling et al., 2004).
        //>>>>>>>>>>>>>>>
        //   xy
        mu = -d_xy * -sp->SVaoSVXY * deg2rad; // Codling et al. (2004). Reference direction of movement last time step is 0.0d0 in this coordinate system (individual's heading is 0-deg)
        if (mu      < -180.0) mu = mu + 360.0;
        if (mu      >  180.0) mu = mu - 360.0;
        if (fabs(mu) >  180.0) {printf("ERROR: Abs(mu) > 180.0\n"); exit(-1);}
        // Subroutine to calculate a random number (angle in radians)
        if (kappa_xy >= kappa_MinForSubroutine) Angle_vonMises = RandomFromvonMises(mu,kappa_xy,Seed);
        //! Subroutine to calculate a uniform random number btw [0,1]; values only + :: kappa = 0 collapses the von Mises distribution to a wrapped uniform distribution
        if (kappa_xy < kappa_MinForSubroutine) RRR = find_random_in_range(0.0, 1.0);
        if (kappa_xy < kappa_MinForSubroutine) Angle_vonMises = (RRR * 2.0 - 1.0) * PI; // kappa = 0 collapses the von Mises distribution to a wrapped uniform distribution
        sp->SVaoSVXY = Angle_vonMises * rad2deg;
        //   z
        mu = -d_z  * -sp->SVaoSVXY * deg2rad; // Codling et al. (2004). Reference direction of movement last time step is 0.0d0 in this coordinate system (horizontal plane is 0-deg)
        if (mu < -90.0) mu = -90.0 - (mu + 90.0);
        if (mu >  90.0) mu =  90.0 - (mu - 90.0);
        if (fabs(mu) >  90.0) {printf("ERROR: Abs(mu) > 90.0\n"); exit(-1);}
        // Subroutine to calculate a random number (angle in radians)
        if (kappa_z >= kappa_MinForSubroutine) Angle_vonMises = RandomFromvonMises(mu,kappa_z,Seed);
        if (kappa_z < kappa_MinForSubroutine) RRR = find_random_in_range(0.0, 1.0);
        if (kappa_z < kappa_MinForSubroutine) Angle_vonMises = (RRR * 2.0 - 1.0) * (PI*0.50); //kappa = 0 collapses the von Mises distribution to a wrapped uniform distribution
        
        // Note that vertical angles have a practical range of 0.0 to |45.0| per JEB paper, 'Swimming speeds and buoyancy
        // compensation of migrating adult chum salmon Oncorhynchus keta revealed by speed/depth/acceleration data logger'
        sp->SVaoMshXYZ[1] = Angle_vonMises * rad2deg;
        
    } else if (sp->sensing_algorithm ==  2) {
        //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        //>>>             (Ornstein-Uhlenbeck).
        //>>>>>>>>>>>>>>>
        RRR =  find_random_in_range(0.0, 1.0);        double Angle_WrappedUniform = (RRR * 2.0 - 1.0) * 180.0;
        //   xy
        // Direction of movement last time step is 0.0d0 in this coordinate system (individual's heading is 0-deg)
        // Unit time step "1.0" used because TimeStepLength makes porting calibration btw projects w/different time steps difficult.
        sp->SVaoSVXY  = 0.0 + lambda_xy * sp->SVaoSVXY * 1.0 + c_xy * Angle_WrappedUniform * sqrt(1.0);
        if (sp->SVaoSVXY < -360.0) sp->SVaoSVXY = fmod(sp->SVaoSVXY,-360.0);
        if (sp->SVaoSVXY >  360.0) sp->SVaoSVXY = fmod(sp->SVaoSVXY, 360.0);
        //   z
        // 0.5 to cut angle range from +/-180 to +/-90
        // Unit time step "1.0" used because TimeStepLength makes porting calibration btw projects w/different time steps difficult.
        Angle_WrappedUniform = Angle_WrappedUniform * 0.50;
        sp->SVaoMshXYZ[1] = sp->SVaoMshXYZ_old[1] + lambda_z  * (sp->SVaoMshXYZ[1] - sp->SVaoMshXYZ_old[1]) * 1.0 + c_z  * Angle_WrappedUniform * sqrt(1.0);
        
        while (sp->SVaoMshXYZ[1] < -90.0) sp->SVaoMshXYZ[1] = -90.0 - (sp->SVaoMshXYZ[1] + 90.0);
        while (sp->SVaoMshXYZ[1] <  90.0) sp->SVaoMshXYZ[1] =  90.0 - (sp->SVaoMshXYZ[1] - 90.0);
        
    } else {
        //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        //>>> Checks.
        //>>>>>>>>>>>>>>>
        printf("Error: Invalid xyz sensing orientation algorithm.");
        exit(-1);
    }
    
    if (sp->SVaoSVXY < -180.0) sp->SVaoSVXY      = sp->SVaoSVXY + 360.0;
    if (sp->SVaoSVXY > 180.0) sp->SVaoSVXY      = sp->SVaoSVXY - 360.0;
    if (sp->SVaoMshXYZ[1] < -90.0) sp->SVaoMshXYZ[1] = -90.0 - (sp->SVaoMshXYZ[1] + 90.0);
    if (sp->SVaoMshXYZ[1] > 90.0) sp->SVaoMshXYZ[1] =  90.0 - (sp->SVaoMshXYZ[1] - 90.0);
    if (sp->SVaoSVXY < -180.0) {printf("Error in Subroutine BehaviorRule (B',BehavChosen[0],'): SVaoSVXY < -180.0\n"); exit(-1);}
    if (sp->SVaoSVXY > 180.0) {printf("Error in Subroutine BehaviorRule (B',BehavChosen[0],'): SVaoSVXY >  180.0\n"); exit(-1);}
    if (sp->SVaoMshXYZ[1] < -90.0) {printf("Error in Subroutine BehaviorRule (B',BehavChosen[0],'): sp->SVaoMshXYZ[1]< -90.0\n"); exit(-1);}
    if (sp->SVaoMshXYZ[1] > 90.0 ) {printf("Error in Subroutine BehaviorRule (B',BehavChosen[0],'): sp->SVaoMshXYZ[1]>  90.0\n"); exit(-1);}
    if (sp->IdSpdRes < sp->speed.length * sp->speed.drift) sp->IdSpdRes = sp->speed.length * sp->speed.drift;
    if (sp->IdSpdRes > sp->speed.length * sp->speed.burst) sp->IdSpdRes = sp->speed.length * sp->speed.burst;
    
    
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //>>> Behavior 4 (vertical).
    //>>>>>>>>>>>>>>>
    if (BehavChosen[1] == 4) {
        // B4 (move to acclimatized depth)
        sp->initialize_flag[3] = 1;
        
        //  ---> Orientation preference.
        double Evid_Value4_old = fabs(sp->sensor[0].r_old.z  - sp->Evid_Acclimatized_old[3]);
        double Evid_Value4   = fabs(sp->sensor[0].r.z      - sp->Evid_Acclimatized[3]  );
        
        // TimeStep-1 ==> not Beh 4 <or> Beh 4 ; depth ==> not improving
        int Up = 1;
        if (sp->sensor[0].r.z > sp->Evid_Acclimatized[3] ) Up = -Up; //// Shallower than acclimatized
        double SVaoMshXYTemp = sp->SVaoMshXYZ_old[1]+ Up * sp->stimuli.vert_increase;  // Vertical change in orientation when not driven by flow vector
        if((sp->sensor[0].r.z > sp->Evid_Acclimatized[3]) && (-sp->sensor[0].FVaoMshXYZ[1] < SVaoMshXYTemp)) {
            // Shallower than acclimatized and likely due to vertical flow vector
            SVaoMshXYTemp = -sp->sensor[0].FVaoMshXYZ[1] * sp->stimuli.flow_drivenl;
        } else if ((sp->sensor[0].r.z < sp->Evid_Acclimatized[3]) && (-sp->sensor[0].FVaoMshXYZ[1] > SVaoMshXYTemp)) {
            // Deeper    than acclimatized and likely due to vertical flow vector
            SVaoMshXYTemp = -sp->sensor[0].FVaoMshXYZ[1] * sp->stimuli.flow_drivenl;
        }
        if ((Up < 0.0)  && (sp->SVaoMshXYZ[1] > SVaoMshXYTemp)) sp->SVaoMshXYZ[1] = SVaoMshXYTemp; // Shallower than acclimatized; don't override a more aggressive response from another (favorable) behavior
        if ((Up > 0.0)  && (sp->SVaoMshXYZ[1] < SVaoMshXYTemp)) sp->SVaoMshXYZ[1] = SVaoMshXYTemp; // Deeper    than acclimatized; don't override a more aggressive response from another (favorable) behavior
        
        // TimeStep-1 ==> Beh 4 ; depth ==> improving
        if ((sp->initialize_flag[3] == 1) && (Evid_Value4 < Evid_Value4_old)) SVaoMshXYTemp = sp->SVaoMshXYZ_old[1] * sp->stimuli.decay_rate; //Coefficients(6,20);
        //  ---> Orientation stochasticity (Codling et al., 2004).
        //   z
        d_z     = sp->sorient_codling.d_z[3];      // TimeStep-1 ==> not Beh 4 <or> Beh 4 ; depth ==> not improving
        kappa_z = sp->sorient_codling.kappa_z[3]; // TimeStep-1 ==> not Beh 4 <or> Beh 4 ; depth ==> not improving
        if ((sp->initialize_flag[3] == 1) && (Evid_Value4 < Evid_Value4_old)) d_z = sp->sorient_codling.slack_d_z[3];     // TimeStep-1 ==> Beh 4 ; depth ==> improving
        if ((sp->initialize_flag[3] == 1) && (Evid_Value4 < Evid_Value4_old)) kappa_z = sp->sorient_codling.slack_kappa_z[3]; // TimeStep-1 ==> Beh 4 ; depth ==> improving
        //  ---> Orientation stochasticity (Ornstein-Uhlenbeck).
        //   z
        lambda_z = sp->sorient_ornstein.lambda_z[3]; // TimeStep-1 ==> not Beh 4 <or> Beh 4 ; depth ==> not improving
        c_z      = sp->sorient_ornstein.c_z[3];     // TimeStep-1 ==> not Beh 4 <or> Beh 4 ; depth ==> not improving
        if ((sp->initialize_flag[3] == 1) && (Evid_Value4 < Evid_Value4_old)) lambda_z = sp->sorient_ornstein.slack_lambda_z[3]; // TimeStep-1 ==> Beh 4 ; depth ==> improving
        if ((sp->initialize_flag[3] == 1) && (Evid_Value4 < Evid_Value4_old)) c_z  = sp->sorient_ornstein.slack_c_z[3];     // TimeStep-1 ==> Beh 4 ; depth ==> improving
        //  ---> Speed.
        double IdSpdResTemp = sp->sensor[0].CondVelM  * sp->speed.speed_coeff1[3]; // TimeStep-1 ==> not Beh 4 <or> Beh 4 ; AcclM ==> not improving
        if ((sp->initialize_flag[3] == 1) && (Evid_Value4 < Evid_Value4_old)) IdSpdResTemp = sp->IdSpdRes_old * sp->speed.value_reward[3]; // TimeStep-1 ==> Beh 4 ; depth ==> improving
        if (IdSpdResTemp > sp->IdSpdRes) sp->IdSpdRes = IdSpdResTemp; // Don't override a more aggressive response from another (favorable) behavior
        
        //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        //>>> Sensing Orientation.
        //>>>>>>>>>>>>>>>
        if (sp->sensing_algorithm == 1) {
            double mu = 0.0, Angle_vonMises = 0.;
            //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //>>> Orientation (Codling et al., 2004).
            //>>>>>>>>>>>>>>>
            
            //   z
            mu = -d_z  * -sp->SVaoMshXYZ[1] * deg2rad; // Codling et al. (2004). Reference direction of movement last time step is 0.0d0 in this coordinate system (horizontal plane is 0-deg)
            if (mu < -90.0) mu = -90.0 - (mu + 90.0);
            if (mu >  90.0) mu =  90.0 - (mu - 90.0);
            if (fabs(mu) >  90.0) {printf("ERROR: Abs(mu) > 90.0\n"); exit(-1);}
            // Subroutine to calculate a random number (angle in radians)
            if (kappa_z >= kappa_MinForSubroutine) Angle_vonMises = RandomFromvonMises(mu,kappa_z,Seed);
            if (kappa_z < kappa_MinForSubroutine) RRR = find_random_in_range(0.0, 1.0);
            if (kappa_z < kappa_MinForSubroutine) Angle_vonMises = (RRR * 2.0 - 1.0) * (PI*0.50); //kappa = 0 collapses the von Mises distribution to a wrapped uniform distribution
            
            // Note that vertical angles have a practical range of 0.0 to |45.0| per JEB paper, 'Swimming speeds and buoyancy
            // compensation of migrating adult chum salmon Oncorhynchus keta revealed by speed/depth/acceleration data logger'
            sp->SVaoMshXYZ[1] = Angle_vonMises * rad2deg;
            
        } else if (sp->sensing_algorithm ==  2) {
            //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //>>>             (Ornstein-Uhlenbeck).
            //>>>>>>>>>>>>>>>
            
            //   z
            RRR = find_random_in_range(0.0, 1.0); // Subroutine to calculate a uniform random number btw [0,1]; values only +
            double Angle_WrappedUniform = (RRR * 2.0 - 1.0) * 90.0;
            // Unit time step "1.0" used because TimeStepLength makes porting calibration btw projects w/different time steps difficult.
            sp->SVaoMshXYZ[1] = sp->SVaoMshXYZ_old[1] + lambda_z  * (sp->SVaoMshXYZ[1] - sp->SVaoMshXYZ_old[1]) * 1.0 + c_z  * Angle_WrappedUniform * sqrt(1.0);
            
            while (sp->SVaoMshXYZ[1] < -90.0) sp->SVaoMshXYZ[1] = -90.0 - (sp->SVaoMshXYZ[1] + 90.0);
            while (sp->SVaoMshXYZ[1] <  90.0) sp->SVaoMshXYZ[1] =  90.0 - (sp->SVaoMshXYZ[1] - 90.0);
            
        } else {
            
            //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //>>> Checks.
            //>>>>>>>>>>>>>>>
            printf("Error: Invalid z sensing orientation algorithm.\n");
            exit(-1);
        }
    }
    if (sp->SVaoMshXYZ[1] < -90.0) sp->SVaoMshXYZ[1] = -90.0 - (sp->SVaoMshXYZ[1]+ 90.0);
    if (sp->SVaoMshXYZ[1] >  90.0) sp->SVaoMshXYZ[1] =  90.0 - (sp->SVaoMshXYZ[1]- 90.0);
    if (sp->SVaoMshXYZ[1] < -90.0) {printf("Error in Subroutine BehaviorRule (B',BehavChosen[1],'): sp->SVaoMshXYZ[1]< -90.0");exit(-1);}
    if (sp->SVaoMshXYZ[1] > 90.0) {printf("Error in Subroutine BehaviorRule (B',BehavChosen[1],'): sp->SVaoMshXYZ[1]>  90.0");exit(-1);}
    if (sp->IdSpdRes < sp->speed.length * sp->speed.drift) sp->IdSpdRes = sp->speed.length * sp->speed.drift;
    if (sp->IdSpdRes > sp->speed.length * sp->speed.burst) sp->IdSpdRes = sp->speed.length * sp->speed.burst;
    
    //****************************************************************************************************
    //*******************************************   Section   ********************************************
    //****************************************************************************************************
    
    if ((sp->sensor[3].found == FALSE) && (sp->sensor[4].found == FALSE)) {  // Left/right sensory points out-of-bounds in 2-D laterally-averaged mesh data
        if (fabs(sp->SVaoSVXY) < 90.0) sp->SVaoSVXY =   0.0;
        if (fabs(sp->SVaoSVXY) > 90.0) sp->SVaoSVXY = 180.0;
        
    }
    
    if ((sp->sensor[5].found == FALSE) && (sp->sensor[6].found == FALSE)) { // Vertical sensory points out-of-bounds in 2-D depth-averaged mesh data
        sp->SVaoMshXYZ[1] = 0.0;
    }
    
    UnstrMeshWallEffect(sp,dt,DEBUG);
    
    //Check SVaoSVXY and SVaoMshXYZ[1].
    if ((sp->SVaoSVXY < -180.0) || (sp->SVaoSVXY >  180.0)) {
        printf("Error: sp->SVaoSVXY calculation\n");
        printf("  sp->SVaoSVXY = ',sp->SVaoSVXY\n");
        exit(-1);
    }
    if ((sp->SVaoMshXYZ[1] < -90.0) || (sp->SVaoMshXYZ[1] > 90.0)) {
        printf("Error: sp->SVaoMshXYZ[1]calculation\n");
        printf("  sp->SVaoMshXYZ[1]= ',SVaoMshXYZ[1]\n");
        exit(-1);
    }
    
    VolMovInOtherFrames(sp,DEBUG);
    return;
}


//*************************************************************************
//*************************************************************************
//****                                                                 ****
//****             Subroutine Decision_UsherMcClelland2001             ****
//****                                                                 ****
//****                 developed by R. Andrew Goodwin                  ****
//****                                                                 ****
//*************************************************************************
//*************************************************************************

void Decision_UsherMcClelland2001(SENSOR_SPECIES *sp,double dt,int *BehavChosen,bool DEBUG) {
    int i,j,Inh;
    
    for (i=0; i<2; i++) BehavChosen[i] = 0; // 1st position => xyz-plane behavior,  2nd position => vertical behavior
    
    // Behavior transition (Usher and McClelland, 2001)
    // Maximum value of Evid_Accumulated is the behavior implemented:
    double Evid_Inhibiting[NBEHAVS_MAX]; for (i=0; i<NBEHAVS_MAX; i++) Evid_Inhibiting[i] = 0.0;
    double Evid_Leak[NBEHAVS_MAX];
    double Evid_Accumulated_old_Inhibiting[NBEHAVS_MAX];
    double Noise[NBEHAVS_MAX];
    double Evid_Accumulated_Max = 0.0;
    double Evid_k                  = sp->trans_usher.Evid_Decay[0];         // Coefficients(8,3)  // Use Coeffs(8, 4-10) if Evid_k differs btw behaviors
    double Evid_InhibitWgt         = sp->trans_usher.Evid_Inhibit[0];       // Coefficients(8,11) // Use Coeffs(8,12-18) if Evid_InhibitWgt differs btw behaviors
    double c_LCA                   = sp->trans_usher.white_noise_weight[0]; // Coefficients(8,19)  // Use Coeffs(8,20-26) if c_LCA differs btw behaviors
    
    // Subroutine to calculate a random number with Mean=0, StdDev=1; values can be +/-
    double RRR_Normal01 = find_random_from_normal(0.0,1.0);
    for (i=0; i<4; i++) {
        
        // --Diff
        if (floor(sp->Evid_Activity[i]) == 0) sp->Evid_Activity[i] = 0.0; // Evidence activity doesn't contribute to evidence for behavior switching until after threshold is exceeded
        Evid_Accumulated_old_Inhibiting[i] = 0.0;
        for (Inh=0; Inh<3; Inh++) { //Inh=1,3  ! =1,NBehavs ! B4 does not inhibit B1-B3, but B1-B3 do inhibit B4
            if (Inh != i) {
                Evid_Accumulated_old_Inhibiting[i] += sp->Info_Accumulator_old[Inh];
            }
        }
        Evid_Inhibiting[i] = Evid_InhibitWgt * Evid_Accumulated_old_Inhibiting[i];
        Evid_Leak[i] = Evid_k * sp->Info_Accumulator_old[i];
        // Unit time step "1.0" because using
        sp->Info_Accumulator[i] = sp->Info_Accumulator_old[i] + (sp->Evid_Activity[i] - Evid_Leak[i] - Evid_Inhibiting[i]) * 1.0 + c_LCA * RRR_Normal01 * sqrt(1.0);//   TimeStepLength complicates porting btw projects.
        Noise[i] = c_LCA * RRR_Normal01 * sqrt(1.0); // To show in v_Debug.txt
        if (sp->Info_Accumulator[i] < 0.0) sp->Info_Accumulator[i] = 0.0;
        
        // --Diff
        if ((sp->Info_Accumulator[i] > Evid_Accumulated_Max) && (i != 3)) {
            BehavChosen[0]          = i;
            Evid_Accumulated_Max = sp->Info_Accumulator[i];
        }
    }
    if (BehavChosen[0] == 0) BehavChosen[0] = 1;
    if (sp->Info_Accumulator[3] >= sp->Info_Accumulator[BehavChosen[0]]) BehavChosen[1] = 4;
    if (sp->IdAttrb[0] != 0) { // Attribute #1: 0=Indv ; 1=Invertebrate ; -1=To be removed from simulation
        BehavChosen[0] = 1;
        BehavChosen[1] = 0;
    }
    if (BehavChosen[0] == 6) BehavChosen[1] = 0;
    if (DEBUG) {
        printf("\n");
        printf("   => IdAttrb = %d\n",sp->IdAttrb[0]);
        printf("   => BehavChosen[0] = %d\n",BehavChosen[0]);
        printf("      BehavChosen[1] [Vertical Override] = %d\n",BehavChosen[1]);
        printf("      sp->Info_Accumulator[0] = %f\n",sp->Info_Accumulator[0]);
        printf("      sp->Info_Accumulator[1] = %f\n",sp->Info_Accumulator[1]);
        printf("      sp->Info_Accumulator[2] = %f\n",sp->Info_Accumulator[2]);
        printf("      sp->Info_Accumulator[3] = %f\n",sp->Info_Accumulator[3]);
        printf("      sp->Info_Accumulator_old[0] = %f\n",sp->Info_Accumulator_old[0]);
        printf("      sp->Info_Accumulator_old[1] = %f\n",sp->Info_Accumulator_old[1]);
        printf("      sp->Info_Accumulator_old[2] = %f\n",sp->Info_Accumulator_old[2]);
        printf("      sp->Info_Accumulator_old[3] = %f\n",sp->Info_Accumulator_old[3]);
        printf("      Evid_Activity[0] = %f\n",sp->Evid_Activity[0]);
        printf("      Evid_Activity[1] = %f\n",sp->Evid_Activity[1]);
        printf("      Evid_Activity[2] = %f\n",sp->Evid_Activity[2]);
        printf("      Evid_Activity[3] = %f\n",sp->Evid_Activity[3]);
        printf("      Evid_Leak[0] = %f\n",Evid_Leak[0]);
        printf("      Evid_Leak[1] = %f\n",Evid_Leak[1]);
        printf("      Evid_Leak[2] = %f\n",Evid_Leak[2]);
        printf("      Evid_Leak[3] = %f\n",Evid_Leak[3]);
        printf("      Evid_k = %f\n",Evid_k);
        printf("      Evid_Inhibiting[0] = %f\n",Evid_Inhibiting[0]);
        printf("      Evid_Inhibiting[1] = %f\n",Evid_Inhibiting[1]);
        printf("      Evid_Inhibiting[2] = %f\n",Evid_Inhibiting[2]);
        printf("      Evid_Inhibiting[3] = %f\n",Evid_Inhibiting[3]);
        printf("      Evid_InhibitWgt  = %f\n",Evid_InhibitWgt);
        printf("      Evid_Accumulated_old_Inhibiting[0]  = %f\n",Evid_Accumulated_old_Inhibiting[0]);
        printf("      Evid_Accumulated_old_Inhibiting[1]  = %f\n",Evid_Accumulated_old_Inhibiting[1]);
        printf("      Evid_Accumulated_old_Inhibiting[2]  = %f\n",Evid_Accumulated_old_Inhibiting[2]);
        printf("      Evid_Accumulated_old_Inhibiting[3]  = %f\n",Evid_Accumulated_old_Inhibiting[3]);
        printf("      Noise[0] = %f\n",Noise[0]);
        printf("      Noise[1] = %f\n",Noise[1]);
        printf("      Noise[2] = %f\n",Noise[2]);
        printf("      Noise[3] = %f\n",Noise[3]);
        printf("      c_LCA = %f\n",c_LCA);
        
    }
    return;
}


//*************************************************************************
//*************************************************************************
//****                                                                 ****
//****                Subroutine Decision_Anderson2002                 ****
//****                                                                 ****
//****                 developed by R. Andrew Goodwin                  ****
//****                                                                 ****
//*************************************************************************
//*************************************************************************
void Decision_Anderson2002(SENSOR_SPECIES *sp, double *Outp_Accumulator,int *BehavChosen,bool DEBUG) {
    
    int i, Beh, Evt, NEvents;
    
    //Behavior transition (Anderson, 2002).
    // Maximum value of Outp_Accumulator is the behavior implemented:
    for (i=0; i<NBEHAVS_MAX; i++) Outp_Accumulator[i] = 0.0;
    double Outp_Accumulator_Max   = 0.0;
    for (i=0; i<2; i++) BehavChosen[i] = 0; // 1st position => xyz-plane behavior,  2nd position => vertical behavior
    
    double Evid_Accumulator[4], Value_Reward[4];
    Evid_Accumulator[0] = sp->trans_anderson.Evid_Accumulator[0]; //Coefficients(2,3)   // Evidence accumulator for Behavior 1 integrating past and present evidence
    Evid_Accumulator[1] = sp->trans_anderson.Evid_Accumulator[1]; //Coefficients(2,4)   //                                   2
    Evid_Accumulator[2] = sp->trans_anderson.Evid_Accumulator[2]; //Coefficients(2,5)   //                                   3
    Evid_Accumulator[3] = sp->trans_anderson.Evid_Accumulator[3]; //Coefficients(2,6)   //                                   4
    Value_Reward[0]      = sp->trans_anderson.value_reward[0]; //Coefficients(2,11)  // Value (intrinsic) of Behavior 1 reward, if successful
    Value_Reward[1]      = sp->trans_anderson.value_reward[1]; //Coefficients(2,12)  //                               2
    Value_Reward[2]      = sp->trans_anderson.value_reward[2]; //Coefficients(2,13)  //                               3
    Value_Reward[3]      = sp->trans_anderson.value_reward[3]; //Coefficients(2,14)  //                               4
    
    for (Beh=0; Beh<4; Beh++) {
        //--Diff
        NEvents = floor(sp->Evid_Activity[Beh]);
        if (floor(sp->Evid_Activity[Beh]) == 0) NEvents = 1;
        if (floor(sp->Evid_Activity[Beh]) == 0) sp->Evid_Activity[Beh] = 0.0;
        if (floor(sp->Evid_Activity[Beh]) > 0)  sp->Evid_Activity[Beh] = 1.0;
        for (Evt=0; Evt<NEvents; Evt++) {
            if (Evt == 1) sp->Info_Accumulator[Beh] = (1.0 - Evid_Accumulator[Beh]) * sp->Evid_Activity[Beh] + Evid_Accumulator[Beh] * sp->Info_Accumulator_old[Beh];
            if (Evt > 1)  sp->Info_Accumulator[Beh] = (1.0 - Evid_Accumulator[Beh]) * sp->Evid_Activity[Beh] + Evid_Accumulator[Beh] * sp->Info_Accumulator[Beh];
        }
        Outp_Accumulator[Beh] = sp->Info_Accumulator[Beh] * Value_Reward[Beh]; // Can add here the cost of pursuing the reward, whether successful or not
        //--Diff
        if ((Outp_Accumulator[Beh] > Outp_Accumulator_Max) && (Beh != 3)) {
            BehavChosen[0] = Beh;
            Outp_Accumulator_Max = Outp_Accumulator[Beh];
        }
    }
    
    if (Outp_Accumulator[3] >= sp->trans_anderson.beh4_threshold) BehavChosen[1] = 4;
    if (sp->IdAttrb[0] != 0) {  // Attribute #1: 0=Indv ; 1=Invertebrate ; -1=To be removed from simulation
        BehavChosen[0] = 1;
        BehavChosen[1] = 0;
    }
    if (BehavChosen[0] == 6) BehavChosen[1] = 0;
    
    if (DEBUG) {
        printf(" \n");
        printf("   => IdAttrb                         = %d\n",sp->IdAttrb[0]);
        printf("   => BehavChosen[0]                     = %d\n",BehavChosen[0]);
        printf("      BehavChosen[1] [Vertical Override] = %d\n",BehavChosen[1]);
        printf("        Outp_Accumulator[0]                 = %f\n",Outp_Accumulator[0]);
        printf("        Outp_Accumulator[1]                 = %f\n",Outp_Accumulator[1]);
        printf("        Outp_Accumulator[2]                 = %f\n",Outp_Accumulator[2]);
        printf("        Outp_Accumulator[3]                 = %f\n",Outp_Accumulator[3]);
        printf("        Outp_Accumulator(5)                 = %f\n",Outp_Accumulator[5]);
        printf("        Outp_Accumulator(6)                 = %f\n",Outp_Accumulator[6]);
        printf("        Outp_Accumulator(7)                 = %f\n",Outp_Accumulator[7]);
        printf("        Outp_Accumulator(8)                 = %f\n",Outp_Accumulator[8]);
        printf("          sp->Evid_Activity[0]            = %f\n",sp->Evid_Activity[0]);
        printf("          sp->Evid_Activity[1]            = %f\n",sp->Evid_Activity[1]);
        printf("          sp->Evid_Activity[2]            = %f\n",sp->Evid_Activity[2]);
        printf("          sp->Evid_Activity[3]            = %f\n",sp->Evid_Activity[3]);
        printf("          sp->Evid_Activity(5)            = %f\n",sp->Evid_Activity[5]);
        printf("          sp->Evid_Activity(6)            = %f\n",sp->Evid_Activity[6]);
    }
    
    return;
}

//*************************************************************************
//*************************************************************************
//****                                                                 ****
//****            Subroutine PreferredDirection_FromVector             ****
//****                                                                 ****
//****                 developed by R. Andrew Goodwin                  ****
//****                                                                 ****
//*************************************************************************
//*************************************************************************

// output :: SVaoSVXY - o Local scalar || sp->SVaoMshXYZ - o Local array
void PreferredDirection_FromVector(SENSOR_SPECIES *sp, int MoveWithAgainstVector, int DEBUG) {
    double SVaoMshXYTemp = 0.0;
    
    if (MoveWithAgainstVector < 0) {
        SVaoMshXYTemp = sp->sensor[0].FVaoMshXYZ[0] - 180.0;
        if (SVaoMshXYTemp <   0.0) SVaoMshXYTemp += 360.0;
        if (SVaoMshXYTemp > 360.0) SVaoMshXYTemp -= 360.0;
    } else {
        SVaoMshXYTemp = sp->sensor[0].FVaoMshXYZ[0];
    }
    sp->SVaoSVXY = SVaoMshXYTemp - sp->SVaoMshXYZ_old[0];
    if (sp->SVaoSVXY   <    0.0) sp->SVaoSVXY +=  + 360.0;
    if (sp->SVaoSVXY   >  360.0) sp->SVaoSVXY -= 360.0;
    if (sp->SVaoSVXY   >  180.0) sp->SVaoSVXY -= 360.0; // Make angles go from 0.0 to 180.0 and 0.0 to -180.0
    if (sp->SVaoSVXY   < -180.0) {printf("Error in Subroutine PreferredDirection_FromVector: sp->*SVaoSVXY < -180.0\n"); exit(-1);}
    if (sp->SVaoSVXY   >  180.0) {printf("Error in Subroutine PreferredDirection_FromVector: sp->*SVaoSVXY >  180.0\n"); exit(-1);}
    sp->SVaoMshXYZ[1] = sp->sensor[0].FVaoMshXYZ[1] * MoveWithAgainstVector;
    
    return;
}


//*************************************************************************
//*************************************************************************
//****                                                                 ****
//****           Subroutine PreferredDirection_FromGradient            ****
//****                                                                 ****
//****                 developed by R. Andrew Goodwin                  ****
//****                                                                 ****
//*************************************************************************
//*************************************************************************

void PreferredDirection_FromGradient(SENSOR_SPECIES *sp, double *CondVarX,int DEBUG) {
    
    if ((CondVarX[0] >= 0.0) && (CondVarX[1] >= 0.0) && (CondVarX[2] >= 0.0) && (CondVarX[3] >= 0.0) && (CondVarX[4] >= 0.0) && (CondVarX[5] >= 0.0) && (CondVarX[6] >= 0.0)) {
    } else {
        printf("Error: PreferredDirection_FromSensoryPts requires all sensory pt values be positive:\n");
        printf("       CondVarX[0] = %f\n",CondVarX[0]);
        printf("       CondVarX[0] = %f\n",CondVarX[1]);
        printf("       CondVarX[0] = %f\n",CondVarX[2]);
        printf("       CondVarX[0] = %f\n",CondVarX[3]);
        printf("       CondVarX[0] = %f\n",CondVarX[4]);
        printf("       CondVarX[0] = %f\n",CondVarX[5]);
        printf("       CondVarX[0] = %f\n",CondVarX[6]);
        exit(-1);
    }
    // output :: SVaoSVXY, SVaoMshXYZ
    PreferredDirection_FromSensoryPts(sp->sensor,CondVarX,DEBUG,&sp->SVaoSVXY,sp->SVaoMshXYZ);
    if (sp->SVaoSVXY < -180.0) {printf("Error in Subroutine PreferredDirection_FromGradient: *SVaoSVXY < -180.0\n"); exit(-1);}
    if (sp->SVaoSVXY >  180.0) {printf("Error in Subroutine PreferredDirection_FromGradient: *SVaoSVXY >  180.0\n"); exit(-1);}
    
    return;
}

//*************************************************************************
//*************************************************************************
//****                                                                 ****
//****                 Subroutine UnstrMeshWallEffect                  ****
//****                                                                 ****
//****                 developed by R. Andrew Goodwin                  ****
//****                                                                 ****
//*************************************************************************
//*************************************************************************                            Equivalence in BehaviorRule.f90:

void UnstrMeshWallEffect(SENSOR_SPECIES *sp, double dt, int DEBUG) {
    
    if (sp->wallEffect.activate_wall_effect == FALSE) return;
    
    int FSO, XYZ;
    
    // Calculate length of each sensory point ray shoot:
    double SPLength[4][NSENSORPT_MAX], SPDistThres[NSENSORPT_MAX];
    SPLength[0][0] = 0.0;
    SPLength[1][0] = 0.0;
    SPLength[2][0] = 0.0;
    for (FSO=1; FSO<NSENSORPT_MAX; FSO++) {
        SPLength[0][FSO] = sp->sensor[FSO].r.x - sp->sensor[0].r.x;
        SPLength[1][FSO] = sp->sensor[FSO].r.y - sp->sensor[0].r.y;
        SPLength[2][FSO] = sp->sensor[FSO].r.z - sp->sensor[0].r.z;
        SPLength[3][FSO] = sqrt(SPLength[0][FSO] * SPLength[0][FSO] + SPLength[1][FSO] * SPLength[1][FSO] + SPLength[2][FSO] * SPLength[2][FSO] );
        SPDistThres[FSO] = sp->sensor[FSO].distance * sp->wallEffect.wall_perm_threshold; //Coefficients(1,6);
    }
    
    // Find sensory point ray shoots coming up short, indicating they're hitting an impermeable wall:
    //  Override volitional movement angle (vertical):
    double DistToTravel_DueToFlow=0, DistToTravel_Allowable=0,DistToTravel_DueToFlowAndBehavior=0;
    if (sp->sensor[5].found == 1 && sp->sensor[6].found == 1) {
        // Vertical sensory points in-bounds (i.e., this is not 2-D depth-averaged mesh data) // Positive ==> up ; negative ==> down
        DistToTravel_DueToFlow = sp->sensor[0].CondVelM * dt * sin( sp->sensor[0].FVaoMshXYZ[1] * deg2rad );
        
        if ((SPLength[3][5] <  SPDistThres[5]) && (SPLength[3][6] >= SPDistThres[6])) {
            
            DistToTravel_Allowable = SPLength[3][5] - sp->wallEffect.wall_min_vdist; //Coefficients(1,7) // Positive ==> up ; negative ==> down (already too close)
            while (1) {
                DistToTravel_DueToFlowAndBehavior = DistToTravel_DueToFlow + sp->IdSpdRes * dt * sin( sp->SVaoMshXYZ[1] * deg2rad); // Positive ==> up ; negative ==> down
                if ((DistToTravel_Allowable > 0.0) && (DistToTravel_Allowable > DistToTravel_DueToFlowAndBehavior)) return; // No action needed, room to move above beyond what individual is going to do via flow/behavior
                if ((DistToTravel_Allowable < 0.0) && (DistToTravel_Allowable > DistToTravel_DueToFlowAndBehavior)) return; // No action needed, individual needs to go down but already sufficiently doing so via flow/behavior
                
                // Action needed, individual either headed up too much or individual not going down enough via flow/behavior:
                if ((sp->SVaoMshXYZ[1] - sp->wallEffect.wall_vert_angle) <= -90.0) break;
                if  ((sp->IdSpdRes + sp->sensor[0].CondVelM * sp->wallEffect.wall_speed_increase) >= sp->speed.length*sp->speed.burst) break;
                
                if ((sp->SVaoMshXYZ[1] - sp->wallEffect.wall_vert_angle) > -90.0) {
                    sp->SVaoMshXYZ[1] -= sp->wallEffect.wall_vert_angle; //Coefficients(1,8)
                } else if ((sp->IdSpdRes + sp->sensor[0].CondVelM * sp->wallEffect.wall_speed_increase) < sp->speed.length*sp->speed.burst) {
                    sp->IdSpdRes += sp->sensor[0].CondVelM * sp->wallEffect.wall_speed_increase; //Coefficients(1,9)
                }
            }
            
            sp->initialize_flag[7] = 1;
            
        } else if  ((SPLength[3][5] >= SPDistThres[5])  && (SPLength[3][6] < SPDistThres[6])) {
            
            DistToTravel_Allowable = -1.0 * (SPLength[3][6] - sp->wallEffect.wall_min_vdist); //Coefficients(1,7))  Positive ==> up (already too close) ; negative ==> down
            while (1) {
                DistToTravel_DueToFlowAndBehavior = DistToTravel_DueToFlow + sp->IdSpdRes * dt * sin( sp->SVaoMshXYZ[1] * deg2rad); // Positive ==> up ; negative ==> down
                if ((DistToTravel_Allowable < 0.0) && (DistToTravel_Allowable < DistToTravel_DueToFlowAndBehavior)) return; // No action needed, room to move below beyond what individual is going to do via flow/behavior
                if ((DistToTravel_Allowable > 0.0) && (DistToTravel_Allowable < DistToTravel_DueToFlowAndBehavior)) return; // No action needed, individual needs to go up but already sufficiently doing so via flow/behavior
                
                // Action needed, individual either headed down too much or individual not going up enough via flow/behavior:
                if ((sp->SVaoMshXYZ[1] - sp->wallEffect.wall_vert_angle) >= 90.0) break;
                if ((sp->IdSpdRes + sp->sensor[0].CondVelM * sp->wallEffect.wall_speed_increase) >= sp->speed.length*sp->speed.burst) break;
                
                if ((sp->SVaoMshXYZ[1] - sp->wallEffect.wall_vert_angle) < 90.0) {
                    sp->SVaoMshXYZ[1] += sp->wallEffect.wall_vert_angle; //Coefficients(1,8)
                } else if  ((sp->IdSpdRes + sp->sensor[0].CondVelM * sp->wallEffect.wall_speed_increase) < sp->speed.length*sp->speed.burst) {
                    sp->IdSpdRes += sp->sensor[0].CondVelM * sp->wallEffect.wall_speed_increase; //Coefficients(1,9)
                }
            }
            sp->initialize_flag[7] = 1;
            
        } else if  ((SPLength[3][5] < SPDistThres[5])  && (SPLength[3][6] < SPDistThres[6])) {
            
            sp->SVaoMshXYZ[1] = sp->sensor[0].FVaoMshXYZ[1];
            if ((sp->SVaoMshXYZ[1] >  90.0) && (sp->SVaoMshXYZ[1] <  90.05)) sp->SVaoMshXYZ[1] =  90.0; // sp->s[0].FVaoMshXYZ(2,1) can leak to  90.0001, which triggers an error in sp->SVaoMshXYZ[1]
            if ((sp->SVaoMshXYZ[1] < -90.0) && (sp->SVaoMshXYZ[1] > -90.05)) sp->SVaoMshXYZ[1] = -90.0; // sp->s[0].FVaoMshXYZ(2,1) can leak to -90.0001, which triggers an error in sp->SVaoMshXYZ[1]
            sp->initialize_flag[7] = 1;
            
        } else if  ((SPLength[3][5] >= SPDistThres[5])  && (SPLength[3][6] >= SPDistThres[6])) {
            //Do Nothing.
        }
    }
    return;
}

//*************************************************************************
//*************************************************************************
//****                                                                 ****
//****                 Subroutine VolMovInOtherFrames                  ****
//****                                                                 ****
//****                 developed by R. Andrew Goodwin                  ****
//****                                                                 ****
//*************************************************************************
//*************************************************************************
// output :: VldSVOrientXYZ, SVvoMshXYZ, SVvoFVXYZ, SVvoSVXYZ, sp->SVaoFVXY
void VolMovInOtherFrames(SENSOR_SPECIES *sp, int DEBUG) {
    
    double SVaoMshXYRad,SVaoMshZRad,IdSpdResXY,SVaoSVXYRad,SVaoMshXYTemp,SVaoFVXYRad;
    
    //Convert the volitional movement angle/velocity to vectors.
    SVaoMshZRad = sp->SVaoMshXYZ[1] * deg2rad;
    IdSpdResXY  = sp->IdSpdRes * cos(SVaoMshZRad);
    if ((sp->SVaoSVXY >=  0.0) && (sp->SVaoSVXY <= 90.0)) {
        SVaoSVXYRad     =  sp->SVaoSVXY  * deg2rad;
        sp->SVvoSVXYZ.x =  IdSpdResXY * cos(SVaoSVXYRad); // (S)wim (V)ector (v)elocity (o)ff the previous (S)wim (V)ector's xy-axis
        sp->SVvoSVXYZ.y =  IdSpdResXY * sin(SVaoSVXYRad); // (S)wim (V)ector (v)elocity (o)ff the previous (S)wim (V)ector's xy-axis
    } else if ((sp->SVaoSVXY >=  90.0) && (sp->SVaoSVXY <= 180.0)) {
        SVaoSVXYRad     = (180.0 - sp->SVaoSVXY) * deg2rad;
        sp->SVvoSVXYZ.x = -IdSpdResXY * cos(SVaoSVXYRad);
        sp->SVvoSVXYZ.y =  IdSpdResXY * sin(SVaoSVXYRad);
    } else if ((sp->SVaoSVXY <=   0.0) && (sp->SVaoSVXY >= -90.0)) {
        SVaoSVXYRad     =  -sp->SVaoSVXY  * deg2rad;
        sp->SVvoSVXYZ.x =  IdSpdResXY * cos(SVaoSVXYRad);
        sp->SVvoSVXYZ.y = -IdSpdResXY * sin(SVaoSVXYRad);
    } else if ((sp->SVaoSVXY <=  -90.0) && (sp->SVaoSVXY >= -180.0)) {
        SVaoSVXYRad     = (180.0 + sp->SVaoSVXY) * deg2rad;
        sp->SVvoSVXYZ.x = -IdSpdResXY * cos(SVaoSVXYRad);
        sp->SVvoSVXYZ.y = -IdSpdResXY * sin(SVaoSVXYRad);
    } else {
        printf("Error: Invalid sp->SVaoSVXY value\n");
        exit(-1);
    }
    
    
    //Determine angle of the volitional movement vector relative to:
    //  Mesh xy-axis:
    sp->SVaoMshXYZ[0] = sp->SVaoMshXYZ_old[0] + sp->SVaoSVXY; // (S)wim (V)ector (a)ngle (o)ff Mesh xy-axis orientation
    if (sp->SVaoMshXYZ[0] < 0.0) {
        sp->SVaoMshXYZ[0] = sp->SVaoMshXYZ[0] + 360.0;
    } else if (sp->SVaoMshXYZ[0] > 360.0) {
        sp->SVaoMshXYZ[0] = sp->SVaoMshXYZ[0] - 360.0;
    }
    SVaoMshXYTemp = sp->SVaoMshXYZ[0];
    if (SVaoMshXYTemp > 180.0) { // Ensure that angles go from 0.0 to 180.0 and 0.0 to -180.0
        SVaoMshXYTemp = SVaoMshXYTemp - 360.0;
    }
    if ((SVaoMshXYTemp >=  0.0) && (SVaoMshXYTemp <= 90.0)) {
        SVaoMshXYRad       = SVaoMshXYTemp  * deg2rad;
        sp->SVvoMshXYZ.x = IdSpdResXY * cos(SVaoMshXYRad); // (S)wim (V)ector (v)elocity (o)ff Mesh xy-axis orientation
        sp->SVvoMshXYZ.y = IdSpdResXY * sin(SVaoMshXYRad); // (S)wim (V)ector (v)elocity (o)ff Mesh xy-axis orientation
    } else if ((SVaoMshXYTemp >=  90.0) && (SVaoMshXYTemp <= 180.0)) {
        SVaoMshXYRad     = (180.0 - SVaoMshXYTemp) * deg2rad;
        sp->SVvoMshXYZ.x = -IdSpdResXY * cos(SVaoMshXYRad);
        sp->SVvoMshXYZ.y =  IdSpdResXY * sin(SVaoMshXYRad);
    } else if ((SVaoMshXYTemp <=  0.0) && (SVaoMshXYTemp >= -90.0)) {
        SVaoMshXYRad       = -SVaoMshXYTemp  * deg2rad;
        sp->SVvoMshXYZ.x = IdSpdResXY * cos(SVaoMshXYRad);
        sp->SVvoMshXYZ.y = -IdSpdResXY * sin(SVaoMshXYRad);
    } else if ((SVaoMshXYTemp <=  -90.0) && (SVaoMshXYTemp >= -180.0)) {
        SVaoMshXYRad       = (180.0 + SVaoMshXYTemp) * deg2rad;
        sp->SVvoMshXYZ.x = -IdSpdResXY * cos(SVaoMshXYRad);
        sp->SVvoMshXYZ.y = -IdSpdResXY * sin(SVaoMshXYRad);
    } else {
        printf("Error: Invalid sp->SVaoMshXYZ[0] value\n");
        exit(-1);
    }
    sp->SVvoMshXYZ.z = sp->IdSpdRes * sin(SVaoMshZRad);
    if ((sp->SVvoMshXYZ.x == 0.0) && (sp->SVvoMshXYZ.y == 0.0)) {
        sp->VldSVOrientXYZ[0] = 0;  // No valid xy-plane orientation of the individual's axis => no individual movement
    } else {
        sp->VldSVOrientXYZ[0] = 1;  // Valid xy-plane orientation of the individual's axis can be determined
    }
    
    //  Water flow vector xy-axis:
    if ((sp->VldSVOrientXYZ[0] == 1) && (sp->VldFVOrientXYZ[0] == 1)) {
        sp->SVaoFVXY = sp->SVaoMshXYZ[0] - sp->sensor[0].FVaoMshXYZ[0] ; // (S)wim (V)ector (a)ngle (o)ff (F)low (V)ector's xy-axis
        if (sp->SVaoFVXY < 0.0) {
            sp->SVaoFVXY = sp->SVaoFVXY + 360.0;
        } else if (sp->SVaoFVXY > 360.0) {
            sp->SVaoFVXY = sp->SVaoFVXY - 360.0;
        }
        if (sp->SVaoFVXY > 180.0) sp->SVaoFVXY = sp->SVaoFVXY - 360.0; // Make it so that angles go from 0.0 to 180.0 and 0.0 to -180.0
    } else {
        sp->SVaoFVXY = 0.0; // (S)wim (V)ector (a)ngle (o)ff (F)low (V)ector's xy-axis
    }
    
    //Determine velocity of the volitional movement vector relative to:
    //  Water flow vector xyz-axis velocity:
    if ((sp->SVaoFVXY >=  0.0) && (sp->SVaoFVXY <= 90.0)) {
        SVaoFVXYRad     =  sp->SVaoFVXY  * deg2rad;
        sp->SVvoFVXYZ.x =  IdSpdResXY * cos(SVaoFVXYRad); // (S)wim (V)ector (v)elocity (o)ff (F)low (V)ector's xy-axis
        sp->SVvoFVXYZ.y =  IdSpdResXY * sin(SVaoFVXYRad); // (S)wim (V)ector (v)elocity (o)ff (F)low (V)ector's xy-axis
    } else if ((sp->SVaoFVXY >=  90.0) && (sp->SVaoFVXY <= 180.0)) {
        SVaoFVXYRad     = (180.0 - sp->SVaoFVXY) * deg2rad;
        sp->SVvoFVXYZ.x = -IdSpdResXY * cos(SVaoFVXYRad);
        sp->SVvoFVXYZ.y =  IdSpdResXY * sin(SVaoFVXYRad);
    } else if ((sp->SVaoFVXY <=   0.0) && (sp->SVaoFVXY >= -90.0)) {
        SVaoFVXYRad     = -sp->SVaoFVXY  * deg2rad;
        sp->SVvoFVXYZ.x =  IdSpdResXY * cos(SVaoFVXYRad);
        sp->SVvoFVXYZ.y = -IdSpdResXY * sin(SVaoFVXYRad);
    } else if ((sp->SVaoFVXY <=  -90.0) && (sp->SVaoFVXY >= -180.0)) {
        SVaoFVXYRad     = (180.0 + sp->SVaoFVXY) * deg2rad;
        sp->SVvoFVXYZ.x = -IdSpdResXY * cos(SVaoFVXYRad);
        sp->SVvoFVXYZ.y = -IdSpdResXY * sin(SVaoFVXYRad);
    } else {
        printf("Error: Invalid sp->SVaoFVXY value\n");
        exit(-1);
    }
    sp->SVvoFVXYZ.x = sp->SVvoMshXYZ.x - sp->sensor[0].v.z; //IndvSensoryVelocity_NP(3,1);
    
    return;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void stimuli_perception_init(STIMULI_PERCEPTION *sp) {
    int i;
    for (i=0; i<NBEHAVS_MAX; i++) {
        sp->Evid_Threshold[i] = 0.0;
        sp->Evid_AcclimRate[i] = 0.0;
    }
    sp->AccIM_min = 0.0;
    sp->vert_increase = 0.0;
    sp->flow_drivenl = 0.0;
    sp->decay_rate = 0.0;
    sp->velM = 0.0;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void pt_interact_init(PT_INTERACT *pt) {
    pt->is_on = false;
    pt->percept_range = 0.0;
    pt->track_neighbor_dist = false;
    pt->ntrack = 0;
    pt->attention_threshold = 0.0;
    pt->max_angle = 0.0;
    pt->cue_threshold = 0.0;
    pt->target_value = 0.0;
    pt->segment_length = 0.0;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void wall_effect_init(WALL_EFFECT *we) {
    we->sp_mean_dist = 0.0;
    we->sp_stdev_dist = 0.0;
    we->downscale = FALSE;
    we->sp_min_dist = 0.0;
    we->activate_wall_effect = FALSE;
    we->wall_perm_threshold = 0.0;
    we->wall_min_vdist = 0.0;
    we->wall_vert_angle = 0.0;
    we->wall_speed_increase = 0.0;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void behavoir_trans_usher_init(BEHAVOIR_TRANS_USHER *be) {
    
    int i;
    for (i=0; i<NBEHAVS_MAX; i++) {
        be->Evid_Decay[i] = 0.0;
        be->Evid_Inhibit[i] = 0.0;
        be->white_noise_weight[i] = 0.0;
    }
    be->evidence_threshold = 0.0;
    be->Evid_Accumulated_t0 = 0.0;
    
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void behavoir_trans_anderson_init(BEHAVOIR_TRANS_ANDERSON *be) {
    
    int i;
    for (i=0; i<NBEHAVS_MAX; i++) {
        be->Evid_Accumulator[i] = 0.0;
        be->value_reward[i] = 0.0;
    }
    be->beh4_threshold = 0.0;
    be->prob_reward = 0.0;
    
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void sspeed_init(SSPEED *sp) {
    
    int i;
    for (i=0; i<NBEHAVS_MAX; i++) {
        sp->speed_coeff1[i] = 0.0;
        sp->value_reward[i] = 0.0;
    }
    sp->length = 0.0;
    sp->drift = 0.0;
    sp->cruise = 0.0;
    sp->sustained = 0.0;
    sp->burst = 0.0;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void sensing_orientation_ornstein_init(SENSING_ORIENTATION_ORNSTEIN *so) {
    int i;
    for (i=0; i<NBEHAVS_MAX; i++) {
        so->lambda_xy[i] = 0.0;
        so->lambda_z[i] = 0.0;
        so->c_xy[i] = 0.0;
        so->c_z[i] = 0.0;
    }
    for (i=0; i<4; i++) {
        so->slack_lambda_xy[i] = 0.0;
        so->slack_lambda_z[i] = 0.0;
        so->slack_c_xy[i] = 0.0;
        so->slack_c_z[i] = 0.0;
    }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void sensing_orientation_codling_init(SENSING_ORIENTATION_CODLING *so) {
    int i;
    for (i=0; i<NBEHAVS_MAX; i++) {
        so->kappa_xy[i] = 0.0;
        so->kappa_z[i] = 0.0;
        so->d_xy[i] = 0.0;
        so->d_z[i] = 0.0;
    }
    for (i=0; i<4; i++) {
        so->slack_kappa_xy[i] = 0.0;
        so->slack_kappa_z[i] = 0.0;
        so->slack_d_xy[i] = 0.0;
        so->slack_d_z[i] = 0.0;
    }
}

//*************************************************************************
//*************************************************************************
//****                                                                 ****
//****                   Subroutine VectorRelation                     ****
//****                                                                 ****
//****                 developed by R. Andrew Goodwin                  ****
//****                                                                 ****
//*************************************************************************
//*************************************************************************

//     I                          VectorType,                              ! Local scalar              => PositionVector=0 ; VelocityVector=1
//     I                          MVComponentsMsh,                         ! Local array              MVComponentsMsh(1:3)
//     I                          VldRVOrientXYZ,                          ! Local array              VldRVOrientXYZ(1:2)
//     O                          VldMVOrientXYZ,                          ! Local array              VldMVOrientXYZ(1:2)
//     I                          RVaoMshXYZ,                              ! Local array              RVaoMshXYZ(1:2)
//     I                          RVvoMshXYZ,                              ! Local array              RVvoMshXYZ(1:3)
//     O                          MVaoRVXY,                                ! Local scalar             => (M)aster (V)ector (a)ngle (o)ff (R)eference (V)ector's xy-axis
//     O                          MVaoMshXYZ,                              ! Local array              MVaoMshXYZ(1:2)
//     O                          MVvoRVXYZ                                ! Local array              MVvoRVXYZ(1:3) => (M)aster (V)ector (v)elocity (o)ff (R)eference (V)ector's xy-axis

void VectorRelation(int VectorType, SVECT MVComponentsMsh, bool *VldRVOrientXYZ, bool *VldMVOrientXYZ, double *RVaoMshXYZ, SVECT RVvoMshXYZ, double *MVaoRVXY, double *MVaoMshXYZ, SVECT *MVvoRVXYZ) {

    double AbsMVaoMshXY=0.,MVvroMshXY=0.,AbsMVaoMshZ=0.,MVvroRVXY=0.,MVaoRVXYRad=0.;
    SVECT MVvoRVMshXYZ;
    
    // Determine angle of the Master Vector relative to the
    // Mesh xy-axis:
    
    if ((MVComponentsMsh.x == 0.0) && (MVComponentsMsh.y == 0.0)) {
        VldMVOrientXYZ[0] = FALSE; // No valid xy-plane orientation of master vector
        MVaoMshXYZ[0] = 0.0;
    } else {
        VldMVOrientXYZ[0] = TRUE; // Valid xy-plane orientation of master vector can be determined
        AbsMVaoMshXY = fabs(atan2(fabs(MVComponentsMsh.y), fabs(MVComponentsMsh.x)));
        if ((MVComponentsMsh.x >= 0.0) && (MVComponentsMsh.y >= 0.0)) {
            MVaoMshXYZ[0] = AbsMVaoMshXY * rad2deg;      // (M)aster (V)ector (a)ngle (o)ff Mesh xy-axis orientation
        } else  if ((MVComponentsMsh.x <= 0.0) && (MVComponentsMsh.y >= 0.0)) {
            MVaoMshXYZ[0] = 180.0 - AbsMVaoMshXY * rad2deg;      // (M)aster (V)ector (a)ngle (o)ff Mesh xy-axis orientation
        } else  if ((MVComponentsMsh.x <= 0.0) && (MVComponentsMsh.y <= 0.0)) {
            MVaoMshXYZ[0] = 180.0 + AbsMVaoMshXY * rad2deg;      // (M)aster (V)ector (a)ngle (o)ff Mesh xy-axis orientation
        } else  if ((MVComponentsMsh.x >= 0.0) && (MVComponentsMsh.y <= 0.0)) {
            MVaoMshXYZ[0] = 360.0 - AbsMVaoMshXY * rad2deg;      // (M)aster (V)ector (a)ngle (o)ff Mesh xy-axis orientation
        }
    }
    
    // Determine angle of the Master Vector vertical component relative to the
    // Mesh xy-plane:
    
    MVvroMshXY = sqrt( MVComponentsMsh.x * MVComponentsMsh.x + MVComponentsMsh.y * MVComponentsMsh.y );
    if ((MVvroMshXY == 0.0) && (MVComponentsMsh.z == 0.0)) {
        VldMVOrientXYZ[1] = FALSE; // No valid vertical orientation of master vector
        MVaoMshXYZ[1] = 0.0;
    } else {
        VldMVOrientXYZ[1] = TRUE;// Valid vertical orientation of master vector can be determined
        AbsMVaoMshZ = fabs(atan2(fabs(MVComponentsMsh.z),fabs(MVvroMshXY)));
        if (MVComponentsMsh.z >= 0.0) {
            MVaoMshXYZ[1] =  AbsMVaoMshZ * rad2deg; // (M)aster (V)ector (a)ngle (o)ff Mesh xy-plane
        } else {
            MVaoMshXYZ[1] = -AbsMVaoMshZ * rad2deg; // (M)aster (V)ector (a)ngle (o)ff Mesh xy-plane
        }
    }
    
    // Determine angle of the Master Vector relative to the
    // Reference Vector's xy-axis:
    
    if ((VldRVOrientXYZ[0]) && (VldMVOrientXYZ[0])) {
        *MVaoRVXY = MVaoMshXYZ[0] - RVaoMshXYZ[0]; // (M)aster (V)ector (a)ngle (o)ff (R)eference (V)ector's xy-axis
        if  (*MVaoRVXY < 0.0) {
            *MVaoRVXY = *MVaoRVXY + 360.0;
        } else if (*MVaoRVXY > 360.0) {
            *MVaoRVXY = *MVaoRVXY - 360.0;
        }
        if (*MVaoRVXY > 180.0) {*MVaoRVXY = *MVaoRVXY - 360.0;} // Make angles go from 0.0 to 180.0 and 0.0 to -180.0
    } else {
        *MVaoRVXY = 0.0;// (M)aster (V)ector (a)ngle (o)ff (R)eference (V)ector's xy-axis
    }
    
    // Determine velocity of the Master Vector relative to the
    // Reference Vector's xyz-axis velocity:
    
    if (VectorType == 0) { // Skip velocity calcs if vector is a 'position vector' (and not a 'velocity vector')
        svect_init(MVvoRVXYZ);
        return;
    }
    
    MVvoRVMshXYZ.x = MVComponentsMsh.x - RVvoMshXYZ.x; // RVvoMshXYZ(XYZ) units are m/s, so calcs here only valid if MVComponentsMsh(XYZ) units are also m/s
    MVvoRVMshXYZ.y = MVComponentsMsh.y - RVvoMshXYZ.y; // RVvoMshXYZ(XYZ) units are m/s, so calcs here only valid if MVComponentsMsh(XYZ) units are also m/s
    MVvoRVMshXYZ.z = MVComponentsMsh.z - RVvoMshXYZ.z; // RVvoMshXYZ(XYZ) units are m/s, so calcs here only valid if MVComponentsMsh(XYZ) units are also m/s
    
    // (M)aster (V)ector (v)elocity (r)esultant (o)ff (R)eference (V)ector's xy-axis
    MVvroRVXY = sqrt( MVvoRVMshXYZ.x * MVvoRVMshXYZ.x + MVvoRVMshXYZ.y * MVvoRVMshXYZ.y );
    if      ((*MVaoRVXY >=  0.0) && (*MVaoRVXY <= 90.0)) {
        MVaoRVXYRad  =          *MVaoRVXY  * deg2rad;
        MVvoRVXYZ->x =  MVvroRVXY * cos(MVaoRVXYRad); // (M)aster (V)ector (v)elocity (o)ff (R)eference (V)ector's xy-axis
        MVvoRVXYZ->y =  MVvroRVXY * sin(MVaoRVXYRad); // (M)aster (V)ector (v)elocity (o)ff (R)eference (V)ector's xy-axis
    } else  if ((*MVaoRVXY >=  90.0) && (*MVaoRVXY <= 180.0)) {
        MVaoRVXYRad  = (180.0 - *MVaoRVXY) * deg2rad;
        MVvoRVXYZ->x = -MVvroRVXY * cos(MVaoRVXYRad);
        MVvoRVXYZ->y =  MVvroRVXY * sin(MVaoRVXYRad);
    } else  if ((*MVaoRVXY <=   0.0) && (*MVaoRVXY >= -90.0)) {
        MVaoRVXYRad  =         -*MVaoRVXY  * deg2rad;
        MVvoRVXYZ->x =  MVvroRVXY * cos(MVaoRVXYRad);
        MVvoRVXYZ->y = -MVvroRVXY * sin(MVaoRVXYRad);
    } else  if ((*MVaoRVXY <=  -90.0) && (*MVaoRVXY >= -180.0)) {
        MVaoRVXYRad  = (180.0 + *MVaoRVXY) * deg2rad;
        MVvoRVXYZ->x = -MVvroRVXY * cos(MVaoRVXYRad);
        MVvoRVXYZ->y = -MVvroRVXY * sin(MVaoRVXYRad);
    } else {
        printf("Error: Invalid MVaoRVXY value\n");
        exit(-1);
    }
    MVvoRVXYZ->z = MVvoRVMshXYZ.z;
    return;
}

//*************************************************************************
//*************************************************************************
//****                                                                 ****
//****          Subroutine PreferredDirection_FromSensoryPts           ****
//****                                                                 ****
//****                 developed by R. Andrew Goodwin                  ****
//****                                                                 ****
//*************************************************************************
//*************************************************************************

//  SPFound_NP,       &   // i Local array              SPFound_NP(1:FSOLIMIT)
//  BVvaSP,           &   // i Local array              BVvaSP(1:FSOLIMIT)
//  *SVaoSVXY,      &   // o Local scalar
//  SVaoMshXYZ_NP,    &   // o Local array              SVaoMshXYZ_NP(1:2)

void PreferredDirection_FromSensoryPts(SENSOR_PT *sp, double *BVvaSP, bool DEBUG, double *SVaoSVXY, double *SVaoMshXYZ) {
    
    //      Integer SPFound_NP(FSOLIMIT),                                  &   // Array
    //              DF,FSO,BstBVvSPXY,BstBVvSPZ,BstBVvSPXYCnt,BstBVvSPXYA, &   // Scalars
    //              BstBVvSPXYB,BstBVvSPXY2,BstBVvSPXY3,BstBVvSPXY4
    //
    //      Real    BVvaSP(FSOLIMIT),                                      &   // Arrays
    //              SVaoMshXYZ[1],                                      &
    //              SenPaoSV[4],                                           &
    //              *SVaoSVXY,                                           &   // Scalars
    //              BstBVvXY,AbsSVaoMshZ,SUMConditions,TstBVvXY
    
    
    int FSO;
    
    //Sensory point angles relative to the direction pointing from tail to head.
    double SenPaoSV[5];
    SenPaoSV[0] =   0.0; // Not applicable
    SenPaoSV[1] =   0.0; // Front
    SenPaoSV[2] = 180.0; // Back
    SenPaoSV[3] =  90.0; // Left side
    SenPaoSV[4] = -90.0; // Right side
    
    
    //Identify cardinal positions with minimum BVvaSP among in-bounds sensory points.
    // XY-plane:
    
    int BstBVvSPXY    = 2; // (B)e(st) (B)ehavior (V)ariable (v)alue (S)ensory (P)oint in (xy)-plane
    // Forward sensory point always in-bounds and the default direction in gradient-tracking
    double BstBVvXY   = BVvaSP[BstBVvSPXY-1]; // (B)e(st) (B)ehavior (V)ariable (v)alue in (xy)-plane
    int BstBVvSPXYCnt = 1; // # of sensory points with the same value of BstBVvXY
    int BstBVvSPXYA   = 4; // Sensory point counter-clockwise of BstBVvSPXY
    int BstBVvSPXYB   = 5; // Sensory point         clockwise of BstBVvSPXY
    double TstBVvXY   = 0.0; // (T)e(st) (B)ehavior (V)ariable (v)alue in (xy)-plane
    int BstBVvSPXY2   = 0; // Location of 2nd     sensory point with same value of BstBVvXY
    int BstBVvSPXY3   = 0; // Location of 3rd     sensory point with same value of BstBVvXY
    int BstBVvSPXY4   = 0; // Location of 4th (+) sensory point with same value of BstBVvXY
    
    for (FSO=2; FSO<5; FSO++ ) { // Look at sensory points 3-5 in the xy-plane
        if (sp[FSO].found == 1) { // Sensory point in-bounds?
            TstBVvXY = BVvaSP[FSO]; // (T)e(st) (B)ehavior (V)ariable (v)alue in (xy)-plane
            if (TstBVvXY > BstBVvXY) {                               // Change "<" to ">" if higher values desired
                BstBVvSPXYCnt = 1;
                BstBVvXY      = TstBVvXY;
                BstBVvSPXY    = FSO;
                BstBVvSPXY2   = 0;
                BstBVvSPXY3   = 0;
                BstBVvSPXY4   = 0;
                if (BstBVvSPXY == 3) { // if location of BstBVvXY is sensory point #3
                    BstBVvSPXYA = 5; //   ...sensory pt to  left (i.e., counter-clockwise)
                    BstBVvSPXYB = 4; //   ...sensory pt to right (i.e.,         clockwise)
                } else if (BstBVvSPXY == 4) { // if location of BstBVvXY is sensory point #4
                    BstBVvSPXYA = 3; //   ...sensory pt to  left (i.e., counter-clockwise)
                    BstBVvSPXYB = 2; //   ...sensory pt to right (i.e.,         clockwise)
                } else if (BstBVvSPXY == 5) { // if location of BstBVvXY is sensory point #5
                    BstBVvSPXYA = 2; //   ...sensory pt to  left (i.e., counter-clockwise)
                    BstBVvSPXYB = 3; //   ...sensory pt to right (i.e.,         clockwise)
                }
            } else if (TstBVvXY == BstBVvXY) {
                BstBVvSPXYCnt = BstBVvSPXYCnt + 1;
                if      (BstBVvSPXYCnt == 2) { // 2 sensory points have same value of BstBVvXY
                    BstBVvSPXY2 = FSO;
                } else if (BstBVvSPXYCnt == 3) { // 3 sensory points have same value of BstBVvXY
                    BstBVvSPXY3 = FSO;
                } else if (BstBVvSPXYCnt > 3) { // >3 sensory points have same value of BstBVvXY
                    BstBVvSPXY4 = FSO;
                } else {
                        printf("Error: Rules #1\n");
                    exit(-1);
                }
            }
        }
    }
    
    // Vertical:
    
    int BstBVvSPZ = 0; // (B)e(st) (B)ehavior (V)ariable (v)alue (S)ensory (P)oint in vertical direction
    if ((sp[5].found == 1) &&  (sp[6].found == 1)) { // Vertical sensory points in-bounds
        if  (BVvaSP[5] > BstBVvXY )  BstBVvSPZ = 6; // Change "<" to ">" if higher values desired
        if ((BVvaSP[6] > BstBVvXY ) && (BVvaSP[6] > BVvaSP[5])) BstBVvSPZ = 7; // Change "<" to ">" if higher values desired
        if  (BVvaSP[5] == BVvaSP[6])  BstBVvSPZ = 0;// In case sensory pts 6 & 7 are best but equal
    }
    
    
    // Volitional angle.
    // XY-plane:
    //   Special conditions:
    while (1) {
        if  ((sp[3].found == 0) &&  (sp[4].found == 0)) { // Left/right sensory points out-of-bounds in 2-D laterally-averaged mesh data
            if (BstBVvSPXY == 2) *SVaoSVXY = 0.0;
            if (BstBVvSPXY == 3) *SVaoSVXY = 180.0;
            break;  // Skip to vertical volitional movement direction calculations
        } else if ((BstBVvSPXY2 != 0) &&  (BstBVvSPXY3 == 0)) { // 2 xy-plane sensory points have same value of BstBVvXY
            if       (((BstBVvSPXY == 2) && (BstBVvSPXY2 == 3)) || ((BstBVvSPXY == 3) && (BstBVvSPXY2 == 2))) {
                *SVaoSVXY = 0.0;
            } else if (((BstBVvSPXY == 4) && (BstBVvSPXY2 == 5)) || ((BstBVvSPXY == 5) && (BstBVvSPXY2 == 4))) {
                *SVaoSVXY = 0.0;
            } else if (((BstBVvSPXY == 3) && (BstBVvSPXY2 == 5)) || ((BstBVvSPXY == 5) && (BstBVvSPXY2 == 3))) { // Don't want to try and average -90 and 180 for SVaoSVXY
                *SVaoSVXY = -135.0;
            } else {
                *SVaoSVXY = (SenPaoSV[BstBVvSPXY-1] +  SenPaoSV[BstBVvSPXY2-1]) / 2.0; // Average of the 2 directions having same BstBVvXY value
            }
            break; // Skip to vertical volitional movement direction calculations
        } else if ((BstBVvSPXY3 != 0) &&  (BstBVvSPXY4 == 0)) { // 3 xy-plane sensory points have same value of BstBVvXY
            if (((BstBVvSPXY  == 2)  &&  (BstBVvSPXY2 == 3) && (BstBVvSPXY3 == 5)) ||
                ((BstBVvSPXY  == 2)  &&  (BstBVvSPXY2 == 5) && (BstBVvSPXY3 == 3)) ||
                ((BstBVvSPXY  == 3)  &&  (BstBVvSPXY2 == 2) && (BstBVvSPXY3 == 5)) ||
                ((BstBVvSPXY  == 3)  &&  (BstBVvSPXY2 == 5) && (BstBVvSPXY3 == 2)) ||
                ((BstBVvSPXY  == 5)  &&  (BstBVvSPXY2 == 2) && (BstBVvSPXY3 == 3)) ||
                ((BstBVvSPXY  == 5)  &&  (BstBVvSPXY2 == 3) && (BstBVvSPXY3 == 2))) { // Don't want to try and average with +180 here for SVaoSVXY
                *SVaoSVXY = -90.0; // Average of the 3 directions yielding the same BstBVvXY value
            } else {
                *SVaoSVXY = (SenPaoSV[BstBVvSPXY-1] + SenPaoSV[BstBVvSPXY2-1] + SenPaoSV[BstBVvSPXY3-1]) / 3.0; // Average of the 3 directions having same BstBVvXY value
            }
            break; // Skip to vertical volitional movement direction calculations
        } else if (BstBVvSPXY4 != 0) {                                  // 4 (all) xy-plane sensory points have same value of BstBVvXY
            *SVaoSVXY = 0.0;
            break; // Skip to vertical volitional movement direction calculations
        }
        
        //   Normal conditions:
        //     SenPaoSV[1] = xy-plane angle for individual's sensory point #2 ~=   0.0  (     Front)
        //     SenPaoSV[2] = xy-plane angle for individual's sensory point #3 ~= 180.0  (      Back)
        //     SenPaoSV[3] = xy-plane angle for individual's sensory point #4 ~=  90.0  ( Left side)
        //     SenPaoSV[4] = xy-plane angle for individual's sensory point #5 ~= -90.0  (Right side)
        
        double SUMConditions = (BVvaSP[BstBVvSPXYA-1]) + (BVvaSP[BstBVvSPXY-1]) + (BVvaSP[BstBVvSPXYB-1]);
        if  (BstBVvSPXY == 2) {
            *SVaoSVXY = (BVvaSP[BstBVvSPXYA-1]) * SenPaoSV[3] + (BVvaSP[BstBVvSPXY-1]) * SenPaoSV[1] + (BVvaSP[BstBVvSPXYB-1]) * SenPaoSV[4];
            *SVaoSVXY /= SUMConditions;
            if (SUMConditions >   1E9) *SVaoSVXY =   0.0;
            if (*SVaoSVXY   < -90.0)   *SVaoSVXY = -90.0;
            if (*SVaoSVXY   >  90.0)   *SVaoSVXY =  90.0;
        } else if (BstBVvSPXY == 3) {
            *SVaoSVXY = (BVvaSP[BstBVvSPXYA-1]) * (360.0 + SenPaoSV[4]) +  (BVvaSP[BstBVvSPXY-1]) * SenPaoSV[2]  +  (BVvaSP[BstBVvSPXYB-1]) * SenPaoSV[3];
            *SVaoSVXY /= SUMConditions;
            if (SUMConditions >   1E9) *SVaoSVXY = 180.0;
            if (*SVaoSVXY   <  90.0)   *SVaoSVXY =  90.0;
            if (*SVaoSVXY   > 270.0)   *SVaoSVXY = 270.0;
            if (*SVaoSVXY   > 180.0)   *SVaoSVXY -= 360.0;
        } else if (BstBVvSPXY == 4) {
            *SVaoSVXY = (BVvaSP[BstBVvSPXYA-1]) * SenPaoSV[2] + (BVvaSP[BstBVvSPXY-1]) * SenPaoSV[3] + (BVvaSP[BstBVvSPXYB-1]) * SenPaoSV[1];
            *SVaoSVXY /= SUMConditions;
            if (SUMConditions >   1E9) *SVaoSVXY =  90.0;
            if (*SVaoSVXY   <   0.0)   *SVaoSVXY =   0.0;
            if (*SVaoSVXY   > 180.0)   *SVaoSVXY = 180.0;
        } else if (BstBVvSPXY == 5) {
            *SVaoSVXY = (BVvaSP[BstBVvSPXYA-1]) *  SenPaoSV[1] + (BVvaSP[BstBVvSPXY-1]) *  SenPaoSV[4] +  (BVvaSP[BstBVvSPXYB-1]) * (SenPaoSV[2] - 360.0);
            *SVaoSVXY /= SUMConditions;
            if (SUMConditions >    1E9) *SVaoSVXY =  -90.0;
            if (*SVaoSVXY   < -180.0)   *SVaoSVXY = -180.0;
            if (*SVaoSVXY   >    0.0)   *SVaoSVXY =    0.0;
        }
    }
    
    // Vertical:
    if ((sp[5].found == 0) &&  (sp[6].found == 0)) {
        SVaoMshXYZ[1] = 0.0;
        if (SVaoMshXYZ[1] < -90.0) SVaoMshXYZ[1] = -90.0;
        if (SVaoMshXYZ[1] >  90.0) SVaoMshXYZ[1] =  90.0;
        return; // Skip to calculating volitional movement speed
    }
    if (BstBVvSPZ != 0) { // Vertical volitional movement vector response is dominant
        double AbsSVaoMshZ = fabs(atan2(fabs(BVvaSP[BstBVvSPZ-1]),fabs(BVvaSP[BstBVvSPXY-1])));
        SVaoMshXYZ[1] = AbsSVaoMshZ * rad2deg;
        if (BstBVvSPZ == 7) SVaoMshXYZ[1] = -SVaoMshXYZ[1];
    } else {                                                               // XY-plane volitional movement vector response is dominant
        SVaoMshXYZ[1] = 0.0;
    }
    if (SVaoMshXYZ[1] < -90.0) SVaoMshXYZ[1] = -90.0;
    if (SVaoMshXYZ[1] >  90.0) SVaoMshXYZ[1] =  90.0;
    
    return;
}
