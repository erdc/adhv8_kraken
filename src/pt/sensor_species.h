#ifndef H_SENSOR_SPECIES_
#define H_SENSOR_SPECIES_

#define NBEHAVS_MAX 8  // the maximum # of biological species behavoirs
#define NIDATTRB 1  // individual attribute :: ! Attribute #1: 0=Individual ; 1=Invertebrate ; -1=To be removed from simulation
#define NSENSORPT_MAX 7  // the maximum # of biological species sensor points

// dependencies :: SVECT

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

typedef struct {
    
    // Evidence threshold that increases WoE
    // Coefficients(6,1-8)
    double Evid_Threshold[NBEHAVS_MAX];
    
    // Rate of acclimatizing to evidence supporting Behavior
    // Coefficients(6,9-16)
    double Evid_AcclimRate[NBEHAVS_MAX];
    
    // Minimum lower-bound for initial acclimatization to AcclM at start-up => Probable Range: 0.4 to 1.0
    // Coefficients(6,17)
    double AccIM_min;
    
    // Beh4 vertical orientation increase when conditions not improving and not driven by flow vector => Range: 0.01 to 45.0
    // Coefficients(6,18)
    double vert_increase;
    
    // likely driven by flow vector => Range: 1.0 (same angle as inverse flow vector) to 2.0 (twice the inverse flow vector angle)
    // Coefficients(6,19)
    double flow_drivenl;
    
    // decay rate when conditions improving => Range: 0.8 (rather fast exponential decay) to 1.0 (no decay)
    // Coefficients(6,20)
    double decay_rate;
    
    // VelM delineating advective from slack flow, which modifies orientation => Range: ? to ?
    // Coefficients(6,21)
    double velM;
    
} STIMULI_PERCEPTION;

// struct methods ------------------------------

void stimuli_perception_init(STIMULI_PERCEPTION *sp);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

typedef struct {
    
    // Activate particle-to-particle interaction => 1=Yes ; 0=No
    // Coefficients(3,1)
    bool is_on;
    
    // Max perceptual range => Example values include (SPH): 500 ; (fish drift-feeding): 0.5x to 2.0x (Coefficients(4,17)*Coefficients(4,21)*TimeStepLength)
    // Coefficients(3,2)
    double percept_range;
    
    // Track nearest neighbor distances (NNDs) for all particles (0=Yes) or only particles with specific attributes (1=Yes; such as attribute #1 that can differentiate btw fish and invertebrate)
    // Coefficients(3,3)
    bool track_neighbor_dist;
    
    // Max # of particles to track NNDs => -1={All simulated individuals} ; >0={# of closest particles only}. This value cannot be zero (put 1 if Coefficients(3,1)=0) ::
    // Coefficients(3,4)
    int ntrack;
    
    // Beh6 threshold of perceived change in range for attention to be paid to a target
    // Coefficients(3,5)
    double attention_threshold;
    
    // max Abs(angle) relative to volitional movement vector for attention to target => Range: 0.0 to 180.0; Probable Range: 0.0 to 90.0
    // Coefficients(3,6)
    double max_angle;
    
    // threshold of consecutive segments w/attention required to cue on target when it begins moving away => Range: non-negative
    // Coefficients(3,7)
    double cue_threshold;
    
    // value of target => Range: 0.0 to 1.0
    // Coefficients(3,8)
    double target_value;
    
    // length (m) of segments the linear path of target movement (from t-1 to t) will be broken into for analysis => Probable Range: 0.01 to 0.5
    // Coefficients(3,9)
    double segment_length;
    
} PT_INTERACT;

// struct methods ------------------------------

void pt_interact_init(PT_INTERACT *pt);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

typedef struct {
    
    // Mean sensory ovoid point distance (m)
    // Coefficients(1,1)
    double sp_mean_dist;
    
    // StdDev of sensory ovoid point distance (m)
    // Coefficients(1,2)
    double sp_stdev_dist;
    
    // Activate downscaling of sensory ovoid point distance based on field variable value; if used, Coefficients(1,1-2) become the max values, and the min values need to be set manually in SensoryPtCreate => 1=Yes ; 0=No
    // Coefficients(1,3)
    bool downscale;
    
    // Min sensory ovoid point distance (m)
    // Coefficients(1,4)
    double sp_min_dist;
    
    // Activate wall effect => 1=Yes ; 0=No
    // Coefficients(1,5)
    bool activate_wall_effect;

    // Wall effect: threshold fraction of sensory point ray shoot below which identifies an impermeable wall => Range: 0.0 to 1.0
    // Coefficients(1,6)
    double wall_perm_threshold;
    
    // Wall effect: min vertical distance (m) from above/below impermeable wall to maintain => Range: 0.0 to 1.0
    // Coefficients(1,7)
    double wall_min_vdist;
    
    // Wall effect: incremental change in vertical angle when impermeable wall response needed => Range: 0.0 to 90.0; Probable Range: 5.0 to 10.0
    // Coefficients(1,8)
    double wall_vert_angle;
    
    // Wall effect: incremental speed increase (FhSpdRes_NP=FhSpdRes_NP+VelM*Coefficients(1,9)) for when vertical angle at max (+/-90.0) but volitional response still inadequate to maintain distance from impermeable wall => Range: >0.0
    // Coefficients(1,9)
    double wall_speed_increase;
    
} WALL_EFFECT;

// struct methods ------------------------------

void wall_effect_init(WALL_EFFECT *we);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

typedef struct {
    
    // Evidence threshold for Beh4 (Activates Behavior 4 so that it is independent of the xy-plane behavior trigger)
    // Coefficients(8,1)
    double evidence_threshold;
    
    // Evid_Accumulated (Value used for 1st time step)
    // Coefficients(8,2)
    double Evid_Accumulated_t0;
    
    // Evid_k for BehX - Decay rate of activity for evidence supporting Behaviors
    // Coefficients(8,3-10)
    double Evid_Decay[NBEHAVS_MAX];
    
    // Evid_InhibitWgt for BehX - Weight of mutual inhibition from other behaviors for Behavior
    // Coefficients(8,11-18)
    double Evid_Inhibit[NBEHAVS_MAX];
    
    // c_LCA for BehX (Weight of white noise)
    // Coefficients(8,19-26)
    double white_noise_weight[NBEHAVS_MAX];
    
    
} BEHAVOIR_TRANS_USHER;

// struct methods ------------------------------

void behavoir_trans_usher_init(BEHAVOIR_TRANS_USHER *be);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

typedef struct {
    
    // Utility threshold for Beh4 (Activates Behavior 4 so that it is independent of the xy-plane behavior trigger)
    // Coefficients(2,1)
    double beh4_threshold;
    
    // Prob_Reward (Value used for 1st time step). Util_Behav1 = Coefficients(2,2) if both (a) Valu_Reward1 = 1.0 and (b) Evid_Accumulator1 = 1.0
    // Coefficients(2,2)
    double prob_reward;
    
    // Evid_Accumulator for BehX (Evidence accumulator integrating past and present evidence)
    // Coefficients(2,3-10)
    double Evid_Accumulator[NBEHAVS_MAX];
    
    // Valu_Reward for BehX (Intrinsic value of Behavior 1 reward, if successful)
    // Coefficients(2,11-18)
    double value_reward[NBEHAVS_MAX];
    
} BEHAVOIR_TRANS_ANDERSON;

// struct methods ------------------------------

void behavoir_trans_anderson_init(BEHAVOIR_TRANS_ANDERSON *be);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

typedef struct {
    
    // Speed coeff #1 for BehX; scales default StdDev (IdBodyLngthVel(2)-IdBodyLngthVel(1)) around mean IdBodyLngthVel(2) => Range: ? to ?
    // Coefficients(4,1-8)
    double speed_coeff1[NBEHAVS_MAX];
    
    // Speed coeff #2 for BehX; scales default StdDev (IdBodyLngthVel(2)-IdBodyLngthVel(1)) around mean IdBodyLngthVel(2) => Range: ? to ?
    // Coefficients(4,9-16)
    double value_reward[NBEHAVS_MAX];
    
    // Individual Length (meters)
    // Coefficients(4,17)
    double length;
    
    // Individual volitional velocity #1 (body lengths/sec), e.g., drift
    // Coefficients(4,18)
    double drift;
    
    // Individual volitional velocity #2 (body lengths/sec), e.g., cruising
    // Coefficients(4,19)
    double cruise;
    
    // Individual volitional velocity #3 (body lengths/sec), e.g., sustained
    // Coefficients(4,20)
    double sustained;
    
    // Individual volitional velocity #4 (body lengths/sec), e.g., burst
    // Coefficients(4,21)
    double burst;
    
} SSPEED;

// struct methods ------------------------------

void sspeed_init(SSPEED *sp);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

typedef struct {
    
    //  lambda_xy for BehX (Sensing ability of the individual; 1=perfect)
    // Coefficients(7,1-8)
    double lambda_xy[NBEHAVS_MAX];
    
    // lambda_z for BehX (nice value to start with: 0.1)
    // Coefficients(7,9-16)
    double lambda_z[NBEHAVS_MAX];
    
    // c_xy for BehX (Orienting ability of the individual; 0=perfect)
    // Coefficients(7,17-24)
    double c_xy[NBEHAVS_MAX];
    
    // c_z for BehX (nice value to start with: 0.2)
    // Coefficients(7,25-32)
    double c_z[NBEHAVS_MAX];
    
    //  lambda_xy for BehX (Sensing ability of the individual; 1=perfect) - slack - (low water velocities)
    // Coefficients(7,33-36)
    double slack_lambda_xy[4];
    
    // lambda_z for BehX (nice value to start with: 0.1) - slack - (low water velocities)
    // Coefficients(7,37-40)
    double slack_lambda_z[4];
    
    // c_xy for BehX (Orienting ability of the individual; 0=perfect) - slack - (low water velocities)
    // Coefficients(7,41-44)
    double slack_c_xy[4];
    
    // c_z for BehX (nice value to start with: 0.2) - slack - (low water velocities)
    // Coefficients(7,45-48)
    double slack_c_z[4];
    
} SENSING_ORIENTATION_ORNSTEIN;

// struct methods ------------------------------

void sensing_orientation_ornstein_init(SENSING_ORIENTATION_ORNSTEIN *so);

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

typedef struct {
    
    // d_xy for BehX (Orienting ability of the individual; 0=perfect)
    // Coefficients(5,1-8)
    double d_xy[NBEHAVS_MAX];
    
    // d_z for BehX (nice value to start with: 0.2)
    // Coefficients(5,9-16)
    double d_z[NBEHAVS_MAX];
    
    //  kappa_xy for BehX (Sensing ability of the individual; 1=perfect)
    // Coefficients(5,17-24)
    double kappa_xy[NBEHAVS_MAX];
    
    // kappa_z for BehX (nice value to start with: 0.1)
    // Coefficients(5,25-32)
    double kappa_z[NBEHAVS_MAX];
    
    // d_xy for BehX (Orienting ability of the individual; 0=perfect) - slack - (low water velocities)
    // Coefficients(5,33-35)
    double slack_d_xy[4];
    
    // d_z for BehX (nice value to start with: 0.2) - slack - (low water velocities)
    // Coefficients(5,37-40)
    double slack_d_z[4];
    
    //  kappa_xy for BehX (Sensing ability of the individual; 1=perfect) - slack - (low water velocities)
    // Coefficients(5,41-44)
    double slack_kappa_xy[4];
    
    // kappa_z for BehX (nice value to start with: 0.1) - slack - (low water velocities)
    // Coefficients(5,45-48)
    double slack_kappa_z[4];
    
} SENSING_ORIENTATION_CODLING;

// struct methods ------------------------------

void sensing_orientation_codling_init(SENSING_ORIENTATION_CODLING *so);


//Evid_Threshold(1) = Coefficients(6,1)   ! Evidence threshold for Behavior 1 that increases WoE
//Evid_Threshold(2)   = Coefficients(6,2) ! Evidence threshold for Behavior 2 that increases WoE
//Evid_Threshold(3)   = Coefficients(6,3) ! Evidence threshold for Behavior 3 that increases WoE
//Evid_Threshold(4) = Coefficients(6,4)   ! Evidence threshold for Behavior 4 that increases WoE

//int IndvNumber           = sp->id;         // Individual ID
//double IdBodyLength      = sp->bodyLength; // Individual Length (meters)
//double IdBodyLngthVel(1) = sp->drift;      // Individual volitional velocity #1 (body lengths/sec), e.g., drift
//double IdBodyLngthVel(2) = sp->cruising;   // Individual volitional velocity #2 (body lengths/sec), e.g., cruising
//double IdBodyLngthVel(3) = sp->sustained;  // Individual volitional velocity #3 (body lengths/sec), e.g., sustained
//double IdBodyLngthVel(4) = sp->burst;      // Individual volitional velocity #4 (body lengths/sec), e.g., burst
//IndvEleva_NPm1 = sp->r.z;                  // elevation at individual location at t - dt
//IndvEleva_NP   = sp->sensor->r[0].z;       // elevation at individual location at current t
//IdAttrb_NP = sp->attribute_old;            // Attribute #1: 0=Individual ; 1=Invertebrate ; -1=To be removed from simulation

// arbitrary but needed value to prevent negative acceleration magnitudes (cjt :: why?)
// IndvSensoryFieldVars_NP(1,FSO,FN2)   ! Pressure at sensory point
// IndvSensoryFieldVars_NP(2,FSO,FN2)   ! Turbulent kinetic energy at sensory point
// IndvSensoryFieldVars_NP(3,FSO,FN2)   ! Acceleration magnitude at sensory point
// IndvSensoryFieldVars_NP(4,FSO,FN2)   ! THS magnitude at sensory point

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// This structure is to propagate particles with sensor points and behavoiral properties
typedef struct {
    
    SENSOR_PT sensor[NSENSORPT_MAX]; // sensor points
    int initialize_flag[NBEHAVS_MAX]; // flag to whether to initialize behavoirs
    int initialize_flag_old[NBEHAVS_MAX]; // flag to whether to initialize behavoirs
    
    int IdAttrb[NIDATTRB]; // Attribute #1: 0=Individual ; 1=Invertebrate ; -1=To be removed from simulation
    int IdAttrb_old[NIDATTRB]; // Attribute #1: 0=Individual ; 1=Invertebrate ; -1=To be removed from simulation
    
    int NumDecisions;
    
    double IdSpdRes; // movement velocity (resultant value)
    double IdSpdRes_old; // movement velocity (resultant value)
    
    // the vertical volitional movement vector angle off of (relative to) the mesh xy-plane itself; in other words, this is the
    // vertical angle off the horizontal plane. This value must range from -90.0 to 90.0.
    double SVaoMshXYZ[2];     // (S)wim (V)ector (a)ngle (o)ff Mesh xy-axis orientation
    double SVaoMshXYZ_old[2];  // (S)wim (V)ector (a)ngle (o)ff Mesh xy-axis orientation
    SVECT  SVvoMshXYZ;     // (S)wim (V)ector (v)elocity (o)ff Mesh xy-axis orientation
    SVECT  SVvoMshXYZ_old; // (S)wim (V)ector (v)elocity (o)ff Mesh xy-axis orientation
    
    double SVaoFVXY;      // (S)wim (V)ector (a)ngle (o)ff (F)low (V)ector's xy-axis
    double SVaoFVXY_old;   // (S)wim (V)ector (a)ngle (o)ff (F)low (V)ector's xy-axis
    SVECT SVvoSVXYZ;      // (S)wim (V)ector (v)elocity (o)ff (F)low (V)ector's xy-axis
    SVECT SVvoFVXYZ;      // (S)wim (V)ector (v)elocity (o)ff (F)low (V)ector's xy-axis
    SVECT SVvoFVXYZ_old;  // (S)wim (V)ector (v)elocity (o)ff (F)low (V)ector's xy-axis
    
    bool VldFVOrientXYZ[2];    // 
    bool VldSVOrientXYZ[2];    //  Valid xy-plane orientation of the individual's axis can be determined
    bool VldSVOrientXYZ_old[2]; //  Valid xy-plane orientation of the individual's axis can be determined
    
    double Info_Accumulator[NBEHAVS_MAX];       // Info_Accumulator(1:NBehavs,FN,NP)
    double Info_Accumulator_old[NBEHAVS_MAX];   // Info_Accumulator(1:NBehavs,FN,NP-1)
    double Evid_Acclimatized[NBEHAVS_MAX];      // Evid_Acclimatized(1:NBehavs,FN,NP)
    double Evid_Acclimatized_old[NBEHAVS_MAX];  // Evid_Acclimatized(1:NBehavs,FN,NP-1)
    double Evid_Acclimatized_new[NBEHAVS_MAX];  // Evid_Acclimatized(1:NBehavs,FN,NP+1)
    double Evid_Threshold[NBEHAVS_MAX];  // Evidence threshold for behaviors that increases WoE
    double Evid_Value[NBEHAVS_MAX];
    double Evid_Activity[NBEHAVS_MAX];
    
    double SVaoSVXY,SVaoSVXY_old;
    
    // Coefficients(9,1):  Decision algorithm: 1=UsherMcClelland2001 ; 2=Anderson2002
    int decision_algorithm;
    
    // Coefficients(9,2):  Sensing orientation algorithm: 1=Codling2004 ; 2=OrnsteinUhlenbeck
    int sensing_algorithm;
    
    STIMULI_PERCEPTION stimuli;
    PT_INTERACT interact;
    WALL_EFFECT wallEffect;
    BEHAVOIR_TRANS_USHER trans_usher;
    BEHAVOIR_TRANS_ANDERSON trans_anderson;
    SSPEED speed;
    SENSING_ORIENTATION_ORNSTEIN sorient_ornstein;
    SENSING_ORIENTATION_CODLING sorient_codling;
    
} SENSOR_SPECIES;

// struct methods ------------------------------

void sensor_species_alloc_init(SENSOR_SPECIES **sp);
void sensor_species_get_behavior(SENSOR_SPECIES *sp, double dt);

#endif

