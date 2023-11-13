#ifndef H_SPARTICLE_
#define H_SPARTICLE_

typedef struct {

    int gID;             // global ID
    int lID;             // local ID
    int ptype;           // type of particle; 0 - general particle, 1 - oyster
    int nodeID;          // the nodeID at time t if the element resides on a node
    int elemID;          // the elemID of the particle at time t
    SVECT r;             // particle position vector
    double w_fall;       // particle fall velocity
    SVECT  v_swim;       // particle swim velocity
    SVECT  v_behavoir;   // particle total behavorial velocity adjustment
    double exp_oxygen;   // particle dissolved oxygen exposure
    double exp_salinity; // particle salinity exposure
    double exp_sunlight; // particle sunlight exposure
    bool isActive;        // is the particle alive/dead?
    double age;          // time the particle has been active in the simulation
    bool isAttached;     // is particle attached to a grid feature
    double float_time;   // how long as the particle been floating for
    int behavoiral_vel_flag; // flag for species velocity adjustment 
    
} SPARTICLE;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// struct methods ++++++++++++++++++++++++++++++++++++++++++

void sparticle_init(SPARTICLE *p, int np);
void sparticle_printScreen(SPARTICLE *p, int np, double time);
void sparticle_set_ptype(SPARTICLE *p, int ptype);
void sparticle_set_nodeID(SPARTICLE *p, int nodeID);
void sparticle_set_elemID(SPARTICLE *p, int elemID);
void sparticle_set_position(SPARTICLE *p, SVECT pos);
void sparticle_set_wfall(SPARTICLE *p, double wfall);
void sparticle_set_v_swim(SPARTICLE *p, SVECT vswim);
void sparticle_set_exp_oxygen(SPARTICLE *p, double exp_oxygen);
void sparticle_set_exp_salinity(SPARTICLE *p, double exp_salinity);
void sparticle_set_exp_sunlight(SPARTICLE *p, double exp_sunlight);
void sparticle_set_activity(SPARTICLE *p, int status);
void sparticle_set_age(SPARTICLE *p, double age);

#endif
