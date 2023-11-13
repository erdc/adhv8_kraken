#include "global_header.h"

static double time_to_adulthood_secs = 3600*24*7*3;    // (~3 weeks)
static double oyster_max_fall_velocity = 0.008;        // meters/sec

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

void soyster_get_fall_velocity(SPARTICLE *p) {
    int i;
    double ratio = 0.;
    ratio = p->age / time_to_adulthood_secs;
    p->w_fall = - ratio * oyster_max_fall_velocity;
    assert(p->w_fall < 1e-6);
}
SVECT soyster_get_analytic_fall_position(SVECT p0, double t) {
    double Lmax = 300; // maximum larval length at time_to_adulthood_sec
    double Lmin = 50;  // minimum larval length at age=0
    double ratio = t / time_to_adulthood_secs;
    SVECT r_fall;
    r_fall.x = p0.x;
    r_fall.y = p0.y;
    r_fall.z = p0.z + (-0.5 * t * t / time_to_adulthood_secs) * oyster_max_fall_velocity;
    return r_fall;
}

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

void soyster_get_swim_velocity(SPARTICLE *p) {
    int i;
    double ratio = p->age / time_to_adulthood_secs;
    double Lmax = 300; // maximum larval length at time_to_adulthood_sec
    double Lmin = 50;  // minimum larval length at age=0
    svect_init(&p->v_swim);
    p->v_swim.z = (0.0089*((Lmax-Lmin)*ratio+Lmin)-0.0076)/1000.0;
    assert(p->v_swim.z > -1e-6);
}
SVECT soyster_get_analytic_swim_position(SVECT p0, double t) {
    double Lmax = 300; // maximum larval length at time_to_adulthood_sec
    double Lmin = 50;  // minimum larval length at age=0
    double ratio = t / time_to_adulthood_secs;
    SVECT r_swim;
    r_swim.x = p0.x;
    r_swim.y = p0.y;
    r_swim.z = (((4.45e-06)*Lmax - (4.45e-06)*Lmin)*t*t + (((8.9e-06)*Lmin - 7.6e-06)*t + p0.z)*time_to_adulthood_secs)/time_to_adulthood_secs;
    return r_swim;
}

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

void soyster_get_behavorial_velocity(SPARTICLE *p) {

    svect_init(&p->v_behavoir);
    
    if (p->behavoiral_vel_flag == -1) return;
    
    // only fall velocity included in behavior adustment
    if (p->behavoiral_vel_flag == 0) {
        soyster_get_fall_velocity(p);
        p->v_behavoir.z = p->w_fall;
        return;
    }
    
    // only swim velocity included in behavior adustment
    if (p->behavoiral_vel_flag == 1) {
        soyster_get_swim_velocity(p);
        p->v_behavoir.z = p->v_swim.z;
        return;
    }
    
    // both fall and swim velocity used with no random fluctuation
    if (p->behavoiral_vel_flag == 2) {
        double r = 0.64;
        soyster_get_fall_velocity(p);
        soyster_get_swim_velocity(p);
        p->v_behavoir.z = p->w_fall * (1 - r) + p->v_swim.z * r;
        return;
    }
    
    // both fall and swim velocity used with random fluctuation
    if (p->behavoiral_vel_flag == 3) {
        double r = find_random_in_range(0.64, 0.83);
        soyster_get_fall_velocity(p);
        soyster_get_swim_velocity(p);
        p->v_behavoir.z = p->w_fall * (1 - r) + p->v_swim.z * r;
        return;
    }
}

SVECT soyster_get_analytic_position(SVECT p0, double t) {
    double Lmax = 300; // maximum larval length at time_to_adulthood_sec
    double Lmin = 50;  // minimum larval length at age=0
    double ratio = t / time_to_adulthood_secs;
    SVECT r;
    r.x = p0.x;
    r.y = p0.y;
    r.z = p0.z + (((2.848e-06)*Lmax - (2.848e-06)*Lmin - 0.0014399999999999999)*t*t + ((5.696e-06)*Lmin - 4.864e-06)*t*time_to_adulthood_secs)/time_to_adulthood_secs;
    return r;
}

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

void soyster_set_wq_exposure(SPARTICLE *p, double oxygen, double salinity, double sunlight) {
    
    p->exp_oxygen = oxygen;
    p->exp_salinity = salinity;
    p->exp_sunlight = sunlight;
    
    if (p->isActive == false) return;
    
    if (p->exp_oxygen > 50)   p->isActive = false;
    if (p->exp_salinity > -0.99 && p->exp_salinity < 30) p->isActive = false;
    if (p->exp_sunlight > 25) p->isActive = false;
    
}

void  soyster_check_wq_exposure(SPARTICLE *p) {
    int flag = ON;
    if (p->exp_oxygen > 50)   flag = OFF;
    if (p->exp_salinity > -0.99 && p->exp_salinity < 30) flag = OFF;
    if (p->exp_sunlight > 25) flag = OFF;
    
    if (flag == ON) {
        return;
    } else {
        p->isActive = false;
        p->behavoiral_vel_flag = 0; // only fall velocity
    }
    
    
}

SVECT soyster_get_analytic_oxygen(SVECT p0, double t) {
    double r = 0.64;
    double tb = 132227.457346548;
    if (t > tb) {
        SVECT p02;
        p02.x = 0.; p02.y = 0.; p02.z = -70.0;
        //return soyster_get_analytic_fall_position(p02,t-tb);
        p02 = soyster_get_analytic_fall_position(p0,t);
        //p02.z *= (1-r);
        p02.z += (138.545195053645 - 70); // p(t=tb) - 70
        if (p02.z < -100) p02.z = -100;
        return p02;
    } else {
        return soyster_get_analytic_position(p0,t);
    }
}

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------

// set this up for MPI later
void soyster_printScreen(SPARTICLE p, double time) {
    int i, myid = 0;
#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif
    printf("OYSTER || myPE: %d || ID: %d(%d) || time: %f || age: %f || {x,y,z}: {%f,%f,%f} || nodeID: %d || elemID: %d || w_fall: %e || w_swim: %e || w_total: %e || exposures{oxy,sal,sun}: {%f,%f,%f} || active: %d || attached: %d || floatTime: %f \n",
           myid,
           p.gID+1,
           p.lID+1,
           time,
           p.age,
           p.r.x,
           p.r.y,
           p.r.z,
           p.nodeID,
           p.elemID,
           p.w_fall,
           p.v_swim.z,
           p.v_behavoir.z,
           p.exp_oxygen,
           p.exp_salinity,
           p.exp_sunlight,
           p.isActive,
           p.isAttached,
           p.float_time);
}

