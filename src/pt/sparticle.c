#include "global_header.h"


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void sparticle_init(SPARTICLE *p, int np) {
    int i;
    for (i=0; i<np; i++) {
        p[i].gID = -1;
        p[i].lID = -1;
        p[i].ptype = -1;
        p[i].nodeID = -1;
        p[i].elemID = -1;
        p[i].behavoiral_vel_flag = -1;
        p[i].r.x = -999999.;
        p[i].r.y = -999999.;
        p[i].r.z = -999999.;
        p[i].w_fall = 0.;
        p[i].v_swim.x = 0.;
        p[i].v_swim.y = 0.;
        p[i].v_swim.z = 0.;
        p[i].v_behavoir.x = 0.;
        p[i].v_behavoir.y = 0.;
        p[i].v_behavoir.z = 0.;
        p[i].exp_oxygen   = -1;
        p[i].exp_salinity = -1;
        p[i].exp_sunlight = -1;
        p[i].isActive = -1;
        p[i].age = 0;
        p[i].isAttached = false;
        p[i].float_time = 0.;
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void sparticle_printScreen(SPARTICLE *p, int np, double time) {
    int i;
    for (i=0; i<np; i++) {
        if (p[i].ptype == 1) {
            soyster_printScreen(p[i],time);
        } else {
            printf("t: %10.5f \t id: %7d \t {%10.5e, %10.5e, %10.5e}\n ",time,i,p[i].r.x,p[i].r.y,p[i].r.z);
        }
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void sparticle_set_type(SPARTICLE *p, int ptype) {
    p->ptype = ptype;
}
void sparticle_set_nodeID(SPARTICLE *p, int nodeID) {
    p->ptype = nodeID;
}
void sparticle_set_elemID(SPARTICLE *p, int elemID) {
    p->ptype = elemID;
}
void sparticle_set_position(SPARTICLE *p, SVECT pos) {
    p->r.x = pos.x;
    p->r.y = pos.y;
    p->r.z = pos.z;
}
void sparticle_set_wfall(SPARTICLE *p, double wfall) {
    p->w_fall = wfall;
}
void sparticle_set_v_swim(SPARTICLE *p, SVECT vswim) {
    p->v_swim = vswim;
}
void sparticle_set_exp_oxygen(SPARTICLE *p, double exp_oxygen) {
    p->exp_oxygen = exp_oxygen;
}
void sparticle_set_exp_salinity(SPARTICLE *p, double exp_salinity) {
    p->exp_salinity = exp_salinity;
}
void sparticle_set_exp_sunlight(SPARTICLE *p, double exp_sunlight) {
    p->exp_sunlight = exp_sunlight;
}
void sparticle_set_activity(SPARTICLE *p, int status) {
    p->isActive = status;
}
void sparticle_set_age(SPARTICLE *p, double age) {
    p->age = age;
}

