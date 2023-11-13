#ifndef H_SOYSTER_
#define H_SOYSTER_

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// struct methods ++++++++++++++++++++++++++++++++++++++++++
void  soyster_check_wq_exposure(SPARTICLE *p);
void soyster_get_fall_velocity(SPARTICLE *p);
SVECT soyster_get_analytic_fall_position(SVECT p0, double t);
void soyster_get_swim_velocity(SPARTICLE *p);
SVECT soyster_get_analytic_swim_position(SVECT p0, double t);
void soyster_get_behavorial_velocity(SPARTICLE *p);
SVECT soyster_get_analytic_position(SVECT p0, double t);
void soyster_printScreen(SPARTICLE p, double time);
void soyster_set_wq_exposure(SPARTICLE *p, double oxygen, double salinity, double sunlight);
SVECT soyster_get_analytic_oxygen(SVECT p0, double t);
#endif
