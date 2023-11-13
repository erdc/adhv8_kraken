#ifndef H_SMODEL_
#define H_SMODEL_

typedef struct {
    
    double dt;   // timestep (s)
    double t0;   // start of simulation (s)
    double tf;   // end of simulation (s)
    double time; // current time of the simulation (s)
    double max_t_abs_error_mag; // max grid absolute error over time
    double max_t_rel_error_mag; // max grid relative error over time
    
    int np_global; // the total number of particles across processors
    int np;        // the total number of particles on PE
    int time_update_flag; // a flag for user chosen time discretization method
    int spool;     // the number of time-steps between writes
    int npoly;     // the number of polygons/reefs used in the model
    bool screen_print_timestep; // write data to screen at each timestep

    bool normalize; // normalize coordinates to [0, 1] \times [0, 1]?
    SGRID *grid;  // the hydrodynamic grid
    SPARTICLE *p; // point to an array of particles
    SPOLYGON *poly; // an array of domain polygons to hold reefs, etc.
    
    // input files
    char file_base[100];
    char vel_filename[100];
    char wq_filename[100];
    char dpl_filename[100];
    FILE *fp_dpl;      // a file pointer to hydro node grid displacements in time
    FILE *fp_vel;      // a file pointer to the hydrodynamic velocity
    FILE *fp_wq;       // a file pointer to water quality field variables
    FILE *fp_salinity; // a file pointer to water salinity
    FILE *fp_doxygen;  // a file pointer to water dissolved oxygen
    FILE *fp_exposure; // a file pointer to sunlight exposure
    
    int time_units_dep; // 0(s), 1(m), 2(h), 3(d)
    int time_units_vel; // 0(s), 1(m), 2(h), 3(d)
    int time_units_wq; // 0(s), 1(m), 2(h), 3(d)
    int time_units_salinity; // 0(s), 1(m), 2(h), 3(d)
    int time_units_deoxygen; // 0(s), 1(m), 2(h), 3(d)
    int time_units_expossure; // 0(s), 1(m), 2(h), 3(d)
    
    // output files
    FILE *fp_out;      // a csv formatted ascii output file of particle data
    FILE *fp_error;    // a file pointer for error printing
    
    SVECT *p0;    // initial particle locations for displacement plotting
    SVECT *pdpl;   // particle displacement from origin
    
    // displacement snap data
    double tL_dpl, tR_dpl;
    double *dpl_tL;   // grid velocity on left snap side
    double *dpl_tR;   // grid velocity on right snap side
    double *dpl; // verticle node displacement at time t
    double *dpl1, *dpl2, *dpl3, *dpl4, *dpl5;  // vertical node displacement at RK times
    
    // velocity snap data
    double tL, tR;
    SVECT *vel_tL;   // grid velocity on left snap side
    SVECT *vel_tR;   // grid velocity on right snap side
    SVECT *vel;      // grid velocity interpolated @ current time
    
    // water quality snap data
    double tL_wq, tR_wq;
    double *oxygen,    *salinity,    *sunlight;
    double *oxygen_tL, *salinity_tL, *sunlight_tL;
    double *oxygen_tR, *salinity_tR, *sunlight_tR;
    
    // analytic functions
    void *get_analytic_p;  //
    void *get_analytic_v;  //
    void *get_analytic_wq; //
    
#ifdef _ADH_HDF5
    HDF5 hdf5; // for binary XDMF output
#endif
    
    
} SMODEL;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// struct methods ++++++++++++++++++++++++++++++++++++++++++
int smodel_init_noRead(SMODEL *mod, SGRID *grid, char *file_base, double t0, double tf, double dt, int np, int *node_ID, int *elem_ID, int spool, int ptype);
int smodel_init(SMODEL *mod);
int smodel_finalize(SMODEL *mod);
int smodel_read_input(SMODEL *mod);
int smodel_run_lockstep(SMODEL *mod);
void smodel_print_ts(SMODEL *mod, int it);

#endif
