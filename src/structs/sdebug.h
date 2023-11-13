#ifndef H_SDEBUG2_
#define H_SDEBUG2_

typedef struct {

    // print these terms to screen for debugging
    int print;
    int newton;
    
    // sw2d
    int load;
    int residual;
    int matrix;
    int rhs;

    // sw3d hvel
    int load_hvel;
    int residual_hvel;
    int matrix_hvel;
    int rhs_hvel;

    // sw3d wvel
    int load_wvel;
    int residual_wvel;
    int matrix_wvel;
    int rhs_wvel;

    // transport
    int rhs_transport2d;
    int rhs_transport3d;

    // physics terms
    int diffusion;
    int advection;
    int friction;
    int supg;
    int temporal;
    int pressure;
    int coriolis;
    int winds;
    int waves;
    int meteor;
    int density;
	int vorticity;

    // boundary condition read
    int readBC;

    // leave these terms out
    int no_diffusion;
    int no_advection;
    int no_friction;
    int no_supg;
    int no_temporal;
    int no_pressure;
    int no_coriolis;
    int no_winds;
    int no_waves;
	int no_hydro;
	double u_vel, v_vel, w_vel, hard_disp;

} SDEBUG;


// these flags are the stop AdH at certain points during calculation
typedef struct {
    int after_init_residual;
    int after_init_residual_hvel;
    int after_init_residual_wvel;

} SSTOP;


/*********************************************************/
/* struct methods -------------------------------------- */

void sdebug_init(SDEBUG *);

#endif
