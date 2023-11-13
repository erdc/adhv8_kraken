#include "global_header.h"

/**************************************************************/
/**************************************************************/

void sdebug_init(SDEBUG *db) {
    
    // print these terms to screen for debugging
    db->print = 0;
    db->newton = 0;

    db->load = 0;
    db->residual = 0;
    db->matrix = 0;
    db->rhs = 0;

    db->load_hvel = 0;
    db->residual_hvel = 0;
    db->matrix_hvel = 0;
    db->rhs_hvel = 0;

    db->load_wvel = 0;
    db->residual_wvel = 0;
    db->matrix_wvel = 0;
    db->rhs_wvel = 0;

    db->rhs_transport2d = 0;
    db->rhs_transport3d = 0;

    db->diffusion = 0;
    db->advection = 0;
    db->friction = 0;
    db->supg = 0;
    db->temporal = 0;
    db->friction = 0;
    db->pressure = 0;
    db->coriolis = 0;
    db->winds = 0;
    db->waves = 0;
    db->meteor = 0;
    db->density = 0;

    db->readBC = 0;

    // leave terms out flag
    db->no_diffusion = OFF;
    db->no_advection = OFF;
    db->no_friction = OFF;
    db->no_supg = OFF;
    db->no_temporal = OFF;
    db->no_pressure = OFF;
    db->no_coriolis = OFF;
    db->no_winds = OFF;
    db->no_waves = OFF;
	db->no_hydro = OFF;
}

/**************************************************************/
/**************************************************************/

void sstop_init(SSTOP *stop) {
    stop->after_init_residual = NO;
    stop->after_init_residual_hvel = NO;
    stop->after_init_residual_wvel = NO;
}


