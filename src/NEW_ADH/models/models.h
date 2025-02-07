#ifndef _H_MODELS_
#define _H_MODELS_

//Integer codes for routine pointers
/* Variable codes */
//Think corey has handle on this and so will go away
#define PERTURB_U 1
#define PERTURB_V 2
#define PERTURB_W 3
#define PERTURB_DPL 4
#define PERTURB_H 5
#define PERTURB_C 6
#define PERTURB_P 7
#define PERTURB_D 8
#define PERTURB_NONE -1


//forward step pointer codes
//see time_loop
//this should go in another file
#define N_TIME_STEPPERS 1
#define FE_NEWTON 0

//total number of resid routines
#define N_RESID_ROUTINES 18
//1d residual routine codes
#define N_1D_RESID_ROUTINES 15
//2d residual routine codes
#define N_2D_RESID_ROUTINES 3
//3d residual routines
#define N_3D_RESID_ROUTINES 0
//indices for residual routines
#define SW2_BC_DISCHARGE 0
#define SW2_BC_ELE 1
#define SW2_BC_FLAPD 2
#define SW2_BC_FLAPU 3
#define SW2_BC_FLUX 4
#define SW2_BC_H 5
#define SW2_BC_HYBRID 6
#define SW2_BC_OUTFLOW 7
#define SW2_BC_SLUICED 8
#define SW2_BC_SLUICEU 9
#define SW2_BC_TAILWATER 10
#define SW2_BC_VEL 11
#define SW2_BC_VEL_ELE 12
#define SW2_BC_WEIRD 13
#define SW2_BC_WEIRU 14
#define SW2 15
#define POISSON 16
#define HEAT 17

//init routine codes (for prepping dependent values before time step)
#define N_INIT_ROUTINES 1
#define SW2_INIT 0

//global array of function pointers to resid routines
int (*fe_resid[N_2D_RESID_ROUTINES])(SMODEL_SUPER *, double *, int, double, int, int, int, int);


//global array of function pointers to init routines
//since not in generic assembly loop, dimension really doesnt matter
int (*fe_init[N_INIT_ROUTINES])(SMODEL_SUPER *);


void set_models(int (*fe_resid[N_RESID_ROUTINES])(SMODEL_SUPER *, double *, int, double, int, int, int, int), int (*fe_init[N_INIT_ROUTINES])(SMODEL_SUPER *));
#endif
