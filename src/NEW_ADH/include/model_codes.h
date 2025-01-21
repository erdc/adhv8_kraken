//Integer codes for routine pointers


/* Variable codes */
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
#define N_TIME_STEPPERS 1
#define FE_NEWTON 0

//residual routine codes
#define N_RESID_ROUTINES 3
#define SW2 0
#define POISSON 1
#define HEAT 2

//init routine codes
#define N_INIT_ROUTINES 1
#define SW2 0


//global array of function pointers to resid routines
int (*fe_resid[N_RESID_ROUTINES])(SMODEL_SUPER *, double *, int, double, int, int, int, int);
//global array of function pointers to init routines
int (*fe_init[N_INIT_ROUTINES])(SMODEL_SUPER *);


