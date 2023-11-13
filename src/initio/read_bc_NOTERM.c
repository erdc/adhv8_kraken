
#include "global_header.h"

// these case will prevent certain terms from being included in the rhs 
// in 3d, this is for the hvel (momentum equations)


void read_bc_NOTERM(SMODEL *mod, char *data) {

  char line[MAXLINE];           /* the input line */
  char *subdata = NULL;         /* the data after the second card is read */
  char *subsubdata = NULL;

  SIO info = *(mod->io);

  switch (parse_card(data, &subdata)) {

  case CARD_HYDRO:
	debug.no_hydro = ON;
	debug.u_vel = read_dbl_field(info, &subdata);
	debug.v_vel = read_dbl_field(info, &subdata);
	debug.w_vel = read_dbl_field(info, &subdata);
    debug.hard_disp = read_dbl_field(info, &subdata);
    mod->max_nsys = 1;
    mod->max_nsys_sq = 1;
	break;

  case CARD_DIFF:
    debug.no_diffusion = ON;
    break;

  case CARD_CONV:
    debug.no_advection = ON;
    break;

  case CARD_SUPG:
    debug.no_supg = ON;
    break;

  case CARD_FRIC:
    debug.no_friction = ON;
    break;

  case CARD_TIME:
    debug.no_temporal = ON;
    break;

  case CARD_PRESS:
    debug.no_pressure = ON;
    break;
  
  case CARD_COR:
    debug.no_coriolis = ON;
    break;

  default:
    break;

  }
}
