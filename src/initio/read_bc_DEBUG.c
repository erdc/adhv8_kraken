
#include "global_header.h"

// these case will screen print

void read_bc_DEBUG(SMODEL *mod, char *data) {

  char line[MAXLINE];           /* the input line */
  char *subdata = NULL;         /* the data after the second card is read */
  char *subsubdata = NULL;

  SIO info = *(mod->io);

  switch (parse_card(data, &subdata)) {

  case CARD_DEBUG:
    debug.print = 1;
    break;
 
  case CARD_READBC:
    printf("Debuggin' bc read: ON\n");
    debug.readBC = 1;
    break;

  case CARD_RHS:
    debug.rhs = 1;
    break;

  case CARD_RHSHV:
    debug.rhs_hvel = 1;
    break;

  case CARD_RHSWV:
    debug.rhs_wvel = 1;
    break;

  case CARD_MAT:
    debug.matrix = 1;
    break;

  case CARD_MATHV:
    debug.matrix_hvel = 1;
    break;

  case CARD_MATWV:
    debug.matrix_hvel = 1;
    break;

  case CARD_NEWTON:
    debug.newton = 1;
    break;

  case CARD_RES:
    debug.residual = 1;
    break;

  case CARD_RESHV:
    debug.residual_hvel = 1;
    break;

  case CARD_RESWV:
    debug.residual_wvel = 1;
    break;

  case CARD_DIFF:
    debug.diffusion = 1;
    break;

  case CARD_CONV:
    debug.advection = 1;
    break;

  case CARD_SUPG:
    debug.supg = 1;
    break;

  case CARD_FRIC:
    debug.friction = 1;
    break;

  case CARD_TIME:
    debug.temporal = 1;
    break;

  case CARD_PRESS:
    debug.pressure = 1;

  case CARD_LOAD:
    debug.load = 1;
    break;

  case CARD_LOADHV:
    debug.load_hvel = 1;
    break;

  case CARD_LOADWV:
    debug.load_wvel = 1;
    break;

  case CARD_WAVE:
    debug.waves = 1;
    break;

  case CARD_WIND:
    debug.winds = 1;
    break;

  default:
    break;

  }
}
