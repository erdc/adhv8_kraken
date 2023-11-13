#include "global_header.h"

#ifdef _ADH_GROUNDWATER
void read_bc_GW(SMODEL * mod, char *data) {
  char line[MAXLINE];           /* the input line */
  char *subdata = NULL;         /* the data after the second card   is read */
  char *subsubdata = NULL;
  CARD card;

  assert(mod);
  SIO *io = mod->io;  // alias

  assert(io && io->bc.fp);

  /** should only be called if GW detected? **/
  assert(mod->flag.GW_FLOW);

  rewind(io->bc.fp);
  while (fgets(line, MAXLINE, io->bc.fp) != NULL) {
    io_save_line(io, io->bc.fp, io->bc.filename, line);
    if (strip_comments(line) <= 1) {
      continue;
    }
    card = parse_card(line, &data);
    switch (card) {
      case CARD_MP:
	read_bc_GW_MP(mod,data);
	break;
      default:
	break;
    }/*end switch*/
  }
  rewind(io->bc.fp);
  io_save_line(io, NULL, "", "");

}


#endif
