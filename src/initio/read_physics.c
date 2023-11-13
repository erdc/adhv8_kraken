/* reads the sizes from the boundary condition file and allocates the memory */

#include "global_header.h"

void read_physics(SMODEL *mod) {
  char line[MAXLINE];           /* the input line */
  char *data = NULL;            /* the data after the first card is read */
  char *subdata = NULL;         /* the data after the second card is read */
  int i, j, k;                  /* loop counters */

  assert(mod->io->bc.fp);
  while (fgets(line, MAXLINE, mod->io->bc.fp) != NULL) {
      io_save_line(mod->io, mod->io->bc.fp, mod->io->bc.filename, line);
      if (strip_comments(line) <= 1) { /* ignore empty ('\n') or comment lines */
          continue;
      }
      switch (parse_card(line, &data)) {
          case CARD_OP:
            switch (parse_card(data, &subdata)) {
                case CARD_SW2:
                  mod->flag.SW2_FLOW = TRUE;
                  break;
                case CARD_SW3:
                  mod->flag.SW3_FLOW = TRUE;
                  break;
   	        case CARD_GW:
		  mod->flag.GW_FLOW = TRUE;
		  break;
                default:
                  break;
            }
            break;
          default:
            break;
      }
  }
  rewind(mod->io->bc.fp);
  io_save_line(mod->io, NULL, "", "");

    if (mod->flag.SW2_FLOW || mod->flag.SW3_FLOW)
        mod->flag.SW_FLOW = ON;
    if (mod->flag.SW2_TRANSPORT || mod->flag.SW3_TRANSPORT)
        mod->flag.TRANSPORT = ON;
    if (mod->flag.SW3_FLOW)
        mod->flag.MG = ON;
#ifdef _ADH_GROUNDWATER
    if (mod->flag.GW_TRANSPORT)
        mod->flag.TRANSPORT = ON;
#endif 

}
