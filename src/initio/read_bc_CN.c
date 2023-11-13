#include "global_header.h"

void read_bc_CN(SMODEL *mod, char *data) {

  char line[MAXLINE];           /* the input line */
  char *subdata = NULL;         /* the data after the second card   is read */
  char *subsubdata = NULL;
  int itrn = 0;

  SIO io = *(mod->io);   // alias
  int nstring = mod->nstring;           // alias

  switch (parse_card(data, &subdata)) {
      // generic constituent
      case CARD_CON:
          itrn = get_transport_id(io, &subdata, mod->ntransport);
          mod->con[itrn].type = CON;
          mod->con[itrn].property[0] = read_dbl_field(io, &subdata);
          break;

      // salt
      case CARD_SAL:
          itrn = get_transport_id(io, &subdata, mod->ntransport);
          mod->con[itrn].type = SAL;
          mod->salinity_id = itrn;
          mod->flag.BAROCLINIC += 1;
          mod->con[itrn].property[0] = read_dbl_field(io, &subdata);
		  
          break;
		  
	  // Temperature and Heat
      case CARD_TMP:
          itrn = get_transport_id(io, &subdata, mod->ntransport);
          mod->con[itrn].type = TMP;
          mod->temperature_id = itrn;
          mod->flag.BAROCLINIC += 10;
          mod->con[itrn].property[0] = read_dbl_field(io, &subdata);
		  printf("\n Temperature is active as transport number %d\n",itrn +1);
          break;

	 // Vorticity
	  case CARD_VOR:
		  itrn = get_transport_id(io, &subdata, mod->ntransport);
		  mod->flag.VORTICITY = ON;
          mod->con[itrn].type = VOR;
		  mod->vorticity_id = itrn;
		  mod->con[itrn].property[0] = read_dbl_field(io, &subdata);
		  mod->con[itrn].property[1] = read_dbl_field(io, &subdata); // As term
		  mod->con[itrn].property[2] = read_dbl_field(io, &subdata); // Ds term
		  if (mod->con[itrn].property[1] <= 0) {
            printf("Vorticity As term is specified as %lf. This is less than 0, defaulting to 5\n",mod->con[itrn].property[1]);
			mod->con[itrn].property[1] = 5.0;
		  }
		  if (mod->con[itrn].property[2] <= 0) {
            printf("Vorticity Ds term is specified as %lf. This is less than 0, defaulting to 0.5\n",mod->con[itrn].property[2]);
			mod->con[itrn].property[2] = 0.5;
		  }
		  break;





      default:
          break;
  }

}
