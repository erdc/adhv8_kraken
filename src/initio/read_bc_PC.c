#include "global_header.h"

void read_bc_PC(SMODEL *mod, char *data) {

  char line[MAXLINE];           /* the input line */
  char *subdata = NULL;         /* the data after the second card   is read */
  char *subsubdata = NULL;

  switch (parse_card(data, &subdata)) {

  case CARD_LVL:
    mod->out_level = read_int_field(*(mod->io), &subdata);
    break;

  case CARD_XDF:
#ifdef _ADH_HDF5
    mod->flag.PC_FILE_XDMF = ON;
#ifdef _MESSG
    if (mod->grid->smpi->myid==0)
#endif
        printf("\nSwitching output type to XDMF ...");

#else
    printf("\nWarning: Attempted to turn on XDMF output, but code was not compiled with HDF5 support.");
    printf("\nIgnoring XDMF card; Default AdH output will be used.");
    //tl_error("Attempted to turn on XDMF output, but code was not compiled with HDF5 support.\n");
#endif
    break;

  default:
    break;
  }
}
