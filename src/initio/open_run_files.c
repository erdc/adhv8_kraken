#include "global_header.h"

#define FND  "FOUND"
#define NFND "NOT FOUND"
#define SPC  " (specified)"
#define DEF  " (default)"

/******************************************************************/
/******************************************************************/
/******************************************************************/

void open_input_file(SFILE_IN *file, char *type, int super)
{
  char msg[MAXLINE];
  int error = FALSE;

  (file->fp) = io_fopen(file->filename, "r", FALSE);
  if (! (file->fp) ) {
    error = TRUE;
  }
  sprintf(msg, " %s%s: %s - %s", type, (super ? (file->specified ? SPC : DEF) : ""), file->filename, (file->fp ? FND : NFND));
  root_print(msg);

  /* error handling */
  if (error) {
    fprintf(stderr, "ERROR: Could not find the file: %s\n",file->filename);
    exit(0);
  }
}

/******************************************************************/
/******************************************************************/
/******************************************************************/

void open_output_file(SFILE *file, char *type, int super) {

  int error = FALSE;

  (file->fp) = io_fopen(file->filename, "w", FALSE);
  if (! (file->fp) ) {
    error = TRUE;
  }
}

/******************************************************************/
/******************************************************************/
/******************************************************************/

