/* closes a data set */

#include "global_header.h"

void print_trailer(FILE * fp_out    /* the output file */
    )
{
    /* writes the end card */
    fprintf(fp_out, "ENDDS\n");
    /* closes the file */
    fclose(fp_out);
    fp_out = NULL;
}
