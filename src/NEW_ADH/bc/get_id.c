#include "adh.h"

int get_id(char **token, int nmax, char *error_descript) {
    
    /* gets and checks the material number */
    int i = get_next_token_int(token);
    if (i > nmax || i < 1) {io_read_error(io, error_descript, TRUE);}
    i--;
    return i;
}
