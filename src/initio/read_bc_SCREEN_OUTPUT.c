#include "global_header.h"

// these case will stop the code

void read_bc_SCREEN_OUTPUT(SMODEL * mod, char *data)
{

    char line[MAXLINE];         /* the input line */
    char *subdata = NULL;       /* the data after the second card is read */
    char *subsubdata = NULL;

    SIO info = *(mod->io);

    switch (parse_card(data, &subdata)) {
        case CARD_RESID:
            screen_output.residuals = 1;
            break;

        case CARD_NLNODE:
            screen_output.worse_node_nonlinear = 1;
            break;

        case CARD_LNODE:
            screen_output.worse_node_linear = 1;
            break;

        case CARD_MERROR:
            screen_output.grid_mass_error = 1;
            break;

        case CARD_ALL:
            screen_output.residuals = 1;
            screen_output.worse_node_nonlinear = 1;
            screen_output.worse_node_linear = 1;
            //screen_output.grid_mass_error = 1;

        default:
            break;
    }
}
