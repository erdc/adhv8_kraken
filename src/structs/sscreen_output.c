#include "global_header.h"

void sscreen_output_init() {
    screen_output.residuals = OFF;
    screen_output.worse_node_nonlinear = OFF;
    screen_output.worse_node_linear = OFF;
    screen_output.grid_mass_error = OFF;
    screen_output.all = OFF;
}
