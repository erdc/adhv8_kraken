#include "global_header.h"
void sfile_output_init(SFILE_OUTPUT *file_output) {

    file_output->bed_velocity = OFF;
    file_output->surface_velocity = OFF;
    file_output->depth_avg_velocity = OFF;
    file_output->pressure = OFF;
    file_output->grid_speed = OFF;
    file_output->wind = OFF;
    file_output->wave = OFF;
    file_output->grid2dm = OFF;

    file_output->adaption = OFF;
    file_output->adapt_grid = OFF;
    file_output->adapt_sw = OFF;
    file_output->adapt_ns = OFF;
    file_output->adapt_con = OFF;
    file_output->adapt_sed = OFF;
#ifdef _ADH_GROUNDWATER
    file_output->adapt_gw  = OFF;
#endif
}

