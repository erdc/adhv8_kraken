#ifndef H_SSCREEN_OUTPUT_
#define H_SSCREEN_OUTPUT_

typedef struct {

    int residuals;
    int worse_node_nonlinear;
    int worse_node_linear;
    int grid_mass_error;
    int all;

} SSCREEN_OUTPUT;


/*********************************************************/
/* struct methods -------------------------------------- */

void sscreen_output_init();

#endif

