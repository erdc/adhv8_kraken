#ifndef H_SFILE_OUTPUT_
#define H_SFILE_OUTPUT_

typedef struct {

    int bed_velocity;
    int surface_velocity;
    int depth_avg_velocity;
    int pressure;
    int grid_speed;
    int wind;
    int wave;
    int hyd_vis; /* GSAVANT */
    int trn_dif; /* GSAVANT */
    int chop;
    int grid2dm;

    // adaption output
    int adaption;
    int adapt_grid;
    int adapt_sw;
    int adapt_ns;
    int adapt_con;
    int adapt_sed;
#ifdef _ADH_GROUNDWATER
    int adapt_gw;
#endif
} SFILE_OUTPUT;


/*********************************************************/
/* struct methods -------------------------------------- */


#endif
