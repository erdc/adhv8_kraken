#ifndef H_SMETEOR_
#define H_SMETEOR_

/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {
    int n;              /* total number of grid points*/
    double dt;          /* time between snaps */
    double tprev;       /* previous time */
    double tnext;       /* next time */
    FILE *fin;          /* file pointer for meteorologicreads */
} SMETEOR_FILE;

/***********************************************************/
/***********************************************************/
/***********************************************************/

typedef struct {
    double x,y;                // station location
    double *node_contrib;      // grid interpolation weights
} SMETEOR_STATION;

/***********************************************************/
/***********************************************************/
/***********************************************************/
// use winds.stress.x
typedef struct {
    SVECT2D stress;
} SWIND;

/***********************************************************/
/***********************************************************/
/***********************************************************/
// use: wave.rads.xx
typedef struct {
    STENSOR2D rads;
    SVECT2D stress;
} SWAVE;

/*********************************************************/
/* struct methods -------------------------------------- */

void smeteor_file_init(SMETEOR_FILE *);
void smeteor_station_init(SMETEOR_STATION *, int);
void smeteor_station_alloc(SMETEOR_STATION **, int);
void smeteor_station_realloc(SMETEOR_STATION **, int, int);
void smeteor_station_free(SMETEOR_STATION *, int);
void smeteor_station_copy(SMETEOR_STATION *, SMETEOR_STATION *, int);
void swind_init(SWIND *);
void swind_alloc(SWIND **, int);
void swind_free(SWIND *, int);
void swave_init(SWAVE *);
void swave_alloc(SWAVE **, int);
void swave_free(SWAVE *, int);
void swave_elem2d_local_init(SWAVE *);
void swind_elem2d_local_init(SWIND *);
void swind_elem2d_get_local(SWIND *, SWIND *, int *);
void swave_elem2d_get_local(SWAVE *, SWAVE *, int *);
void swave_realloc(SWAVE **, int, int);
void swind_realloc(SWIND **, int, int);
//void swave_renumber(SWAVE *, int, int *, int *);
//void swind_renumber(SWIND *, int, int *, int *, SVECT2D *);

//SVECT2D swave_elem2d_get_local_stress(int ndim, int CSTORM_WSID, SWAVE *waves, SELEM_2D elem2d, int *grid2d_nodeIDs);
//SVECT2D swind_elem2d_get_local_stress(int ndim, SWIND *winds, SELEM_2D *elem2d, int *grid2d_nodeIDs, double elem_depth_avg, double gravity, double density, double     windatt, int wind_flag, int WIND_LIBRARY, int NWS);


#endif

