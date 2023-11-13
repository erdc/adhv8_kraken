//////////////////////////////////////////////////////////////////////////////////////////////////
// FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. //
#ifndef H_SWINDLIB_
#define H_SWINDLIB_


/* Note: Entire file added by Gajanan. */
typedef struct {
  int nws;                  /* nws type */
  int ics;                  /* ics = 1 (Cartesian) for AdH necessarily as of July 2015. */
  int NScreen;
  int MyProc;
  int ScreenUnit;

  /* Timing information */
  double statim;            /* wind start time */
  double wtiminc;           /* wind time increment */
  double wtime_prev;        /* previous wind time */
  double wtime_curr;        /* current wind time */

  /* nws dependent information */
    /* grid info */
    int nwlon;                /* number of longitudes (x) */
    int nwlat;                /* number of latitudes (y) */
    double wlonmin;           /* minimum longitude value (x_min) */
    double wlatmax;           /* maximum latitude value (y_max) */
    double wloninc;           /* longitude increment */
    double wlatinc;           /* latitude increment */

    /* time conversion */
    int iwyr;
    int irefyr;
    int irefmo;
    int irefday;
    int irefhr;
    int irefmin;
    int iwtime;
    double refsec;
    double wreftim;
    double wtimed;

  /* ics dependent information; currently not in use since ICS = 1 */
//  double slam0;               /* ADCIRC fort.15 input parameter */
//  double sfea0;               /* ADCIRC fort.15 input parameter */

  /* Actual wind variables in x and y direction. Windx, Windy, Pressure */
  double *prevx;              /* 1d array of size model->nnodes */
  double *prevy;              /* 1d array of size model->nnodes */
  double *prevp;              /* 1d array of size model->nnodes */
  double *currx;              /* 1d array of size model->nnodes */
  double *curry;              /* 1d array of size model->nnodes */
  double *currp;              /* 1d array of size model->nnodes */

} SWINDLIB;

/*********************************************************/
/* struct methods -------------------------------------- */

void swindlib_alloc(SWINDLIB **windlib);
void swindlib_init(SWINDLIB *windlib, int nnodes);
void swindlib_printScreen(SWINDLIB *windlib);
void swindlib_free(SWINDLIB *windlib, int nnodes);

// ************************************************************//
// ************************************************************//
// local prototypes *******************************************//

void AdH_direct_nwsget(double *wvnx, double *wvny, double *press, int nnodes, SIO *io);
void AdH_ERA_nwsget(double *wvnx, double *wvny, double *press, SGRID *grid, SIO * io);

#ifdef WIND_GNU_NM
#define    windlib_nws3get(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s) \
        __wind_MOD_nws3get(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s)

#define    windlib_nws4get(a, b, c, d, e, f) \
        __wind_MOD_nws4get(a, b, c, d, e, f)

#define    windlib_nws6get(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q) \
        __wind_MOD_nws6get(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q)

#define    windlib_nws7get(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q) \
        __wind_MOD_nws7get(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q)

#define    windlib_timeconv(a, b, c, d, e, f, g, h, i, j) \
        __wind_MOD_timeconv(a, b, c, d, e, f, g, h, i, j)

void __wind_MOD_nws3get(double *x, double *y, double *slam, double *sfea, double *wvnx, double *wvny,
                        int *iwtime, int *iwvyr, double *wtimed, int *nnodes,
                        int *nwlon, int *nwlat, double *wlatmax, double *wlonmin,
                        double * wlatinc, double *wloninc, int *ics, int *nscreen, int *screenunit);

void __wind_MOD_nws4get(double *wvnx, double *wvny, double *press, int *nnodes, double *rho, double *g);

void __wind_MOD_nws6get(double *x, double *y, double *slam, double *sfea, double *wvnx, double *wvny,
                        double *press, int *nnodes, int *nwlon, int *nwlat, double *wlatmax,
                        double *wlonmin, double * wlatinc, double *wloninc, int *ics,
                        double *rho, double *g);

void __wind_MOD_nws7get(double *x, double *y, double *slam, double *sfea, double *wvnx, double *wvny,
                        double *press, int *nnodes, int *nwlon, int *nwlat, double *wlatmax,
                        double *wlonmin, double * wlatinc, double *wloninc, int *ics,
                        double *rho, double *g);

void __wind_MOD_timeconv(int *iyr, int *imo, int *iday, int *ihr, int *imin, double *isec,
                         double *timesec, int *MyProc, int *NScreen, int *ScreenUnit);
#endif
#ifdef WIND_INTEL_NM
#define  windlib_nws3get(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s) \
        wind_mp_nws3get_(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s)

#define  windlib_nws4get(a, b, c, d, e, f) \
        wind_mp_nws4get_(a, b, c, d, e, f)

#define  windlib_nws6get(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q) \
        wind_mp_nws6get_(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q)

#define  windlib_nws7get(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q) \
        wind_mp_nws7get_(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q)

#define  windlib_timeconv(a, b, c, d, e, f, g, h, i, j) \
        wind_mp_timeconv_(a, b, c, d, e, f, g, h, i, j)

void wind_mp_nws3get_(double *x, double *y, double *slam, double *sfea, double *wvnx, double *wvny,
                      int *iwtime, int *iwvyr, double *wtimed, int *nnodes,
                      int *nwlon, int *nwlat, double *wlatmax, double *wlonmin,
                      double * wlatinc, double *wloninc, int *ics, int *nscreen, int *screenunit);

void wind_mp_nws4get_(double *wvnx, double *wvny, double *press, int *nnodes, double *rho, double *g);

void wind_mp_nws6get_(double *x, double *y, double *slam, double *sfea, double *wvnx, double *wvny,
                      double *press, int *nnodes, int *nwlon, int *nwlat, double *wlatmax,
                      double *wlonmin, double * wlatinc, double *wloninc, int *ics,
                      double *rho, double *g);

void wind_mp_nws7get_(double *x, double *y, double *slam, double *sfea, double *wvnx, double *wvny,
                      double *press, int *nnodes, int *nwlon, int *nwlat, double *wlatmax,
                      double *wlonmin, double * wlatinc, double *wloninc, int *ics,
                      double *rho, double *g);

void wind_mp_timeconv_(int *iyr, int *imo, int *iday, int *ihr, int *imin, double *isec,
                       double *timesec, int *MyProc, int *NScreen, int *ScreenUnit);
#endif
#ifdef WIND_PGI_NM
#define  windlib_nws3get(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s) \
        wind_nws3get_(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s)

#define  windlib_nws4get(a, b, c, d, e, f) \
        wind_nws4get_(a, b, c, d, e, f)

#define  windlib_nws6get(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q) \
        wind_nws6get_(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q)

#define  windlib_nws7get(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q) \
        wind_nws7get_(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q)

#define  windlib_timeconv(a, b, c, d, e, f, g, h, i, j) \
           wind_timeconv_(a, b, c, d, e, f, g, h, i, j)

void wind_nws3get_(double *x, double *y, double *slam, double *sfea, double *wvnx, double *wvny,
                      int *iwtime, int *iwvyr, double *wtimed, int *nnodes,
                      int *nwlon, int *nwlat, double *wlatmax, double *wlonmin,
                      double * wlatinc, double *wloninc, int *ics, int *nscreen, int *screenunit);

void wind_nws4get_(double *wvnx, double *wvny, double *press, int *nnodes, double *rho, double *g);

void wind_nws6get_(double *x, double *y, double *slam, double *sfea, double *wvnx, double *wvny,
                      double *press, int *nnodes, int *nwlon, int *nwlat, double *wlatmax,
                      double *wlonmin, double * wlatinc, double *wloninc, int *ics,
                      double *rho, double *g);

void wind_nws7get_(double *x, double *y, double *slam, double *sfea, double *wvnx, double *wvny,
                      double *press, int *nnodes, int *nwlon, int *nwlat, double *wlatmax,
                      double *wlonmin, double * wlatinc, double *wloninc, int *ics,
                      double *rho, double *g);
void wind_timeconv_(int *iyr, int *imo, int *iday, int *ihr, int *imin, double *isec,
                    double *timesec, int *MyProc, int *NScreen, int *ScreenUnit);
#endif

/***********************************************************/
/***********************************************************/
/***********************************************************/

// ABOVE LINES ADDED BY GAJANAN                                                                 //
//////////////////////////////////////////////////////////////////////////////////////////////////
#endif
