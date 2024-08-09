#ifndef H_SSERIES_
#define H_SSERIES_

/***********************************************************/
/***********************************************************/
/* stores wind time-series data */
typedef struct {
  double time;       /* time */
  double *area;      /* an array of nvalue areas */
  double *slope;     /* an array of nvalue slope */
  double *value;     /* an array of nvalue values */
} SSERIES_ENTRY;

/***********************************************************/
/***********************************************************/

/* an xy time-series node for a linked list */
typedef struct sts {
  int id;                   /* series id */
  int size;                 /* the number of items in the series */
  int type;                 /* the type of series */
  int nvalues;              /* the number of values on entry line */
  int nnodes;               /* the number of 2d grid nodes for node contribution */
  int mat_id;               /* material id associated with series (for water sourcing in 2d) */
  double infact;            /* conversion factor to seconds */
  double outfact;           /* conversion factor from seconds */
  double tol;               /* the tolerance of the series */
  double *ivalue;           /* current interpolation values */
  SMETEOR_STATION *station; /* meteorologic station tied to the series if there is one */
  SSERIES_ENTRY *entry;     /* an array of series entries */
  struct sts *next;         /* next item in linked list */
  struct sts *prev;         /* previous item in linked list */

} SSERIES;

/*********************************************************/
/* struct methods -------------------------------------- */


void sseries_entry_init(SSERIES_ENTRY *);
void sseries_init(SSERIES *);
void sseries_alloc(SSERIES **, int, int, int);
void sseries_free(SSERIES *);
void sseries_cpy(SSERIES *, SSERIES *);
void sseries_copy_entry(SSERIES_ENTRY *, SSERIES_ENTRY, int);
void sseries_check(SSERIES);
void sseries_checklist(int, SSERIES *);
void sseries_printScreen(SSERIES, int);
void sseries_printScreen_list(int, SSERIES *);
SSERIES *sseries_create_list(SSERIES *, SSERIES **, SSERIES **);
SSERIES *sseries_add(SSERIES *, SSERIES **, SSERIES **, bool);
SSERIES *sseries_search(int, SSERIES **, SSERIES *);
int sseries_delete(int, SSERIES **, SSERIES **);
void sseries_free_list(SSERIES **);
double sseries_get_value(int, SSERIES *, int);
void sseries_set_values(SSERIES *, double);
void sseries_setall_ts_values(SSERIES *, double);
void sseries_set_values_average(SSERIES *, double, double);
void sseries_setall_ts_valuesAVG(SSERIES *, double, double);
SSERIES *sseries_extract(SSERIES *, int);
int sseries_read_id(SIO io, char **data, int);
//int sseries_set_type(SMODEL_SUPER *, char **, int);
int sseries_get_interval(SSERIES, double);
int series_get_type_count(SSERIES *, int);
void sseries_set_meteor_stations(SSERIES *, SGRID *, int);
void series_list_update_values(SSERIES *series_head, double t_prev);
double series_list_get_value_from_mat(SSERIES *, int);

/***********************************************************/
/***********************************************************/
/***********************************************************/

#endif

