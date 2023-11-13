
#ifndef H_SSTR_VALUE_
#define H_SSTR_VALUE_

// dependencies :: 

// /***********************************************************/
// /***********************************************************/
// /***********************************************************/

typedef struct
{
  /* if you add another parameter here, update initialize_input_data_type() */
  int bc_flag;          /* the type of boundary condition */
  int src_flag;         /* the type of source */
  int isigma;           /* flux specification */
  int iu_0;         /* solution specification */
  int iq;           /* the strength of the source */
  int ici;          /* concentration of the source */
  int ivx;          /* the boundary specification series number for x */
  int ivy;          /* the boundary specification series number for y */
  int ivz;          /* the boundary specification series number for z */
  double ceq_0;         /* the boundary equilibrium concentration */
} INPUT_DATA;           /* input data structure */

typedef struct
{
  /* if you add another parameter here, update initialize_friction_type() */
  double manningsn;     /* Manning's n roughness */
  double skfricoef;     /* skin friction coefficient */
  double dragcoef;      /* drag coefficient */
  double eqrheight;     /* equivalent roughness height */
  double rheight;       /* roughness height (K-E) */
  double hghtstem;      /* undeflected stem height */
  double diamstem;      /* stem diameter */
  double densstem;      /* stem density */
  double icethick;      /* thickness of the ice */
  double icedense;      /* density of the ice */
  double bedrhght;      /* bed roughness height */
  double icerhght;      /* ice roughness height */
  int mng_flag;         /* flags YES if Manning's n roughness is input */
  int erh_flag;         /* flags YES if roughness height is input */
  int sav_flag;         /* flags YES for submerged aquatic vegetation */
  int urv_flag;         /* flags YES for unsubmerged rigid vegetation */
  int icemors;          /* 0 = stationary ice, 1 = moving ice */
} FRICTION;


typedef struct {
  INPUT_DATA flow;      /* the flow data */
  INPUT_DATA heat;      /* the thermal data */
  INPUT_DATA ol_flow;       /* the overland flow data */
  INPUT_DATA ch_flow;       /* the channel flow data */
  INPUT_DATA pressure;      /* the pressure data */
  INPUT_DATA displacement;  /* storing area for boundary conditions related to displacement */
  INPUT_DATA bed;       /* designation of the bed */
  INPUT_DATA *trans;        /* the transport data */
  INPUT_DATA *sed;          /* sediment input data */
  INPUT_DATA Qs;        /* Rdiation for Heat GSAVANT */
  INPUT_DATA Td;       /* Dew POint Temp GSAVANT */
  FRICTION fterms;      /* friction terms */
  double roughness;     /* the roughness for overland or channel flow */
  double conveyance;        /* the total conveyance for edge string discharge */
  double height;        /* March 18 2004 roughness height */
  double ref_height;        /* reference height (NS hydrostatic pressure) */
  double ref_pressure;      /* reference pressure (NS hydrostatic pressure) */
  int ps_flag;          /* the post processing flag */
  int string_type;      /* tells the type of the string (NDS, EGS, FCS, etc) */
  int ice_string;
  int phys_flag;        /* flag for setting strings and secondary physics */
  int flux_flag;        /* flag for printing flux values across string */
  int sed_div_flag;     /* flag for defining a sediment diversion flux adjustment (for 2D only) */
  int link;         /* links one string to another */
  int weir_num;       /* weir Associated with a string */
  int flap_num;     /* Flp gate associated with a string*/
  int sluice_num;   /* Sluice gate associated wiht a string */
  double total_area;    /* total area of all 2d faces associated with this string */
} STR_VALUE;            /* string values */

/*********************************************************/
/* struct methods -------------------------------------- */

void sstr_value_init(STR_VALUE **, int, int, int);
void sstr_value_free(STR_VALUE *, int, int, int);
void sstr_value_check(STR_VALUE, int);
void sstr_value_checkall(STR_VALUE *, int);
void sstr_value_set_physics(STR_VALUE *str, int nstring, SFLAGS flag);
void sstr_value_set_roughness(STR_VALUE *, int, int, int,  double, double);
//void sstr_value_set_total_area(STR_VALUE *, SGRID *);
void initialize_friction_type(FRICTION *);
void initialize_input_data_type(INPUT_DATA *);

#endif

