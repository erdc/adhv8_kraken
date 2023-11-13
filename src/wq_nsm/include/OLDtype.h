/* this is the standard structure file for the ADH model */

typedef struct
{
  double value;			/* data value */
  int rank;			/* data rank  */
} DOUBLE_INT;			/* double, integer structure */

typedef struct
{
  int next;			/* the next member */
  int data;			/* the data */
} INT_LIST;			/* linked list element for integer data */

typedef struct
{
  double x, y, z;		/* the coordinates of the vector */
} VECT;				/* a vector quantity */

typedef struct
{
	double con0, con1, con2;
} CON_RES;
typedef struct
{
  double x_eq, y_eq, z_eq, c_eq;	/* the particular equation x,y,z, or continuity */
} DOF_4;

typedef struct
{
  double x_eq, y_eq, c_eq;	/* the particular equation x,y or continuity */
} DOF_3;

typedef struct
{
  double x_eq, y_eq;		/* the particular equation x or y */
} DOF_2;

typedef struct
{
  double xx, xy, xz, yy, yz, zz;	/* the components of the tensor */
} TENSOR;

typedef struct
{
  double xx, xy, yy;		/* the components of the tensor */
} TENSOR2D;

typedef struct
{
  double xx, xy, yx, yy;            /* the components of the tensor */
} TENSOR2DAI;

typedef struct
{
  double x, y;			/* the coordinates of the vector */
} VECT2D;			/* a vector quantity */

typedef struct
{
  int item0, item1, item2;
} TUPLE_3;

typedef struct
{
  double djac;			/* the jacobian of the element */
  VECT grad_shp[NDPRELM];	/* the gradients of the shape functions */
  int nodes[NDPRELM];		/* the nodes in the element */
  int levels[NDPRELM];		/* the node levels in the element */
  int global_num;
  int levels_forced[NDPRELM];   /* level of secondary refinement.  I think this is no longer needed */
} ELEM_3D;			/* the components of an element */

typedef struct
{
  double djac;			/* the jacobian of the element */
  double area3d;		/* the 3d area of the face */
  VECT2D grad_shp[NDPRFC];	/* the gradients of the shape functions */
  VECT nrml;			/* the normal to the face */
  CON_RES node_res[NDPRFC];
  int nodes[NDPRFC];		/* the nodes in the element */
  int levels[NDPRFC];		/* the node levels in the element */
  int global_num;
  int elem3d;			/* this is the 3d element that has a face coincident with this 2d element */
  int levels_forced[NDPRFC];    /* level of secondary refinement.  I think this is no longer needed */
  int sd;
} ELEM_2D;			/* the components of an element */

typedef struct
{
  double djac;			/* the jacobian of the element */
  double grad_shp[NDPREDG];	/* the gradients of the shape functions */
  int nodes[NDPREDG];		/* the nodes in the element */
  int levels[NDPREDG];		/* the node levels in the element */
  int elem2d;			/* this is the 2d element that has an edge that matches this 1d element */
  VECT2D nrml;			/* this is the 2d outward normal to the element this is connected to,  */
  /* it is in the x-y plane */
} ELEM_1D;			/* the components of an element */

typedef struct
{
  TENSOR ev;			/* the eddy viscosities */
  double refine_tolerance, unrefine_tolerance;	/* the refinement tolerances */
  int equation_type;		/* is equation type NONCONSERVATIVE or CONSERVATIVE */
  int tur_model_hor;            /* The utilized turbulence model GSavant */
  int tur_model_ver;
  double smag_coeff;
} NAV_STOKES_MAT;

typedef struct
{
  TENSOR2D ev;			/* the eddy viscosities for a 2d shallow water model */
  double refine_tolerance, unrefine_tolerance;	/* the refinement tolerances */
  int eev_flag;					/* flag inidicating Y/N for elemental determination of viscosity */
  int bed_disp_flag;                 		/* flag inidicating Y/N for sediment bed displacement */
  int wind_flag;                                /* flag for wind variable calculations (cjt) */
  int vor_flag;                                 /* flag indicating Y/N for vorticity calculations */
  double eev_coef;				/* coefficent used in viscosity calculation if eev_flag = ON */
  double windatt;				/* wind attenuation */
  double fraction;      			/* Used with FRC flag GSAVANT */
  int EEV_MODE;                 		/* EEV Computation Flag */
  int d_flag;                                   /* INdicate whether density coupling is off */
} SHALLOW_WATER_MAT;

typedef struct
{
  TENSOR k;			/* the hydraulic conductivities */
  double s_s;			/* the storage coefficient */
  double porosity;		/* porosity */
  double water_vol;		/* the volume of water of a volumetric source */
  double refine_tolerance;	/* the refinement tolerance */
  double residual_sat;		/* residual saturation */
  double vangen_alpha;		/* alpha from vangenuchten equations */
  double vangen_n;		/* vangenuchten exponent (assume mualem) */
  double brooks_lambda;		/* exponent for brooks-corey */
  double brooks_pd;		/* displacement pressure for brooks corey */
  int ikr;			/* index to the pressure - relative conductivity series */
  int isat;			/* index to the pressure - saturation curve */
  double tortuosity;		/* tortuosity (only used for diffusion transport) */
  double d_l;			/* longitudinal dispersivity */
  double d_t;			/* transverse dispersivity */
} GROUNDWATER_MAT;

typedef struct
{
  TENSOR thk;			/* the thermal conductivities */
  double sh;			/* the specific heat of the solid */
  double sg;			/* the specific gravity of the solid material (not bulk) */
  double ht_source;		/* the heat source */
  double refine_tolerance;	/* the refinement tolerance */
  double porosity;		/* porosity */
  double tortuosity;		/* tortuosity (only used for diffusion transport) */
  double d_l;			/* longitudinal dispersivity */
  double d_t;			/* transverse dispersivity */
  double albedo;		/* albedo for reflection of sw radiation */
  double emissivity;		/* emissivity */
  double thermal_model_param1;	/* first parameter in the thermal mixture model */
  double thermal_model_param2;	/* second parameter in the thermal mixture model */
  int iq;			/* index to the well temperature series */
  int thermal_model_flag;	/* flag used to select among the thermal mixture models */
} HEAT_TRANSPORT_MAT;

typedef struct
{
  double d_m;			/* the molecular diffusion coefficient */
  double source;		/* the strength of the volumetric source */
  double d_l, d_t;		/* the longitudinal and transverse dispersion coefficients */
  double *react;		/* linear reaction coefficient */
  double rd;			/* retardation coefficient */
  double tortuosity;		/* tortuosity */
  double refine_tolerance;	/* the refinement tolerance */
} TRANSPORT_MAT;

typedef struct
{
  double range;                 /* marsh porosity range */
  double kappamin;              /* marsh porosity minimum area coefficient */
  double offset;                /* marsh porosity elevation offset */
  int offset_flag;              /* flag to determine type of elevation offset (absolute or relative) */
  int mp_flag;                  /* flag to turn marsh porosity on or off */
} MARSH_POROSITY;

typedef struct
{
  double alpha0,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6;
  double beta1,beta2,beta3,beta4;
  double K1,K2,K3,K4;
  double mu, rho;
  double sigma1,sigma2,sigma3,sigma4,sigma5;
  double Kl,Kn,Kp;
  double Pn;
  double lambda0,lambda1,lambda2;
  double min_do_conc;
} NSM_PAR;

typedef struct
{
  int max_lev;			/* the maximum number of refinement levels in the material */
  NAV_STOKES_MAT nsm;		/* navier stokes materials for eddy viscosities */
  GROUNDWATER_MAT gwm;		/* material structure for groundwater flow */
  HEAT_TRANSPORT_MAT htm;	/* material structure for heat transport */
  SHALLOW_WATER_MAT swm;	/* shallow water materials for eddy viscosities */
  TRANSPORT_MAT *tm;		/* pointer to transport materials */
  MARSH_POROSITY marshpor;      /* marsh porosity terms */
  double coriolis;		/* Coriolis effect */
  NSM_PAR nps;
} MATERIAL;

typedef struct
{
  double thk;			/* thermal conductivity (assumed isotropic for fluids) */
  double sh;			/* specific heat of the solid */
  double viscosity;		/* viscosity of the fluid */
  double sg;			/* specific gravity of the fluid */
} FLUID_PROP;

typedef struct
{
  FLUID_PROP water;		/* properties for water */
  FLUID_PROP gas;		/* properties for gas (air) */
} FLUID_MAT;

typedef struct
{
  int type;			/* type of transported constituent */
  double property[NPRCN];	/* properties of the consituent like particle diameter, etc. */
} CONSTITUENT;

typedef struct
{
  double manningsn;		/* Manning's n roughness */
  double skfricoef;		/* skin friction coefficient */
  double dragcoef;		/* drag coefficient */
  double eqrheight;		/* equivalent roughness height */
  double rheight;		/* roughness height (K-E) */
  double hghtstem;		/* undeflected stem height */
  double diamstem;		/* stem diameter */
  double densstem;		/* stem density */
  double icethick;		/* thickness of the ice */
  double icedense;		/* density of the ice */
  double bedrhght;		/* bed roughness height */
  double icerhght;		/* ice roughness height */
  int mng_flag;			/* flags YES if Manning's n roughness is input */
  int erh_flag;			/* flags YES if roughness height is input */
  int sav_flag;			/* flags YES for submerged aquatic vegetation */
  int urv_flag;			/* flags YES for unsubmerged rigid vegetation */
  int icemors;			/* 0 = stationary ice, 1 = moving ice */
} FRICTION;

typedef struct
{
  int bc_flag;			/* the type of boundary condition */
  int src_flag;			/* the type of source */
  int isigma;			/* flux specification */
  int iu_0;			/* solution specification */
  int iq;			/* the strength of the source */
  int ici;			/* concentration of the source */
  int ivx;			/* the boundary specification series number for x */
  int ivy;			/* the boundary specification series number for y */
  int ivz;			/* the boundary specification series number for z */
  double ceq_0;			/* the boundary equilibrium concentration */
} INPUT_DATA;			/* input data structure */

typedef struct
{
  int hyd_flag; /* type of hydraulic structure */
  int up;       /* the upstream string number */
  int down;     /* the downstream string number */
  int iflow;    /* the series number with the flow */
  double value; /* the concentration for each constitute */
} HYDRAULIC_STRUCTURE;

typedef struct
{
  int nodeu,egsu,upnode;                                    /* upstream node and edge string */
  int noded,egsd,dnode;                                    /* downstream node and edge string*/
  double cd,z1,z2,w1,w2,sl1,sl2,hw;
  double fluxu,fluxd;                           /* flux through edge string */
}WEIR_STRUCT_C; 

typedef struct
{
  int nodeu,egsu,upnode;                                    /* upstream node and edge string */
  int noded,egsd,dnode;                                    /* downstream node and edge string*/
  double A, B, C, D, E, F;
  double GG;   /* Length of Flap Gate*/
  double fluxu,fluxd;                           /* flux through edge string */
  double density, length, width;
  int calc;
}FLAP_STRUCT_C; 


typedef struct
{
  INPUT_DATA flow;		/* the flow data */
  INPUT_DATA heat;		/* the thermal data */
  INPUT_DATA ol_flow;		/* the overland flow data */
  INPUT_DATA ch_flow;		/* the channel flow data */
  INPUT_DATA pressure;		/* the pressure data */
  INPUT_DATA displacement;	/* storing area for boundary conditions related to displacement */
  INPUT_DATA bed;		/* designation of the bed */
  INPUT_DATA *trans;		/* the transport data */
  FRICTION fterms;		/* friction terms */
  double roughness;		/* the roughness for overland or channel flow */
  double conveyance;		/* the total conveyance for edge string discharge */
  double height;		/* March 18 2004 roughness height */
  double ref_height;		/* reference height (NS hydrostatic pressure) */
  double ref_pressure;		/* reference pressure (NS hydrostatic pressure) */
  int ps_flag;			/* the post processing flag */
  int string_type;		/* tells the type of the string */
  int ice_string;
  int phys_flag;		/* flag for setting strings and secondary physics */
  int flux_flag;		/* flag for printing flux values across string */
  int sed_div_flag;		/* flag for defining a sediment diversion flux adjustment (for 2D only) */
  int link;			/* links one string to another */
  int weir_num;       /* weir Associated with a string */
  int flap_num;     /* Flp gate associated with a string*/
  double out_flow;  /* Specified outflow */
  double in_flow;
  int dpl_brch;    /* Nodal Displacement due to breach  GSAVANT */
  int breach_flag;
  double total_area;
} STR_VALUE;			/* string values */

typedef struct
{
  double time;			/* time */
  double value;			/* value */
  double y_value;		/* direction */
  double slope;			/* slope */
} XY_ITM;			/* a series data point */

typedef struct
{
  int size;			/* # of items in the XY series */
  int type;			/* the type of series */
  double infact;		/* conversion factor to seconds */
  double outfact;		/* conversion factor from seconds */
  double *start;		/* start time */
  double *end;			/* end time */
  double *dt;			/* interval */
  double value;			/* current value of the interpolant */
  double x_coord;		/* coordinate of the wind station for the xy series */
  double y_coord;
  int flag;                     /* flag for shears or x,y speeds (WU or TEETER) for wind */
  double x_wind;		/* wind term in first column (x directions or speeds )*/
  double y_wind;                /* wind term in second colum (y directions or speeds */
  double x_stress;
  double y_stress;
  double tol;			/* the tolerance of the series */
  XY_ITM *sers;			/* the XY series */
} XY_SERS;

typedef struct
{
  int julian;			/* julian day */
  double time;			/* time */
  double air_temp;		/* air temperature, C */
  double rel_humidity;		/* relative humidity */
  double pressure;		/* barometric pressure UNITS */
  double wind;			/* windspeed x meters above ground */
  double solar;			/* incoming shortwave solar radiation, J/cm2-s */
  double longwave;		/* downwelling longwave radiation, J/cm2-s */
  double precip;		/* precipitation */
  double zenith;		/* solar zenith */
  double azimuth;		/* solar azimuth */
} MET_ITM;			/* structure for an individual met data entry */

typedef struct
{
  int size;			/* # of time entries in the met data */
  int interp_type;		/* can use piecewise linears or step interpolation (later) */
  MET_ITM *entry;		/* the met data entries */
} MET_DATA;

typedef struct
{
  double *value;		/* the values of the sparse vector */
  int *index;			/* the indices of the sparse vector */
  int size;			/* the number of entries in the sparse vector */
  int max_size;			/* the number of allocated entries in the sparse vector */
} SPARSE_VECT;			/* a sparse vector */

typedef struct
{
  int dm;			/* flag tells if the dm card has been read */
  int tvs;			/* flag tells if the tvs card has been read */
  int dpl;			/* flag tells if the dpl card has been read */
  int dpt;			/* flag tells if the dpt card has been read */
  int *rct;			/* flag tells if the rct card has been read */
  int rd;			/* flag tells if the rd card has been read */
  int tor;			/* flag tells if the tor card has been read */
  int trt;			/* flag tells if the trt card has been read */
} TRN_MAT_FLAGS;		/* transport material flags */

typedef struct
{
  int type;			/* type (broad-crested, ogee, etc.) */
  double hght;			/* height w.r.t. local elevation */
  double wdth;			/* width */
  double len;			/* length */
  double cyoo;			/* flow across the crest */
  double cdis;			/* user input coefficient of discharge */
  int segs;			/* number of 1D elements along the weir */
  int ustr;			/* associated upstream node string */
  int dstr;			/* associated downstream node string */
} WEIR;				/* it's a weir! */

typedef struct
{
  int max_lev;			/* flag tells if the max level card has been read */
  int ss;			/* flag tells if the ss card has been read */
  int sat;			/* flag tells if the sat card has been read */
  int k;			/* flag tells if the k card has been read */
  int kr;			/* flag tells if the kr card has been read */
  int ev;			/* flag tells if the ev card has been read */
  int evs;			/* flag tells if the evs card has been read */
  int por;			/* flag tells if the por card has been read */
  int fvs;			/* flag tells if the fvs card has been read */
  int frt;			/* flag tells if the frt card has been read */
  int nce;			/* flag tells if the nce card has been read */
  int nrt;			/* flag tells if the nrt card has been read */
  int srt;			/* flag tells if the srt card has been read */
  TRN_MAT_FLAGS *trans;		/* transport material flags */
} MAT_FLAGS;			/* material flags */

typedef struct
{
  int flag[MATFLAG_MAX];	/* flag for each material element */
  int **trans;			/* matrix for flags for each transport element */
  int **rct;			/* matrix for flags for each reaction element */
} MAT_ELEM;			/* keeps flags for material element */

/* the TRANSPORT_IDENTITY will contain an indicator for each specific variable you are to route */
/* if it is a general constituent you won't need to add it here */
typedef struct
{
  int tke;			/* turbulent kinetic energy variable */
  int tds;			/* turbulent dissipation */
  int vor;			/* vorticity */
  int tmp;			/* temperature */
  int sal;			/* salinity */
} TRANSPORT_IDENTITY;

struct FACE_LIST_ITEM
{
  int nd1;			/* the first node on the face */
  int nd2;			/* the second node on the face */
  int nd3;			/* the third node on the face */
  int ie1;			/* the first element on the face */
  int ie2;			/* the second element on the face */
  struct FACE_LIST_ITEM *next;	/* the next item in the linked list */
};
typedef struct FACE_LIST_ITEM FACE_LIST_ITEM;	/* special line for linked list structure - refers to itself */

struct ELEM3D_LIST_ITEM
{
  int nd1;			/* the first node */
  int nd2;			/* the second node */
  int nd3;			/* the third node */
  int nd4;			/* the fourth node */
  int ielem;			/* the element number */
  struct ELEM3D_LIST_ITEM *next;	/* the next item in the linked list */
};				/* a linked list item for 3D elements */
typedef struct ELEM3D_LIST_ITEM ELEM3D_LIST_ITEM;	/* special line for linked list structure - refers to itself */

struct ELEM2D_LIST_ITEM
{
  int nd1;			/* the first node */
  int nd2;			/* the second node */
  int nd3;			/* the third node */
  int ielem;			/* the element number */
  struct ELEM2D_LIST_ITEM *next;	/* the next item in the linked list */
};				/* a linked list item for 2D elements */
typedef struct ELEM2D_LIST_ITEM ELEM2D_LIST_ITEM;	/* special line for linked list structure - refers to itself */

struct ELEM1D_LIST_ITEM
{
  int nd1;			/* the first node */
  int nd2;			/* the second node */
  int ielem;			/* the element number */
  struct ELEM1D_LIST_ITEM *next;	/* the next item in the linked list */
};				/* a linked list item for 1D elements */
typedef struct ELEM1D_LIST_ITEM ELEM1D_LIST_ITEM;	/* special line for linked list structure - refers to itself */

typedef struct
{
  int begin;			/* the beginning of the band */
  int end;			/* the end of the band */
  int size;			/* the allocated size of the band */
  double *value;		/* the values */
} BAND_VECT;			/* a band storage sparse vector */

struct EDGE_LIST_ITEM
{
  int nd1;			/* the first node on the edge */
  int nd2;			/* the second node on the edge */
  int new_node;			/* the new node */
  int rank;			/* the rank of the edge */
  struct EDGE_LIST_ITEM *next;	/* the next item in the linked list */
};				/* a linked list item for edges */
typedef struct EDGE_LIST_ITEM EDGE_LIST_ITEM;	/* special line for linked list structure - refers to itself */

/* structure used to rank the edges */
typedef struct
{
  int number;			/* which edge */
  double length;		/* the length of the edge */
} EDGE_RANK;			/* structure to sort the edges */

struct NODE_LIST_ITEM
{
  struct NODE_LIST_ITEM *next;	/* the next item in the linked list */
  int rnode;			/* relative node number on owning processor */
  int sd;			/* owning processor */
  int local;			/* local node number on myid */
};				/* a linked list item for global node numbers */
typedef struct NODE_LIST_ITEM NODE_LIST_ITEM;	/* special line for linked list structure - refers to itself */

struct ELEM_REF_LIST_ITEM
{
  struct ELEM_REF_LIST_ITEM *next;	/* the next item in the linked list */
  int ielem;			/* the element number */
};				/* a linked list item for element numbers */
typedef struct ELEM_REF_LIST_ITEM ELEM_REF_LIST_ITEM;	/* special line for linked list structure - refers to itself */

/*!
   Every Send/Receive uses a "message buffer" structure that stores info about 
   the send/receive. 

   size - Allocated Size in Bytes of *buffer. 
   It should only be non-zero if *buffer is specifically allocated with its 
   "own" memory. For instance, in a send, it may be necessary
   to gather pieces of non-contiguous information. In that case the *buffer
   is allocated its own memory, pieces are packed into it, in some order, and then 
   a contiguous chunk of information is ready to send. If *buffer is acting as 
   a simple pointer to (other) contiguous memory, then size needs to be zero. Any 
   previously allocated memory to *buffer should be freed first. 

   nitem - Number of items being sent/received. Needed by the MPI call.

   type - Type of data being sent (MPI_INT, MPI_DOUBLE, etc.)

   sd - The "other" processor with which processor myid is communicating with

   pos - An offset within *buffer, used for packing/unpacking mixed-type messages
   where (MPI_PACK) is called to pack the buffer.

   *buffer - Is the buffer that will be sent/received. It is a (void *) because
   many different types of information may be sent (int, double, packed/mixed, etc.) 
   It may be its own allocated space, or it may be a pointer to some other 
   (separately allocated memory).
 */
typedef struct
{
  int size;			/* the size of the buffer - in bytes */
  int nitem;			/* the number of items in the buffer */
  int type;			/* the type of buffer */
  int sd;			/* the processor the message is being send to */
  int pos;			/* current position in the buffer */
  void *buffer;			/* the buffer */
} MESSG_BUFFER;			/* message buffer */

typedef struct
{
  int size;			/* the size of the key */
  int *key;			/* the key */
} MESSG_KEY;			/* message key */

typedef struct
{
  int rnode;			/* the relative node number - local node number on 
				   owning processor */
  int sd;			/* the owning processor */
  int global_num;		/* global node number */
} GLOBAL_NODE;			/* the global node number */

typedef struct
{
  double shpn[3];
  int nd[3];
} CSTRCT;

typedef struct
{
  int num;
  int num_spnt;
  int *spnt;
  int *indx;
} WCOL;

typedef struct __LIGHT
{
  int alpha;  /* series containing angles between x-axis and light ray */
  int theta;  /* series containing angles between x-y plane and light ray */
  int intensity; /* light intensity */
  double min_intens;  /* minimum light intensity */
  double attenuation; /* light attenuation */
} LIGHT;

typedef struct __ICM_INPUT
{
  double *kback;
  double *kfs;
  double *koc;
  double *wss;
  double *wsl;
  double *wsr;
  double *wsa[3];
  double kadpo4;
  double kadtox;
} ICM_INPUT;

typedef struct __ICM_DATA
{
  int ntags;			/* number of ICM element groups */
  int active[NICMVARS];		/* indicates which ICM variables are activated */
  int c_link[NICMVARS];		/* links ICM vars. to ADH constituents */
  int allnodes[NICMVARS];	/* set initial value to all nodes */
  int *node_icm;		/* ICM element group tag (nnode) */
  double run_length;		/* tfinal - tinitial */
  double vdata[NICMVARS];	/* nodal variable data */
  double *init[NICMVARS];	/* ICM varible initial value (ntags) */
  double *derive_q;		/* quantities derived by ICM (nnode * NICMDRQ) */
  double *dcdt;			/* reaction rates (nnode * NICMVARS) */
} ICM_DATA;

typedef struct
{
  char sdate[MAXLINE];		/* the initial date (days) for the DSS time window */
  char stime[MAXLINE];		/* the initial time (minutes) for the DSS time window */
  char edate[MAXLINE];		/* the final date (days) for the DSS time window */
  char etime[MAXLINE];		/* the final time (minutes) for the DSS time window */
  int tinit;			/* the initial time (seconds) for the DSS time window */
  int tfin;			/* the final time (seconds) for the DSS time window */
  int nsers;			/* the number of DSS series */
  int node_str;			/* the node string containing the output stations */
  int nstation;			/* the number of nodes in the station node string */
  int itime;			/* the current step in time (0 - ntimes-1) */
  int *stnodes;			/* the nodes of the output stations */
  int *times;			/* output times (output_sers.size) */
  float **values;		/* values associated with each time (nstation * output_sers.size) */
} HEC_DATA;

/* Structures for SW3 COLUMNS */

struct CENT_LIST_ITEM
{
  VECT2D vect;   		/* the coordinates of the point */
  int index;                    /* an index for this column (e.g., column id) */
  struct CENT_LIST_ITEM *next;	/* the next item in the linked list */
};				/* a linked list item for centroids (points) */
typedef struct CENT_LIST_ITEM CENT_LIST_ITEM;	/* special line for linked list structure - refers to itself */

struct ID_LIST_ITEM
{
  int id;			/* the id of the member */
  struct ID_LIST_ITEM *next;	/* the next item in the linked list */
  struct ID_LIST_ITEM *prev;	/* the next item in the linked list */
};
typedef struct ID_LIST_ITEM ID_LIST_ITEM;	/* special line for linked list structure - refers to itself */

struct MIDPT_LIST_ITEM
{
  VECT vect;                    /* the coordinates of the point */
  int index;                    /* an index for this column (e.g., column id) */
  int node1;               /* this midpoint lies on the edge formed by
                              node1 and node2 */
  int node2;
  int surf_node1;          /* the two surface nodes that lie above node1 and node2 */
  int surf_node2;
  int vertical;            /* 0 = not a vertical edge; 1 = this is a vertical edge */
  int column1;             /* the indices of the 2 columns that share this edge */
  int column2;

  int elem_upper[2];          /* the element that lies above this edge/midpt in each of
                                 the two columns */
  int elem_lower[2];          /* the element that lies below this edge/midpt in each of
                                 the two columns */

  double value[5];            /* data associated with this midpoint */
                              /* value[0] = pressure at midpt;
                               * value[1] = pressure from positive perturbation at the "left" surface node
                               * value[2] = pressure from positive perturbation at the "right" surface node
                               * value[3] = pressure from negative perturbation at the "left" surface node
                               * value[4] = pressure from negative perturbation at the "right" surface node
                               */
  struct MIDPT_LIST_ITEM *next;	/* the next item in the linked list */
  struct MIDPT_LIST_ITEM *prev;	/* the next item in the linked list */
};				/* a linked list item for centroids (points) */
typedef struct MIDPT_LIST_ITEM MIDPT_LIST_ITEM;	

typedef struct edges_entry
{
    struct edges_entry *next;
      int edge_number;
        int nodes[2];
} EDGES_ENTRY;

typedef struct
{
  int nedges;                   /* Number of Edges in the 2D mesh on this PE */
  int num_columns;
  int num_surf_nodes;
  int num_midpts;
  int *surface_elements;       /* An array of the surface element id for each
                                  column */
  int *surface_elements2d;     /* An array of the 2D surface element id for each
                                                                    column */
  int *bottom_elements2d;     /* An array of the bottom 2D element id for each
                                                                    column */
  int *bottom_nodes;         /* a flag indicating whether node is bottom or not cjt */
  int *bottom_elems;         /* a flag indicating whether the elem is on bottom or not cjt */
  CENT_LIST_ITEM **vertical_hash; /* Hash table for xy-coords of points;
                                       used to build list of nodes that lie
                                       on a vertical line */
  ID_LIST_ITEM **vertical_list;   /* Linked list of nodes that lie on
                                       vertical segment */
  CENT_LIST_ITEM **column_hash;   /* Hash table for columns of 3D elements */
  ID_LIST_ITEM **column_list;     /* Linked list of 3D elements that lie in columns */
  CENT_LIST_ITEM **column_hash2d;   /* Hash table for columns of 2D elements */
  ID_LIST_ITEM **column_list2d;     /* Linked list of 2D elements that lie in columns */
  CENT_LIST_ITEM **midpt_hash;
  MIDPT_LIST_ITEM **midpt_list;
  ID_LIST_ITEM **sidewall_list;    /* Linked list of 2D sidewall elements */

  EDGES_ENTRY **edges_hash;       /* Hash table for edges */

  ELEM_1D *edges;                 /* Edges Data Structure List */
} Mesh_Info;

/* structs for wind and wave file reading (cjt) */
typedef struct {
  int n;                        /* total number of wind nodes */
  double dt;                    /* time between wind snaps */
  double tprev;                 /* previous wind time */
  double tnext;                 /* next wind time */
  FILE * fin;                   /* file pointer for wind reads */
} WIND_INFO;

/* struct that stores tidal series info (cjt) */
typedef struct {
        int total_bcs;                  /* total number of xy-series to create */
        int xyseries;                   /* the xy-series associated with this edge string */
        int edge_string;                /* edge string associated with tidal series */
        int nlines;                     /* number of lines in tidal series */
        double to;                      /* initial time of tidal series */
        double tf;                      /* final time of tidal series */
        double dt;                      /* tidal series time increment */
        double *eta;                    /* vector[nlines] of water surface elevations */
        double *time;                   /* vector[nlines] of times */
} TIDAL_SERIES;

/* struct to store nodal short-wave information (cjt) */
/* note :: if anything is added to this, the following must be changed/added to:
        include/define.h (NWVAR)
        winds/wave_routines.c
        ps_print_ts_swave.c
*/
typedef struct {
        double Sxx;             /* xx - wave stress */
        double Sxy;             /* xy - wave stress */
        double Syy;             /* yy - wave stress */
        double height;          /* wave height (either rms or mo) */
        double period;          /* wave period */
        double angle;           /* wave angle (radians) */
        double number;          /* wave number */
        double breaking;        /* describes the degree of wave breaking */
        double speed;           /* wave phase speed */
        double ediss_break;     /* wave energy dissipation from breaking */
        double ediss_frict;     /* wave energy dissipation from friction */
        double TauX;            /* wave force -> (-dSxx/dx - dSxy/dy) */
        double TauY;            /* wave force -> (-dSyy/dy - dSyx/dx) */
} SWAVE;
