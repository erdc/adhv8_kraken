/* constants */
#define EARTH_ROTAT 7.2722E-5	/* Earth's rotation in 1/s (coriolis: added 6-02) */
#define DEGREE2RAD 0.017453292519943
#define PI 3.141592653589793238462643383279502884197169399375105820974944592308
#define SQUAREROOTTWO 1.4142135623730950488016887242096980785696718753769480732

#define DZERO 0.0e+0
#define D8_TH 1.25e-1
#define DQUAR 2.5e-1
#define DHALF 0.5e+0
#define DONE 1.0e+0
#define DTWO 2.0e+0
#define DFOUR 4.0e+0
#define NONE -1.0e+0
#define SMALL2 1.0e-2
#define SMALL3 1.0e-3
#define SMALL4 1.0e-4
#define SMALL5 1.0e-5
#define SMALL6 1.0e-6
#define SMALL7 1.0e-7
#define SMALL8 1.0e-8
#define SMALL10 1.0e-10
#define SMALL12 1.0e-12
#define SMALL14 1.0e-14
#define SMALL16 1.0e-16
#define EPSILON 1.0e-3
#define SMALL DBL_EPSILON	/* the smallest double */
#define NOT_QUITE_SMALL FLT_EPSILON	/* the smallest single precision number */

/* FOR GRID READ */
#define NEWIO_HASHSIZE 10069

#define one_2  0.5e+0
#define one_3  0.333333333333333333333333
#define one_4  2.5e-1
#define one_5  0.200000000000000000000000
#define one_6  0.166666666666666666666667
#define one_8  0.125000000000000000000000
#define one_9  0.111111111111111111111111
#define one_12 0.083333333333333333333333
#define one_18 0.055555555555555555555556
#define one_20 0.050000000000000000000000
#define one_24 0.041666666666666666666667
#define one_30 0.033333333333333333333333
#define one_36 0.027777777777777800000000
#define one_45 0.022222222222222222222222
#define one_48 0.020833333333333332000000
#define one_60 0.016666666666666666666667
#define one_72 0.013888888888888888888888
#define one_120 0.00833333333333333000000
#define one_144 0.00694444444444444444444
#define one_180 0.00555555555555555555556
#define one_360 0.00277777777777778000000
#define one_720 0.00138888888888889000000
#define one_1440 0.0006944444444444440000


/* SEDIMENT ALGORITHM CONSTANTS FOR READING PARAMETERS FROM BC FILE */
#define COHSET 4
#define WINDWAVE 10
#define NONCOENT 1 /* (cjt) -- mdw SW3_FLOW svn merge */
#define DIVCOEF 3

/* input line definitions */
#define MAXLINE 300		/* the maximum number of characters on an input line */

/* geometric definitions */
#define TRIANGLE 1
#define QUADRILATERAL 2
#define TETRAHEDRON 3
#define PRISM 4

// type of grid
#define UNSTRUCTURED 2
#define COLUMNAR 3

#define BODY -1
#define SURFACE 0
#define BED 1
#define SIDEWALL 2
#define COLUMN 3
#define BODY2D 4

#define NDONSEG 2       /* nodes on line segment */
#define NDONTRI 3       /* nodes on triangle */
#define NDONTRIQUAD 6
#define NDONQUAD 4       /* nodes on quadrilateral */
#define NDONQUADQUAD 8
#define NDONTET 4         /* nodes on tetrahedron */
#define NDONPRISM 6       /* nodes on prism */
#define NDONPYRAMID 5     /* nodes on a pyramid */
#define NDONTETQUAD 10    /* nodes on a tetrahedron with midpoints for quadratic baasis */
#define NDONPRISMQUAD 15  /* nodes on a triangular prism with midpoints for quadratic baasis */

#define MAX_NNODES_ON_ELEM1D 2  /* max number of nodes on all possible 2d elements (line here) */
#define MAX_NNODES_ON_ELEM2D 4  /* max number of nodes on all possible 2d elements (quads here) */
#define MAX_NNODES_ON_ELEM3D 6  /* max number of nodes on all possible 3d elements (triprisms here) */
#define MAX_NNODES_ON_ELEM3D_QUAD 15  /* max number of quadratic nodes on all possible 3d elements (triprisms here) */

#define MAX_QUAD_ORDER 4 /* the max integrand order for quadrature */
#define MAX_QUAD_POINTS 18 /* max number of quadrature points for an element */

#define NFCPRELM 4		/* the number of faces per element */
#define NDPRELM 4		/* the number of nodes per element */
#define NDPRFC 3		/* the number of nodes per face */
#define NEDGEPRELM 6		/* the number of edges per element */
#define NEDGEPRFC 3		/* the number of edges per face */
#define NDPREDG 2		/* the number of nodes per edge */
#define NDPRELM_MINUS_1 3	/* the number of nodes per element -1 */
#define NDPRFC_MINUS_1 2	/* the number of nodes per face -1 */

/* matrix sizes */
#define NELEMAT 16		/* NDPRELM*NDPRELM */
#define NFCMAT 9		/* NDPRFC*NDPRFC */
#define NEDGEMAT 4		/* NEDGEPRELM*NEDGEPRELM */
#define NFCMAT2 12		/* NDPRELM*NDPRFC */
#define NPRCN 4			/* The number of properties per constituent */

/* memory increment for small static arrays */
#define MEM_INC 128

/* the amount of reduction in time step for failure of the nonlinear iteration */
#define DT_REDUCE_FACTOR 0.25	/* time step reduction for failure */
#define DT_ENLARGE_FACTOR 2.0	/* time step expansion for success */

/* the default tolerance for a series */
#define XY_DEFAULT_TOL 1.0E-1

/* Constituent numbers */
#define TKE 1			/* turbulent kinetic energy */
#define TDS 2			/* turbulent dissipation */
#define CLA 3			/* Clay sediment */
#define SLT 4			/* Silt sediment */
#define SND 5			/* Sand sediment */
#define CON 6			/* Any old constituent */
#define VOR 7			/* Vorticity Constituent */
#define TMP 8			/* Temperature Constituent */
#define SAL 9			/* Salinity Constituent */
#define NOO 10          /* WQ Constituent */
#define NOOO 11          /* WQ Constituent */
#define NHHHH 12          /* WQ Constituent */
#define OrN 13           /* WQ Constituent */
#define OP 14           /* WQ Constituent */
#define DP 15           /* WQ Constituent */
#define POOOO 16          /* WQ Constituent */
#define ALG 17          /* WQ Constituent */
#define CBOD 18         /* WQ Constituent */
#define DO 19           /* WQ Constituent */

/* Define Dimension Units */
#define MKS 1			/* Meters-Kilograms-Seconds and expects Centigrade */
#define FPS 2			/* Feet-Pound-Seconds and expects Fahrenheit */

/* other definitions */
#define NDIM 3			/* the maximum number of dimensions */
#define ON 1			/* on */
#define OFF 0			/* off */
#define YES 1			/* yes */
#define NO -3			/* no */
#define TRUE 1			/* true */
#define FALSE 0			/* false */
#define AGREE 0			/* agreement in a string comparison */
#define USED 1			/* in use */
#define UNUSED 0		/* not in use */
#define ONE 1			/* 1 for binary io */
#define SOLV_TOL DBL_EPSILON	/* the linear solver tolerance - usually single precision zero - 
				   NOTE:  this number MUST be significantly larger than SMALL */
#define TET_PSI_AT_CENTROID 0.25	/* the linear basis functions on a tet evaluated at the centroid */
#define SOLV_UPDATE_FACTOR 0.01	/* the factor for the update in the solver */
#define HASHSIZE 25013
#define CHAR_SIZE 12

/* definitions of names for various solves */
//#define HYD 1
//#define TRN 2
//#define BLT 3
//#define SLT 4
//#define GRD 5
//#define HVEL 6
//#define WVEL 7
//#define SALT 8
//#define CLAY 9
//#define SAND 10
//#define TEMP 11
//#define DIFF 12
//#define NS 13
//#define GW 14

/* post processing dataset flags */
typedef enum
  {
    UNSET_PS_FLAG = -1,    /* Unassigned flag */
    PS_FLAG_HEAD,          /* heads */
    PS_FLAG_SAT,           /* saturations */
    PS_FLAG_OLHEAD,        /* overland flow heads */
    PS_FLAG_CHHEAD,        /* channel heads */
    PS_FLAG_PRS,           /* pressures */
    PS_FLAG_VEL,           /* velocities */
    PS_FLAG_DPL,           /* grid displacements */
    PS_FLAG_BED_DPL,       /* bed displacements */
    PS_FLAG_OLVEL,         /* overland velocities */
    PS_FLAG_GSP,           /* grid speed */
    PS_FLAG_ERR,           /* error indicator */
    PS_FLAG_FLX,           /* nodal flux */
    PS_FLAG_HEAT,          /* the temperatures */
    PS_FLAG_UGD,           /* x-vel gradient */
    PS_FLAG_VGD,           /* y-vel gradient */
    PS_FLAG_WGD,           /* z-vel gradient */
    PS_FLAG_WS_POT,        /* potential (SW3) */
    PS_FLAG_VEL_COR,       /* velocity correction (SW3) */
    PS_FLAG_POLHEAD,       /* old overland flow heads */
    PS_FLAG_POLVEL,        /* old overland velocities */
    PS_FLAG_SHSTR,         /* shear stresses */
    PS_FLAG_2D_SC,         /* 2D scatter plot */
    PS_FLAG_3D_SC,         /* 3D scatter plot */
    PS_FLAG_DEPTH,         /* depth */
    PS_FLAG_PDEPTH,        /* old depth */
    PS_FLAG_BSH,           /* bed shear stress */
    PS_FLAG_VORT,          /* SW2 vorticity */
    PS_FLAG_SFLX,          /* surface heat flux */
    PS_FLAG_PRS_HOT,       /* pressure hotstart creation (NS) */
    PS_FLAG_VEL_HOT,       /* velocity hotstart creation (NS) */
    PS_FLAG_DEP_HOT,       /* depth hotstart creation (SW2) */
    PS_FLAG_PDP_HOT,       /* previous depth hotstart creation (SW2) */
    PS_FLAG_OVL_HOT,       /* velocity hotstart creation (SW2) */
    PS_FLAG_POV_HOT,       /* previous velocity hotstart creation (SW2) */
    PS_FLAG_DIS_HOT,       /* displacement hotstart creation (SW2) */
    PS_FLAG_DSP_HOT,       /* displacement hotstart creation (NS) */
    PS_FLAG_GSP_HOT,       /* grid speed hotstart creation (NS) */
    PS_FLAG_HED_HOT,       /* head hotstart creation (GW) */
    PS_FLAG_TMP_HOT,       /* temperature hotstart creation (HT) */
    PS_FLAG_SMR,           /* SW2 sed mass balance */
    PS_FLAG_ICM_VAR,       /* ICM variable */
    PS_FLAG_ICM_VAR_HOT,   /* ICM variable hotstart creation */
    PS_FLAG_ICM_DRQ,       /* ICM derived quantity */
    PS_FLAG_SWAVE_HEIGHT,  /* cjt wind-waves output (WAVE) */
    PS_FLAG_SWAVE_FORCES,  /* cjt wind-wave forces (WAVE) */
    PS_FLAG_SWIND_FORCES,  /* cjt wind forces */
    PS_FLAG_BED_ELV,       /* bed elev */
    PS_FLAG_ERR_HYDRO,     /* error indicator-hydro part */
    PS_FLAG_SUS,           /* suspended load vector */
    PS_FLAG_PDPL,          /* previous time step displacement */
    PS_FLAG_PVEL,          /* previous time step velocity */
    PS_FLAG_CON,           /* concentration (SW2 and GW) */
    PS_FLAG_CON_HOT,       /* concentration hotstart creation (SW2) */
    PS_FLAG_PCON,          /* previous concentration (SW2) */
    PS_FLAG_PCON_HOT,      /* previous concentration hotstart creation (SW2) */
    PS_FLAG_BLT,           /* bed layer thickness (SW2) */
    PS_FLAG_BLT_HOT,       /* bed layer thickness hotstart creation (SW2) */
    PS_FLAG_BLD,           /* bed layer distribution (SW2) */
    PS_FLAG_BLD_HOT,       /* bed layer distribution hotstart creation (SW2) */
    PS_FLAG_ALT,           /* active layer thickness (SW2) */
    PS_FLAG_ALT_HOT,       /* active layer thickness hotstart creation (SW2) */
    PS_FLAG_ALD,           /* active layer distribution (SW2) */
    PS_FLAG_ALD_HOT,       /* active layer distribution hotstart creation (SW2) */
    PS_FLAG_CBP,           /* Cohesive bed layer properties (SW2) */
    PS_FLAG_CBP_HOT,       /* Cohesive bed layer properties hotstart creation (SW2) */
    PS_FLAG_ERR_CON,       /* transport residual error (SW3) */
    PS_FLAG_WINDS,         /* wind stress output file */        
    PS_FLAG_WAVES,         /* wind wave stress output file */
    PS_FLAG_DAVG_VEL,      /* depth avg velocity file (sw3d) */
    PS_FLAG_SURF_VEL,      /* surface velocity file (sw3d) */
    PS_FLAG_BOTT_VEL,      /* bottom velocity file (sw3c) */
    PS_FLAG_BED,           /* sediment output file */
    PS_FLAG_BED_FLUX,       /* sediment error file */
    PS_FLAG_ALAYER,       /* sediment bedload vector */
    PS_FLAG_BED_LAYER,       /* sediment suspendend load vector */
    PS_FLAG_SL_GRAIN,
    PS_FLAG_BL_GRAIN,
	PS_FLAG_HYV,       /* Flag for Hydro viscosity GSAVANT */
	PS_FLAG_TRD,       /* Flag for Hydro diffusivity GSAVANT */
    /** groundwater **/
    PS_FLAG_PHEAD,    /* flag for the pressure heads */
    PS_FLAG_THEAD,    /* flag for the total heads */
    PS_FLAG_DENS,    /* flag for the density */
	/* Add new flags above here */
    /* This must be last */
    NUM_PS_FLAG_TYPES      /* number of ps flags */
  } PS_FLAG_TYPE;

/* types of physics flags - primarily used to set string bcs and secondary physics */
#define NS_FLAG 1		/* navier-stokes flow */
#define GW_FLAG 2		/* groundwater flow */
#define HT_FLAG 3		/* heat flow */
#define OL_FLAG 4		/* overland flow */
#define CH_FLAG 5		/* channel flow */
#define SW2_FLAG 6		/* shallow water 2D flow */
#define SW3_FLAG 7		/* shallow water 3D flow */
#define USR_FLAG 8		/* sediment on this face */
#define SED_FLAG 9		/* sediment on this face */

/* material matrices sizes and attributes (see initio/init_read_bc.c */
#define MPFLAG_MAX 8		/* number of MP flags (so far) */
#define MATFLAG_MAX 11		/* number of material type flags */
#define MISCFLAG_MAX 13		/* number of misc. flags */
#define TRANSFLAG_MAX 4		/* number of required transport type flags */
#define PHYS_TYPES 8		/* number of physics types */

/* material attributes */
#define ML  0
#define EV  1
#define EVS 2
#define FVS 3
#define FRT 4
#define NCE 5
#define NRT 6
#define SRT 7
#define EEV 8
#define DPL 9
#define DPT 10
#define SS  11
#define SAT 12
//Mark Causes issues with petsc? changing from K to K1
#define K1   13
#define KR  14
#define POR 15
#define TOR 16

/* transport material attributes */
#define DF  0
#define TRT 1
#define TVS 2
#define RD  3

/* MP flags */
#define MU  0
#define G   1
#define RHO 2
#define TMN 3
#define TCN 4
#define U0  5
#define L0  6
#define R0  7

/* Misc. flags */
#define MIT 0
#define NIT 1
#define T0  2
#define IDT 3
#define TF  4
#define INC 5
#define PRE 6
#define FTYPE 7
#define SDI 8
#define SST 9
#define NTL 10 
#define ITL 11
#define FIN 12          /* External file read (cjt) */

/* Hydraulic Structures Flag */
#define HST_FLOW_THROUGH 1

/* an unflaged normal item */
#define NORMAL -1		/* normal */

/* refinement flags */
#define REFINE 1                /* refine the element */
#define FINE 0                  /* leave the element alone */
//#define UNREFINE -1             /* unrefine the element */
#define NOUNR 0                 /* refined element */
#define CONSV 1                 /* merged element at t(n-1) */
//#define UNREF 2                 /* merged element at t(n) */ 

/* boundary condition string types */
#define BCT_DIR 1		/* a dirichlet boundary condition */
#define BCT_ROB 2		/* a robin boundary condition */
#define BCT_NEU 3		/* a neumann boundary condition */
#define BCT_OUTFLOW 4		/* an out flow boundary condition - for transport */
#define BCT_FLUX 5		/* total flux boundary conditions - for transport */
#define BCT_NO_FLUX 6		/* no flux bc */
#define BCT_PRS_DIR 7		/* Pressure Dirichlet boundary condition */
#define BCT_VEL_DIR 8		/* Velocity Dirichlet boundary condition */
#define BCT_PRS_NEU 9		/* Pressure Neumann boundary condition */
#define BCT_VEL_NEU 10		/* Velocity Neumann boundary condition */
#define BCT_BED 11		/*  This is the sediment bed and invokes several boundary conditions */
#define BCT_WATER_SOURCE 12
#define BCT_FLUX_COUPLE 13
#define BCT_CEQ 36		/* Equilibrium Concentration boundary condition */
#define BCT_HSP_DIR 41          /* hydrostatic pressure boundary */
#define BCT_HYBRID_INTERNAL 37
#define BCT_HYBRID_EXTERNAL 38  /* gkc 2d-3d FLUX coupling, for flagging interface edges/surfaces */

// overland flow
#define BCT_ZERO_DEPTH_GRAD 50
#define BCT_CRITICAL_DEPTH 51

/* well types */
#define EXTRACTION_WELL 12	/* an extraction well */
#define INJECTION_WELL 13	/* an injection well */

#define BCT_FREE_DIR 14		/* set free surface boundary condition, pressure is defined */
#define BCT_DPL_DIR  15		/* this means that the displacement is read/defined */
#define BCT_VEL_PRS_DIR 16	/* Velocity and pressure Dirichlet boundary condition */

#define BCT_DIR_INT  17		/* a dirichlet boundary condition that is calculated rather than read */
#define BCT_OVH_NEU  18		/* this is a boundary that uses head all the time and momentum flux */
#define BCT_FRS 19		/* this is the free surface used with face fluxes in SW3 */
#define BCT_LID_DFT 20		/* the draft is set via pressure in SW2 */
#define BCT_LID_ELV 21		/* the culvert elevation is set in SW2 */
#define BCT_LID_DEP 22		/* the culvert depth is set in SW2 */
#define BCT_SPL_NEU 23		/* the velocity out the model is set to the critical velocity */
#define BCT_DIS_NEU 24		/* Total Discharge Neumann boundary condition */
#define BCT_OUT_NEU 25      /* Outflow from inside the mesh */
#define BCT_IN_NEU 26       /* Inflow linked to outflow */
#define BCT_WRSU 27       /* weir u/s flag*/
#define BCT_WRSD 28       /* Weir D/S flag*/
#define BCT_WEIRD 29
#define BCT_WEIRU 30
#define BCT_SDR_NEU 31  /* NB SDR Condition Stage-Discharge GS*/
#define BCT_FLAPD 32   /* Flap gate Strings 32 - 35 */
#define BCT_FLAPU 33
#define BCT_FLPU 34
#define BCT_FLPD 35
#define BCT_CEQ 36		/* Equilibrium Concentration boundary condition */
#define POINT_SOURCE 37
#define BCT_NEU_LOAD 38 /* transport load */
#define BCT_EQT 39 /* Equilibrium temperature boundary GSAVANT */
#define BCT_DB_VEL 40 /* cjt :: for new EGS created during DB OVL */
#define BCT_SLSU 41   /* Sluice gate strings 41 - 44 */
#define BCT_SLSD 42
#define BCT_SLUICEU 43
#define BCT_SLUICED 44


/* equation types */
#define CONSERVATIVE 0		/* equations are conservative */
#define NONCONSERVATIVE 1	/* equations are nonconservative */

/* string types */
#define STR_NODE 1		/* a node string */
#define STR_EDGE 2		/* a edge string */
#define STR_FACE 3		/* a face string */
#define STR_MID  4		/* a mid string */
#define STR_E2F  5		/* an edge string destined to be a face string */
#define STR_ICE  6		/* an ice string */

/* initialization definitions */
#define UNSET_INT -3		/* an integer variable that has not been set */
#define UNSET_FLT -9999999.9 //3.0		/* a float variable that has not been set */

/* refinement tolerances - these are checked against element errors 
   the element errors are already scaled by user defined tolerances */
#define REF_TOL 1.0		/* refinement tolerance after scaling */
#define UNREF_TOL 0.1		/* percent of refine tolerance used to define unrefine tolerance */

/* allocation block for the sparse vectors */
#define SPV_BLOCK 75 // was 15
#define BV_BLOCK 15

/* communication stuff */
#ifdef _MESSG
#define MESSG_REQ_INC 20	/* the increment to allocate message requests */
#define MESSG_INT 2		/* integer message */
#define MESSG_DOUBLE 3		/* double message */
#define MESSG_PACKED 4		/* a packed message */
#define NMPI_FLAG 10		/* total number of asynchronous communication flags per pe */
#define TAG_UPDATE 990		/* message tag for update communications */
#define TAG_NODE_OUT 991	/* message tag for nodes out in repartitioning routines */
#define TAG_NODE_NUM 992	/* message tag for new node numbers in repartitioning routines */
#define TAG_NODE_DATA 993	/* message tab for node data in repartitioning routines */
#endif

/* type of element for element level communication */
#define COMM_3D_ELEM_LEVELS 3
#define COMM_2D_ELEM_LEVELS 2
#define COMM_1D_ELEM_LEVELS 1

/* post processing flags */
#define PS_FLUX 1		/* compute the fluxes through the face string */

/* linked list definitions - these are tied very closely to the 
   initialization of the object size array in tl_list_alloc */
#define EDGE_LIST 0		/* the list of edges */
#define FACE_LIST 1		/* the list of faces with the adjacent elements */
#define ELEM1D_LIST 2		/* the list of 1d elements */
#define ELEM2D_LIST 3		/* the list of 2d elements */
#define ELEM3D_LIST 4		/* the list of 3d elements */
#ifdef _MESSG
#define NODE_LIST 5		/* the list of global node numbers and their local counterparts */
#define ELEM_REF_LIST 6		/* the list of refined elements being sent to other processors */
#define NLIST 7			/* the number of linked lists */
#else
#define NLIST 5
#endif

/* the types of series */
#define ANY_SERIES 0		    /* flag to not check series type */
#define TIME_SERIES 1		    /* the series represents a time series */
#define OUTPUT_SERIES 2		    /* the series represents an output time series */
#define CONSTITUITIVE_SERIES 3  /* the series represents a constituitive equation */
#define WIND_SERIES 4           /* the series represents a wind series */
#define DT_SERIES 5             /* time-step series */
#define WAVE_SERIES 6           /* wave stress series */

/* time units */
#define SECONDS	0
#define MINUTES	1
#define HOURS	2
#define DAYS	3
#define WEEKS	4
#define MONTHS	5
#define YEARS	6

/* time conversion direction */
#define TO 0
#define FROM 1

/* ice friction coefficient flags */
#define BED_CO 1		/* bed coefficient */
#define ICE_CO 2		/* ice coefficient */
#define TOT_CO 3		/* total (bed+ice) coefficient */

/* model to use to calculate combined wave-current shear stress */
#define GM79   1
#define F84    2
#define HT91   3
#define DSK88  4
#define DATA13 5
#define DATA2  6

/* shifts for multiple equation indexing */
#define EQN4_ID11  0
#define EQN4_ID12  1
#define EQN4_ID13  2
#define EQN4_ID14  3
#define EQN4_ID21  4
#define EQN4_ID22  5
#define EQN4_ID23  6
#define EQN4_ID24  7
#define EQN4_ID31  8
#define EQN4_ID32  9
#define EQN4_ID33 10
#define EQN4_ID34 11
#define EQN4_ID41 12
#define EQN4_ID42 13
#define EQN4_ID43 14
#define EQN4_ID44 15

#define EQN3_ID11  0
#define EQN3_ID12  1
#define EQN3_ID13  2
#define EQN3_ID21  3
#define EQN3_ID22  4
#define EQN3_ID23  5
#define EQN3_ID31  6
#define EQN3_ID32  7
#define EQN3_ID33  8

#define EQN2_ID11  0
#define EQN2_ID12  1
#define EQN2_ID21  2
#define EQN2_ID22  3

#ifdef _ADH_ICM
/* ICM variable information */
#define NUM_ICM_DRQ 10
#define ICM_FACT 1000.0
typedef enum
  {
    /* Do not re-order these, insert new vars at end */ 
    UNSET_ICM_VAR = -1,
    ICM_VAR_SLN,	/* Salinity */
    ICM_VAR_TEM,	/* Temperature */
    ICM_VAR_SSI,	/* Suspended Solids */
    ICM_VAR_AG1,	/* Algal Group 1 */
    ICM_VAR_AG2,	/* Algal Group 2 */
    ICM_VAR_AG3,	/* Algal Group 3 */
    ICM_VAR_LDC,	/* Labile Dissolved Organic Carbon */
    ICM_VAR_RDC,	/* Refractory Dissolved Organic Carbon */
    ICM_VAR_LPC,	/* Labile Particulate Organic Carbon */
    ICM_VAR_RPC,	/* Refractory Particulate Organic Carbon */
    ICM_VAR_NH4,	/* Ammonium */
    ICM_VAR_NO3,	/* Nitrate */
    ICM_VAR_LDN,	/* Labile Dissolved Organic Nitrogen */
    ICM_VAR_RDN,	/* Refractory Dissolved Organic Nitrogen */
    ICM_VAR_LPN,	/* Labile Particulate Organic Nitrogen */
    ICM_VAR_RPN,	/* Refractory Particulate Organic Nitrogen */
    ICM_VAR_PO4,	/* Phosphate */
    ICM_VAR_LDP,	/* Labile Dissolved Organic Phosphorus */
    ICM_VAR_RDP,	/* Refractory Dissolved Organic Phosphorus */
    ICM_VAR_LPP,	/* Labile Particulate Organic Phosphorus */
    ICM_VAR_RPP,	/* Refractory Particulate Organic Phosphorus */
    ICM_VAR_COD,	/* Chemical Oxygen Demand */
    ICM_VAR_DOX,	/* Dissolved Oxygen */
    ICM_VAR_PTH,	/* Pathogen */
    ICM_VAR_TOX,	/* Contaminant */
    /* Add new vars above here and update the companion name arrays in 
       init_read_icm_xy and print_xms_pointers */
    /* This must be last */
    NUM_ICM_VARS  /* number of ICM variable types */
  } ICM_VAR;
#endif


/*Individual Debug Flags */
/* For example "PRINT_MATRIX" is currently level 16 of DEBUG_FE, but if a change were required */
/* it can be made in this one location rather than several others */

#define FE_MATRIX 	15	/* fe level 15 */
#define FE_RESID 	14	/* fe level 14 */
#define GRID_ENSPRINT 	15	/* grid level 15, prints ensight mesh exclusively */
#define GRID_GMSPRINT	14	/* grid level 14, prints gms mesh exclusively */

/* (cjt) */
#define NWVAR           13      /* number of wave variables on a node */

/*LP XDMF print control flags */
#define NODE_CENTERED 0
#define ELEM_3D_CENTERED 1
#define ELEM_2D_CENTERED 2

#define SCALAR_DATA 0
#define VECTOR2D_DATA 1
#define VECTOR3D_DATA 2

//Mark adding max Nodal NVAR allowed, used to create sparsity
#define MAX_NVAR 5
#define MAX_NNODE 6
#define MAX_ELEM_DOF 30 
