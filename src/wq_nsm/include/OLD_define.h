/* this is the standard definition file for the adh model */

/* standard header files */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <stdlib.h>
#include <stddef.h>
#include <limits.h>

#ifdef _MPI
#include <mpi.h>
#endif

/* constants */
#define EARTH_ROTAT 7.2722E-5	/* Earth's rotation in 1/s (coriolis: added 6-02) */
#define PI 3.141592653589793238462643383279502884197169399375105820974944592308
#define SQUAREROOTTWO 1.4142135623730950488016887242096980785696718753769480732
#define ONETHIRD 0.333333333333333333333333333333333333333333333333333333333333
#define ONESIXTH 0.166666666666666666666666666666666666666666666666666666666666
#define ONE_24TH 0.041666666666666666666666666666666666666666666666666666666666
#define ONE_72ND 0.013888888888888888888888888888888888888888888888888888888888
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

/* SEDIMENT ALGORITHM CONSTANTS FOR READING PARAMETERS FROM BC FILE */
#define COHSET 4
#define WINDWAVE 10
#define NONCOENT 1 /* (cjt) -- mdw SW3_FLOW svn merge */
#define DIVCOEF 3

/* input line definitions */
#define MAXLINE 150		/* the maximum number of characters on an input line */
#define MAXCARD 3		/* the maximum number of characters in a card */

/* geometric definitions */
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
#define SOLV_TOL FLT_EPSILON	/* the linear solver tolerance - usually single precision zero - 
				   NOTE:  this number MUST be significantly larger than SMALL */
#define TET_PSI_AT_CENTROID 0.25	/* the linear basis functions on a tet evaluated at the centroid */
#define SOLV_UPDATE_FACTOR 0.01	/* the factor for the update in the solver */
#define HASHSIZE 25013
#define CHAR_SIZE 12

/* definitions of names for various solves */
#define HYD 1
#define TRN 2
#define BLT 3
#define SLT 4
#define GRD 5
#define MVG 6
#define VRT 7

/* post processing flags - positive flags indicate transport quantities */
#define PS_FLAG_HEAD -1		/* flag for the heads */
#define PS_FLAG_SAT -2		/* flag for saturations */
#define PS_FLAG_OLHEAD -3	/* flag for overland flow heads */
#define PS_FLAG_CHHEAD -4	/* flag for channel heads */
#define PS_FLAG_PRS -5		/* flag for pressures */
#define PS_FLAG_VEL -6		/* flag for velocities */
#define PS_FLAG_DPL -7		/* flag for displacements */
#define PS_FLAG_OLVEL -8	/* flag for overland velocities */
#define PS_FLAG_GSP -9		/* flag for grid speed */
#define PS_FLAG_ERR -10		/* flag for error indicator */
#define PS_FLAG_FLX -11		/* flag for nodal flux */
#define PS_FLAG_HEAT -12	/* flag for the temperatures */
#define PS_FLAG_UGD -13		/* flag for x-vel gradient */
#define PS_FLAG_VGD -14		/* flag for y-vel gradient */
#define PS_FLAG_WGD -15		/* flag for z-vel gradient */
#define PS_FLAG_WS_POT -16	/* flag for potential is SW3 */
#define PS_FLAG_VEL_COR -17	/* flag for velocity correction in SW3 */
#define PS_FLAG_POLHEAD -18	/* flag for old overland flow heads */
#define PS_FLAG_POLVEL -19	/* flag for old overland velocities */
#define PS_FLAG_ALT -20		/* flag for active layer thickness */
#define PS_FLAG_ALD -21		/* flag for active layer distribution */
#define PS_FLAG_BED -22		/* flag for bedload vector */
#define PS_FLAG_SHSTR -23	/* flag for shear stresses */
#define PS_FLAG_2D_SC -24	/* flag for 2D scatter plot */
#define PS_FLAG_3D_SC -25	/* flag for 3D scatter plot */
#define PS_FLAG_DEPTH -26	/* flag for depth */
#define PS_FLAG_PDEPTH -27	/* flag for old depth */
#define PS_FLAG_BSH -28		/* flag for bed shear stress */
#define PS_FLAG_VORT -29	/* flag for SW2 vorticity */
#define PS_FLAG_SFLX -30	/* flag for surface heat flux */
#define PS_FLAG_PRS_HOT -31	/* flag for pressure hotstart creation (NS) */
#define PS_FLAG_VEL_HOT -32	/* flag for velocity hotstart creation (NS) */
#define PS_FLAG_DEP_HOT -33	/* flag for depth hotstart creation (SW2) */
#define PS_FLAG_PDP_HOT -34	/* flag for previous depth hotstart creation (SW2) */
#define PS_FLAG_OVL_HOT -35	/* flag for velocity hotstart creation (SW2) */
#define PS_FLAG_POV_HOT -36	/* flag for previous velocity hotstart creation (SW2) */
#define PS_FLAG_DIS_HOT -37	/* flag for displacement hotstart creation (SW2!) */
#define PS_FLAG_CON_HOT -38	/* flag for concentration hotstart creation (SW2) */
#define PS_FLAG_PCN_HOT -39	/* flag for previous concentration hotstart creation (SW2) */
#define PS_FLAG_DSP_HOT -40	/* flag for displacement hotstart creation (NS!) */
#define PS_FLAG_GSP_HOT -41	/* flag for grid speed hotstart creation (NS) */
#define PS_FLAG_HED_HOT -42	/* flag for ....head hotstart creation (GW) */
#define PS_FLAG_TMP_HOT -43	/* flag for temperature hotstart creation (HT) */
#define PS_FLAG_SMR -44		/* flag for SW2 sed mass balance *//* these must be in consecutive order!! - dss */
#define PS_FLAG_SLN_HOT -45	/* flag for salinity hotstart creation (ICM) */
#define PS_FLAG_TEM_HOT -46	/* flag for temperature hotstart creation (ICM) */
#define PS_FLAG_SSI_HOT -47	/* flag for suspended solids hotstart creation (ICM) */
#define PS_FLAG_AG1_HOT -48	/* flag for algal group 1 hotstart creation (ICM) */
#define PS_FLAG_AG2_HOT -49	/* flag for algal group 2 hotstart creation (ICM) */
#define PS_FLAG_AG3_HOT -50	/* flag for algal group 3 velocity hotstart creation (ICM) */
#define PS_FLAG_LDC_HOT -51	/* flag for labile dissolved organic carbon hotstart creation (ICM) */
#define PS_FLAG_RDC_HOT -52	/* flag for refractory dissolved organic carbon hotstart creation (ICM) */
#define PS_FLAG_LPC_HOT -53	/* flag for labile particulate organic carbon hotstart creation (ICM) */
#define PS_FLAG_RPC_HOT -54	/* flag for refractory particulate organic carbon hotstart creation (ICM) */
#define PS_FLAG_NH4_HOT -55	/* flag for ammonium hotstart creation (ICM) */
#define PS_FLAG_NO3_HOT -56	/* flag for nitrate hotstart creation (ICM) */
#define PS_FLAG_LDN_HOT -57	/* flag for labile dissolved organic nitrogen hotstart creation (ICM) */
#define PS_FLAG_RDN_HOT -58	/* flag for refractory dissolved organic nitrogen hotstart creation (ICM) */
#define PS_FLAG_LPN_HOT -59	/* flag for labile particulate organic nitrogen hotstart creation (ICM) */
#define PS_FLAG_RPN_HOT -60	/* flag for refractory particulate organic nitrogen hotstart creation (ICM) */
#define PS_FLAG_PO4_HOT -61	/* flag for phosphate hotstart creation (ICM) */
#define PS_FLAG_LDP_HOT -62	/* flag for labile dissolved organic phosphorus hotstart creation (ICM) */
#define PS_FLAG_RDP_HOT -63	/* flag for previous refractory dissolved organic phosphorus hotstart creation (ICM) */
#define PS_FLAG_LPP_HOT -64	/* flag for labile particulate organic phosphorus hotstart creation (ICM) */
#define PS_FLAG_RPP_HOT -65	/* flag for refractory particulate organic phosphorus hotstart creation (ICM) */
#define PS_FLAG_COD_HOT -66	/* flag for chemical oxygen demand hotstart creation (ICM) */
#define PS_FLAG_DOX_HOT -67	/* flag for dissolved oxygen hotstart creation (ICM) */
#define PS_FLAG_PTH_HOT -68	/* flag for pathogen hotstart creation (ICM) */
#define PS_FLAG_TOX_HOT -69	/* flag for contaminant hotstart creation (ICM) */

#define PS_FLAG_SWAVE_HEIGHT -70 /* cjt flag for wind-waves output (WAVE) */
#define PS_FLAG_SWAVE_FORCES -71 /* cjt flag for wind-wave forces (WAVE) */
/* ****************************************** */
#define PS_FLAG_TAG_HOT -70	/* flag for nodal material hotstart creation (ICM) */


#define PS_FLAG_BED_ELV -71 /* flag for bed elev */
#define PS_FLAG_ERR_HYDRO -72	/* flag for error indicator-hydro part */
#define PS_FLAG_SUS -73     /* flag for suspended load vector */
#define PS_FLAG_PDPL -74 /* flag for previous time step displacement */
#define PS_FLAG_PVEL -75 /* flag for previous time step velocity */

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
#define K   13
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
#define UNREFINE -1             /* unrefine the element */
#define NOUNR 0                 /* refined element */
#define CONSV 1                 /* merged element at t(n-1) */
#define UNREF 2                 /* merged element at t(n) */ 

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
#define BCT_CEQ 36		/* Equilibrium Concentration boundary condition */

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
#define BCT_FLAPD 32
#define BCT_FLAPU 33
#define BCT_FLPU 34
#define BCT_FLPD 35
#define BCT_CEQ 36		/* Equilibrium Concentration boundary condition */
#define POINT_SOURCE 37
#define BCT_NEU_LOAD 38 /* transport load */

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
#define UNSET_FLT -3.0		/* a float variable that has not been set */

/* refinement tolerances - these are checked against element errors 
   the element errors are already scaled by user defined tolerances */
#define REF_TOL 1.0		/* refinement tolerance after scaling */
#define UNREF_TOL 0.1		/* percent of refine tolerance used to define unrefine tolerance */

/* allocation block for the sparse vectors */
#define SPV_BLOCK 15
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
#define ANY_SERIES 0		/* flag to not check series type */
#define TIME_SERIES 1		/* the series represents a time series */
#define OUTPUT_SERIES 2		/* the series represents an output time series */
#define CONSTITUITIVE_SERIES 3	/* the series represents a constituitive equation */
#define WIND_SERIES 4		/* the series represents a wind series */

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

/* ICM variable information */
#define NICMDRQ 10
#define NICMVARS  25
#define ICMFACT 1000.0
#define SLN  0			/* Salinity */
#define TEM  1			/* Temperature */
#define SSI  2			/* Suspended Solids */
#define AG1  3			/* Algal Group 1 */
#define AG2  4			/* Algal Group 2 */
#define AG3  5			/* Algal Group 3 */
#define LDC  6			/* Labile Dissolved Organic Carbon */
#define RDC  7			/* Refractory Dissolved Organic Carbon */
#define LPC  8			/* Labile Particulate Organic Carbon */
#define RPC  9			/* Refractory Particulate Organic Carbon */
#define NH4 10			/* Ammonium */
#define NO3 11			/* Nitrate */
#define LDN 12			/* Labile Dissolved Organic Nitrogen */
#define RDN 13			/* Refractory Dissolved Organic Nitrogen */
#define LPN 14			/* Labile Particulate Organic Nitrogen */
#define RPN 15			/* Refractory Particulate Organic Nitrogen */
#define PO4 16			/* Phosphate */
#define LDP 17			/* Labile Dissolved Organic Phosphorus */
#define RDP 18			/* Refractory Dissolved Organic Phosphorus */
#define LPP 19			/* Labile Particulate Organic Phosphorus */
#define RPP 20			/* Refractory Particulate Organic Phosphorus */
#define COD 21			/* Chemical Oxygen Demand */
#define DOX 22			/* Dissolved Oxygen */
#define PTH 23			/* Pathogen */
#define TOX 24			/* Contaminant */

/* defines the card values for the input */
/* these numbers are calculate using base 37 where digits are 
   defined by (' ',0),(A-Z,1-26),(0-9,27-36)
   the first digit is i1, second i2, and third i3
   the calculated number then is found by
   1369*i1 + 37*i2 + i3;
   for example CARD_TRN
   i1 = T = 20
   i2 = R = 18
   i3 = N = 14
   1369 * 20 + 37 * 18 + 14 = 28060
 */
#define CARD_ADP  1533
#define CARD_AG1  1656
#define CARD_AG2  1657
#define CARD_AG3  1658
#define CARD_ALB  1815
#define CARD_APH  1969
#define CARD_ATF  2115
#define CARD_ATT  2129	/* wind attenuation card (cjt) */
#define CARD_AWB  2222
#define CARD_BCH  2857 /* Breach CARD */
#define CARD_BD   2886
#define CARD_BED  2927
#define CARD_BLK  3193
#define CARD_BNS  3275
#define CARD_BRH  3412
#define CARD_BT	  3478
#define CARD_BTS  3497
#define CARD_CBA  4182
#define CARD_CBM  4194
#define CARD_CBN  4195
#define CARD_CF   4329
#define CARD_CFW  4352
#define CARD_CH   4403
#define CARD_CIC  4443
#define CARD_CLA  4552
#define CARD_CN   4625
#define CARD_COD  4666
#define CARD_CON  4676
#define CARD_COR  4680
#define CARD_CPA  4700
#define CARD_CPM  4712
#define CARD_CPN  4713
#define CARD_CSV  4832
#define CARD_C0   5106
#define CARD_DB   5550
#define CARD_DBG  5557
#define CARD_DEN  5675
#define CARD_DF   5698
#define CARD_DIS  5828
#define CARD_DOX  6055
#define CARD_DPL  6080
#define CARD_DPT  6088
#define CARD_DRG  6149
#define CARD_DSF  6185
#define CARD_DSS  6198
#define CARD_DS0  6206
#define CARD_DTL  6228
#define CARD_EEV  7052
#define CARD_EFL  7079
#define CARD_EFS  7086
#define CARD_EFT  7087
#define CARD_EGS  7123
#define CARD_ELM  7302   /* Card for ELAM format print */
#define CARD_ELS  7308   /* Card for ELAM nodes */
#define CARD_EMS  7345
#define CARD_END  7367
#define CARD_EQ   7474
#define CARD_ERH  7519
#define CARD_EV   7659
#define CARD_EVS  7678
#define CARD_EXT  7753
#define CARD_E3T  7975
#define CARD_E4T  8012
#define CARD_FCS  8344
#define CARD_FET  8419
#define CARD_FGT  8493
#define CARD_FIN  8561  /* External file reading (cjt) */
#define CARD_FLI  8667
#define CARD_FLP  8674
#define CARD_FLT  8678
#define CARD_FLW  8681
#define CARD_FLX  8682
#define CARD_FNI  8741
#define CARD_FNS  8751
#define CARD_FPS  8825
#define CARD_FR   8880
#define CARD_FRC  8883
#define CARD_FRS  8899
#define CARD_FRT  8900
#define CARD_FUT  9011
#define CARD_FVS  9047
#define CARD_G    9583
#define CARD_GNS 10120
#define CARD_GW  10434
#define CARD_HFX 11198
#define CARD_HID 11289
#define CARD_HOT 11527
#define CARD_HS  11655
#define CARD_HSP 11671
#define CARD_HT  11692
#define CARD_HTS 11716
#define CARD_IAC 12361
#define CARD_ICE 12437
#define CARD_ICM 12445
#define CARD_IDT 12489
#define CARD_IMN 12816
#define CARD_INC 12842
#define CARD_INJ 12849
#define CARD_INS 12858
#define CARD_INT 12859
#define CARD_IP  12913
#define CARD_IR  12987
#define CARD_IRH 12995
#define CARD_ITL 13073
#define CARD_JUL 14479
#define CARD_K   15059
#define CARD_KR  15725
#define CARD_LA  16465
#define CARD_LAT 16485
#define CARD_LDC 16579
#define CARD_LDE 16581
#define CARD_LDH 16584
#define CARD_LDN 16590
#define CARD_LDP 16592
#define CARD_LID 16765
#define CARD_LPC 17023
#define CARD_LPN 17034
#define CARD_LPP 17036
#define CARD_LVL 17254  /* Output Level */
#define CARD_L0  17427
#define CARD_MDS 17964
#define CARD_MEO 17997
#define CARD_MET 18002
#define CARD_MIT 18150
#define CARD_MKS 18223
#define CARD_ML  18241
#define CARD_MNG 18322
#define CARD_MNS 18334
#define CARD_MP  18389
#define CARD_MPA 18390
#define CARD_MPM 18402
#define CARD_MTS 18556
#define CARD_MU  18574
#define CARD_MUC 18577
#define CARD_MU0 18601
#define CARD_MV  18611
#define CARD_MVI 18620
#define CARD_MVS 18630
#define CARD_NB  19240
#define CARD_NBE 19245
#define CARD_NBL 19252
#define CARD_NC  19277
#define CARD_NCE 19282
#define CARD_NCP 19293
#define CARD_ND  19314
#define CARD_NDA 19315
#define CARD_NDM 19327
#define CARD_NDP 19330
#define CARD_NDS 19333
#define CARD_NFS 19407
#define CARD_NF2 19417
#define CARD_NH4 19493
#define CARD_NIT 19519
#define CARD_NO3 19751
#define CARD_NRT 19852
#define CARD_NS  19869
#define CARD_NSE 19874
#define CARD_NSF 19875
#define CARD_NSM 19882
#define CARD_NTL 19918
#define CARD_NUT 19963
#define CARD_NVA 19981
#define CARD_NVM 19993
#define CARD_OB  20609
#define CARD_OC  20646
#define CARD_OF  20757
#define CARD_OFF 20763
#define CARD_OHD 20835
#define CARD_OIC 20871
#define CARD_OL  20979
#define CARD_OLD 20983
#define CARD_OP  21127
#define CARD_OS  21238
#define CARD_OTW 21298
#define CARD_OUT 21332
#define CARD_OVH 21357
#define CARD_OVL 21361
#define CARD_PC  22015
#define CARD_PHS 22219
#define CARD_POR 22477
#define CARD_PO4 22490
#define CARD_PRA 22571
#define CARD_PRE 22575
#define CARD_PRS 22589
#define CARD_PTH 22652
#define CARD_RAY 24704
#define CARD_RCT 24773
#define CARD_RD  24790
#define CARD_RDA 24791
#define CARD_RDC 24793
#define CARD_RDN 24804
#define CARD_RDP 24806
#define CARD_RHO 24953
#define CARD_RKE 25054
#define CARD_RPC 25237
#define CARD_RPN 25248
#define CARD_RPP 25250
#define CARD_RTL 25394
#define CARD_R0  25641
#define CARD_R00 25668
#define CARD_R01 25669
#define CARD_R02 25670
#define CARD_R03 25671
#define CARD_R04 25672
#define CARD_R05 25673
#define CARD_R06 25674
#define CARD_R07 25675
#define CARD_R08 25676
#define CARD_R09 25677
#define CARD_R10 25705
#define CARD_R11 25706
#define CARD_R12 25707
#define CARD_R13 25708
#define CARD_R14 25709
#define CARD_R15 25710
#define CARD_R16 25711
#define CARD_R17 25712
#define CARD_R18 25713
#define CARD_R19 25714
#define CARD_R20 25742
#define CARD_R21 25743
#define CARD_R22 25744
#define CARD_R23 25745
#define CARD_R24 25746
#define CARD_R25 25747
#define CARD_R26 25748
#define CARD_R27 25749
#define CARD_R28 25750
#define CARD_R29 25751
#define CARD_SAL 26060
#define CARD_SAT 26068
#define CARD_SAV 26070
#define CARD_SBA 26086
#define CARD_SBM 26098
#define CARD_SBN 26099
#define CARD_SDI 26168
#define CARD_SDR 26177
#define CARD_SDV 26181
#define CARD_SED 26200
#define CARD_SFA 26234
#define CARD_SFC 26236
#define CARD_SFM 26246
#define CARD_SFN 26247
#define CARD_SGG 26277
#define CARD_SGS 26289
#define CARD_SGW 26293
#define CARD_SHG 26314
#define CARD_SHS 26326
#define CARD_SHW 26330
#define CARD_SI  26344
#define CARD_SLT 26475
#define CARD_SMB 26494
#define CARD_SND 26533
#define CARD_SOC 26569
#define CARD_SP  26603
#define CARD_SPL 26615
#define CARD_SRC 26680
#define CARD_SRT 26697
#define CARD_SS  26714
#define CARD_SSI 26723
#define CARD_SST 26734
#define CARD_STD 26755
#define CARD_STR 26769 /* wind stress transform type (cjt) */
#define CARD_STH 26759
#define CARD_STL 26763
#define CARD_SW  26862
#define CARD_SW2 26891
#define CARD_SW3 26892
#define CARD_TC  27491
#define CARD_TCN 27505
#define CARD_TDS 27547
#define CARD_TEM 27578
#define CARD_TF  27602
#define CARD_THD 27680
#define CARD_THK 27687
#define CARD_THT 27696
#define CARD_TID 27717  /* tidal boundary condition card (cjt) */
#define CARD_TKE 27792
#define CARD_TKG 27794
#define CARD_TKM 27800
#define CARD_TKW 27810
#define CARD_TMN 27875
#define CARD_TMP 27877
#define CARD_TOR 27953
#define CARD_TOX 27959
#define CARD_TPG 27979
#define CARD_TRI 28055
#define CARD_TRL 28058
#define CARD_TRN 28060
#define CARD_TUR 28175
#define CARD_TRT 28066
#define CARD_TSC 28086
#define CARD_TUT 28177
#define CARD_TVS 28213
#define CARD_T0  28379
#define CARD_URV 29437
#define CARD_USR 29470
#define CARD_U0  29748
#define CARD_VEL 30315
#define CARD_VOR 30691
/* #define CARD_WAT 31544   no longer needed (cjt) */
#define CARD_WAV 31546  /* Material card for wave stress calculations (cjt) */
#define CARD_WDR 31653
#define CARD_WDT 31655
#define CARD_WER 31690
#define CARD_WL  31931
#define CARD_WMD 31972
#define CARD_WND 32009  /* Material card for wind stress calculations (cjt) */
#define CARD_WRS 32172
/* #define CARD_WV  32301  no longer used (cjt) */
#define CARD_WWS 32357
#define CARD_XYC 33784
#define CARD_XYT 33801
#define CARD_XY1 33809
#define CARD_XY2 33810

/*Individual Debug Flags */
/* For example "PRINT_MATRIX" is currently level 16 of DEBUG_FE, but if a change were required */
/* it can be made in this one location rather than several others */

#define FE_MATRIX 	15	/* fe level 15 */
#define FE_RESID 	14	/* fe level 14 */
#define GRID_ENSPRINT 	15	/* grid level 15, prints ensight mesh exclusively */
#define GRID_GMSPRINT	14	/* grid level 14, prints gms mesh exclusively */

/* (cjt) */
#define NWVAR           13      /* number of wave variables on a node */
