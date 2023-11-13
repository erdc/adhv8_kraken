#ifndef __CARDS_H__
#define __CARDS_H__

/* Define input card values */
/* Calculate integer value equivalence of file card using base 37, where individual
   characters A-Z and a-z are converted to 1-26 and 0-9 are converted to 27-36 and
   combined:

   card_value = card[0] + (37 * card[1]) + (37^2 * card[2]) + ...

   For example: "XY1" converts to:

          'X'           'Y'           '1'
   (37^0 * 24) + (37^1 * 25) + (37^2 * 28) = 24 + (37 * 25) + (37^2 * 28) = 39281
    
   This formulation allows for any length card.
 */

#define MAX_CARD_LENGTH 6 /* this defines the integer size used in parse_card */
typedef enum {

    UNSET_CARD = -1,
    NO_CARD = 0,

    CARD_NC     = 125,
    CARD_TS     = 723,

    // Superfile
    CARD_BC     = 113,
    CARD_GEO    = 20727,    // 2D GEO CARD
    CARD_GGEO   = 766906,   // 3D GEO CARD
    CARD_HOT    = 27943,
    CARD_FACE   = 257415,
    CARD_NMOD   = 223642,
    CARD_MODEL  = 22749241,
    CARD_LINK   = 576694,
    CARD_FLUX   = 1244871,
    CARD_FULL   = 625047,
    CARD_SMODEL = 841721936,

    // geo cards
    CARD_ND     = 162,
    CARD_E3T    = 28495,
    CARD_E4T    = 28532,
    CARD_E6P    = 21539,
    CARD_NFS    = 26247,
    CARD_FCS    = 26128,
    CARD_GRID   = 215606,
    CARD_TRI    = 13007,
    CARD_QUAD   = 204775,
    CARD_TET    = 27585,
    CARD_PRISM  = 25339503,
    CARD_SCALEX = 1674235108,
    CARD_SCALEY = 1743579065,
    CARD_SCALEZ = 1812923022,
    CARD_NOCOL  = 23254403,
    CARD_COL    = 16986,

    // hotstart cards
    CARD_NAME   = 271113,
    CARD_ENDDS  = 35817670,
    CARD_BEGSCL = 838722144,
    CARD_BEGVEC = 218526812,

    // bc cards - series
    CARD_SERIES = 1327386711,
    CARD_WRITE  = 10396875,
    CARD_AWRITE = 384684376,
    CARD_DT     = 744,
    CARD_BCS    = 26124,
    CARD_WIND   = 222134,
    CARD_WAVE   = 283443,
    CARD_SOURCE = 353283345,
    CARD_CONSTI = 662560964, //constitutive relations
    // bc cards - strings
    CARD_NDS    = 26173,
    CARD_EGS    = 26275,
    CARD_MDS    = 26172,    
    CARD_MTS    = 26764,
    
    // bc cards 
    CARD_TC     = 131,
    CARD_OP     = 607,
    CARD_IP     = 601,
    CARD_PC     = 127,
    CARD_OC     = 126,
    CARD_OFF    = 8451,
    CARD_MP     = 605,
    CARD_DB     = 78,
    CARD_NB     = 88,
    CARD_FR     = 672,
    CARD_END    = 5999, 
    CARD_MG     = 272,

    // bc cards - Operational parameters
    CARD_TRN    = 19852,
    CARD_INC    = 4634,
    CARD_BLK    = 15505,
    CARD_NS     = 717,
    CARD_NS2    = 40418,
    CARD_NS3    = 41787,
    CARD_SW     = 870,
    CARD_SW2    = 40571,
    CARD_SW3    = 41940,
    CARD_WAV    = 30178,   /* Material card for wave stress calculations */ 
    CARD_WND    = 6017,    /* Material card for wind stress calculations */
    CARD_TEM    = 18002,
    CARD_PRE    = 7527,
    CARD_TPG    = 10195,
    CARD_DIF    = 8551, /* Card designating diffusive wave equation set */
    CARD_EOS    = 26571, /* Card for equation of state */
    CARD_GW     = 858, /* Card for groundwater physics */

    // bc cards - Iteration parameters
    CARD_MIT    = 27726,
    CARD_NIT    = 27727,
    CARD_ITL    = 17177,
    CARD_NTL    = 17182,

    // bc cards - time control
    CARD_T0     = 1019,
    CARD_TF     = 242,
    CARD_IDT    = 27537,
    CARD_NDP    = 22066,

    // bc cards - print control
    CARD_LVL    = 17254,   /* Output Level */

    // bc cards - materials
    CARD_ML     = 457,
    CARD_COR    = 25200,
    CARD_EVS    = 26830,
    CARD_EEV    = 30308,
    CARD_SRT    = 28065,
    CARD_MU     = 790,
    CARD_MUC    = 4897,
    CARD_G      = 7,
    CARD_RHO    = 20849,
    CARD_DTL    = 17172,
    CARD_ATT    = 28121,   /* wind attenuation card */
    CARD_STR    = 25401,   /* wind stress transform type */
    CARD_TUR    = 25439,   /* Turbulence card */
    CARD_HYDCON = 999086181,
    CARD_SUCK   = 562086,
    CARD_SATDEP = 1119104165,
    // groundwater specific
    CARD_FRT    = 28052,
    CARD_SS     = 722,
    CARD_K      = 11,
    CARD_POR    = 25213,
    CARD_FVS    = 26831,
    CARD_RSD    = 6197,
    CARD_VGA    = 1650,
    CARD_VGN    = 19447,
    CARD_VGP    = 22185,
    CARD_VGX    = 33137,
    CARD_BCL    = 16541,
    CARD_BCP    = 22017,
    CARD_BCE    = 6958,
    CARD_BCX    = 32969,
    CARD_KR     = 677,
    CARD_SAT    = 27436,
    CARD_TOR    = 25217,
    CARD_DPL    = 17024,
    CARD_DPT    = 27976,
    CARD_SSA    = 2091,
    CARD_BUL    = 17207,
    // bc cards - dirichlet boundaries
    CARD_VEL    = 16635,
    CARD_PRS    = 26693,
    CARD_FRS    = 26683,
    CARD_OVL    = 17257,
    CARD_THD    = 5792,
    CARD_FLW    = 31937,
    CARD_PSI    = 13040,
    // bc cards - natural boundaries
    CARD_FLXNML = 857234025,
    CARD_DIS    = 26348,
    CARD_OTW    = 32242,
    CARD_OVH    = 11781,
    CARD_BED    = 5663,
    CARD_OF     = 237,


	// bc cards Structures
	CARD_WER    = 24850,
	CARD_WRS    = 26700,
	CARD_FLP    = 22354,
	CARD_FGT    = 27645,
	CARD_SLS    = 26474,
	CARD_SLUICE = 352827357,

    // bc cards - friction
    CARD_MNG    = 10114,
    CARD_ERH    = 11623,

    // bc cards - transport constituents
    CARD_DF     = 226,
    CARD_TRT    = 28066,
    CARD_CN     = 521,
    CARD_CON    = 19724,
    CARD_SAL    = 16484,
	CARD_TMP    = 22405,
	CARD_VOR    = 25219,

    // PDE TERMS CARDS (for leaving terms out)
    CARD_NOTERM = 935487553,
    CARD_DIFF   = 312469,
    CARD_CONV   = 1134090,
    CARD_SUPG   = 377271,
    CARD_TIME   = 271415,
    CARD_FRIC   = 164952,
    CARD_PRESS  = 36578993,
	CARD_HYDRO  = 29030578,

    // DEBUG PRINT CARDS
    CARD_DEBUG  = 14185767,
    CARD_NEWTON = 999972559,
    CARD_RES  = 26214,
    CARD_RESHV = 41662980,
    CARD_RESWV = 42422775,
    CARD_LOAD   = 204548,
    CARD_LOADHV = 1540764890,
    CARD_LOADWV = 1568877305,
    CARD_MAT = 27430,
    CARD_MATHV = 41664196,
    CARD_MATWV = 42423991,
    CARD_RHS = 26325,
    CARD_RHSHV = 41663091,
    CARD_RHSWV = 42422886,
    CARD_READBC = 211984377,

    // STOP THE CODE CARD
    CARD_STOP = 831742,

    // SCREEN OUTPUT OPTIONS
    CARD_SOUT = 1042383,
    CARD_RESID = 7978735,
    CARD_ALL = 16873,
    CARD_NLNODE = 354995848,
    CARD_LNODE = 9594482,
    CARD_MERROR = 1277240235,

    // TEST CASE OUTPUT
    CARD_TEST = 1039276,
    CARD_XCOR = 932424,
    CARD_YCOR = 932425,
    CARD_PTEST = 38453228,
    CARD_SLOSH = 15976693,
    CARD_WNDH = 411241,
    CARD_SALT = 1029544,
    CARD_LOCK = 561857,
    CARD_WD = 171,
    CARD_DAM = 17838,
	CARD_FRONT = 38213569,
    CARD_CSLOSH = 591137644,
    CARD_TIDE = 259094,
    CARD_DWEH = 412924,

////////////////////////////////////////////////////////////////////////////////////////////////////
// FOLLOWING LINES ADDED BY GAJANAN [ gkc June 2016 ]. These are for 2D-3D coupling.              //
    CARD_SLSH23 = 2135101077,
    CARD_TIDE23 = 2134928473,
// ABOVE LINES ADDED BY GAJANAN                                                                   //
////////////////////////////////////////////////////////////////////////////////////////////////////

    // SEDLIB
    CARD_SEDLIB = 156168879,
    CARD_GRAIN = 26696173,
    CARD_NBL = 16516,
    CARD_SBA = 1462,
    CARD_CBA = 1446,
    CARD_CLA = 1816,
    CARD_SND = 6013,
    CARD_NSF = 8931,
    CARD_NCP = 22029,
	CARD_SP     = 611,
	CARD_HID    = 5817,
	CARD_NSE    = 7562,
	CARD_NBE    = 6933,
	CARD_CSV    = 30824,
	CARD_SBM    = 17890,
    CARD_NDM    = 17959,

//////////////////////////////////////////////////////////////////////////////////////////////////
// FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. //
    CARD_WNDLIB = 156169216,
// ABOVE LINES ADDED BY GAJANAN                                                                 //
//////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
// FOLLOWING LINES ADDED BY GAJANAN AND COREY [ gkc cjt May 2016 ]. These are for 2D-3D coupling. //
        CARD_NUMIFC = 219751302,    // gkc April 2016 for coupling
        CARD_INTFCE = 352674093,    // gkc April 2016 for coupling
        CARD_HY = 933,    // cjt 2d-3d coupling
        CARD_INT = 27907, // cjt 2d-3d coupling
        CARD_CPL = 17023,           // cjt 2d-3d coupling
        CARD_EXT = 28273,           // gkc 2d-3d flux coupling
        CARD_SUPIFC = 219755414,    // gkc 2d-3d flux coupling
        CARD_IFCSM  = 25330838,     // gkc 2d-3d flux coupling
// ABOVE LINES ADDED BY GAJANAN AND COREY                                                         //
////////////////////////////////////////////////////////////////////////////////////////////////////

	// Water quality
	CARD_WQNSM = 25346318,
	CARD_NSM = 18514,
	CARD_NO2 = 40270,
	CARD_NO3 = 41639,
	CARD_NH4 = 42749,
	CARD_ON = 533,
	CARD_DP = 596,
	CARD_PO4 = 43010,
	CARD_ALG = 10028,
	CARD_CBOD = 223224,
	CARD_DO = 559,
	CARD_ALP0 = 1389980,
	CARD_ALP1 = 1440633,
	CARD_ALP2 = 1491286,
	CARD_ALP3 = 1541939,
	CARD_ALP4 = 1592592,
	CARD_ALP5 = 1643245,
	CARD_ALP6 = 1693898,
	CARD_BET1 = 1445851,
	CARD_BET2 = 1496504,
	CARD_BET3 = 1547157,
	CARD_BET4 = 1597810,
	CARD_K1 = 1047,
	CARD_K2 = 1084,
	CARD_K3 = 1121,
	CARD_K4 = 1158,
	CARD_GRO = 21208,
	CARD_SIG1 = 1428219,
	CARD_SIG2 = 1478872,
	CARD_SIG3 = 1529525,
	CARD_SIG4 = 1580178,
	CARD_SIG5 = 1630831,
	CARD_KL = 445,
	CARD_KN = 529,
	CARD_KP = 603,
	CARD_PN = 534,
	CARD_LAM0 = 1385477,
	CARD_LAM1 = 1436130,
	CARD_LAM2 = 1486783,
	CARD_MDO = 20696,
	CARD_MTMP = 828998,

    // OPTIONAL FILE OUTPUT
    CARD_FOUT = 1042370,
    CARD_SED = 5680,
    CARD_ADAPT = 38295186,
    CARD_BEDVEL = 842618318,
    CARD_SURVEL = 842638093,
    CARD_AVGVEL = 842623053,
    CARD_GSPEED = 287022512,
	CARD_VIS = 26366, /* GSAVANT */
	//CARD_DIF = 8551,  /* GSAVANT */
	CARD_CHOP = 831282,
    CARD_XDF = 8386, /*added by GAJANAN and MWF for XDMF*/

    // utility cards
    CARD_BIN = 19501

} CARD;

#endif
