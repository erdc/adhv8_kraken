#pragma once
#ifndef NSMV7_0_H
#define NSMV7_0_H



#ifdef _USRDLL
#define IMPORT __declspec(dllimport)
#else
#define IMPORT
#endif



#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif





typedef long boolean;


 
enum  NSM_SEDIMENT_SORBING_ENUM
{
    NSM_SEDIMENT_IS_SORBING_Enum,
    NSM_SEDIMENT_IS_NOT_SORBING_Enum
};




typedef struct 
{    
    unsigned long numberSediments;

	double  bindingEffectivenessCoefficient;	// De: DOC binding


    // Following are arrays of size: numberSediments

    int*  tag;                              // Array of ints for host use - NOT USED by NSM

    enum NSM_SEDIMENT_SORBING_ENUM*     sorbing;     // Specifies if sediment class is capable of sorbing

	double*  concentration;			        // Aquatic : M/L^3

 	double*  FOC;							// Fraction Organic Carbon

	double*  particleInteractionParameter;	// DiToro's Vx parameter 

    double*  P;                             // Spread fraction - internal

} NSMSedimentEnvironment;






// NSMSedimentEnvironment   P R O T O T Y P E S


EXTERN IMPORT  NSMSedimentEnvironment*  NSMSedimentEnvironment_Create(unsigned long numberSediments);


EXTERN IMPORT  void   NSMSedimentEnvironment_Init(
                    NSMSedimentEnvironment* p, 
                    unsigned long numberSediments);


EXTERN IMPORT  NSMSedimentEnvironment*  NSMSedimentEnvironment_Clone(NSMSedimentEnvironment* p);


EXTERN IMPORT  NSMSedimentEnvironment*  NSMSedimentEnvironment_Delete(NSMSedimentEnvironment* p);









// Library structures

 struct NSMListNode
{
	void* item;

	struct NSMListNode* next;
}; 

typedef struct NSMListNode NSMListNode;









// P A R T I C L E   D E S C R I P T I O N 

typedef struct 
{
	double  aquaticConcentration;			// Aquatic : M/L^3

	double  FOC;							// Particle's Fraction Organic Carbon

	double  particleInteractionParameter;	// DiToro's Vx parameter

} NSMCSchemaParticle;



// NSMCSchemaParticle    P R O T O T Y P E S


EXTERN IMPORT  NSMCSchemaParticle*  NSMCSchemaParticle_Create();


EXTERN IMPORT  void  NSMCSchemaParticle_Init(NSMCSchemaParticle*);


EXTERN IMPORT  NSMCSchemaParticle*  NSMCSchemaParticle_Clone(NSMCSchemaParticle*);


EXTERN IMPORT  NSMCSchemaParticle*  NSMCSchemaParticle_Delete(NSMCSchemaParticle*);










// P A R T I C L E - D I S T R I B U T I O N    D E S C R I P T I O N 

typedef struct
{
	unsigned long  numberOfParticles;

	double  bindingEffectivenessCoefficient;	// De: DOC binding

	NSMListNode* particleList;

}  NSMCSchemaParticleDistribution;




// NSMCSchemaParticleDistribution P R O T O T Y P E S 


EXTERN IMPORT  NSMCSchemaParticleDistribution*  NSMCSchemaParticleDistribution_Create();


EXTERN IMPORT  void  NSMCSchemaParticleDistribution_Init(NSMCSchemaParticleDistribution*);


EXTERN IMPORT  unsigned long  NSMCSchemaParticleDistribution_GetCount(NSMCSchemaParticleDistribution*);


EXTERN IMPORT  NSMCSchemaParticle*  NSMCSchemaParticleDistribution_GetNext(
                    NSMCSchemaParticleDistribution*, 
                    NSMCSchemaParticle*);


EXTERN IMPORT  unsigned long  NSMCSchemaParticleDistribution_Add(
                    NSMCSchemaParticleDistribution*, 
                    NSMCSchemaParticle*);


EXTERN IMPORT  NSMCSchemaParticleDistribution*  NSMPSCSchemaParticleDistribution_Delete(NSMCSchemaParticleDistribution*);













// M E T E O R L O G I C A L    D E S C R I P T I O N s

 struct NSMMetParameter
{
	boolean specified;
	double value;
}; 
typedef  struct NSMMetParameter  NSMMetParameter;




struct NSMMetEnvironment
{
	// Time parameters;
	long  year;

	long  month;

	long  day;

	double  julianDay;				// Current time in days

	double  timeStepInDays;			// Current time step in days


	// Site specific
	double  latitudeDegrees;		// degrees - no minutes or seconds

	char  latitudeDirection;		// 'N' for north or 'S' for south


	double  longitudeDegrees;		// degrees - no minutes or seconds

	char  longitudeDirection;		// 'E' for east or 'W' for west


	double  elevationMSL;			// Elevation above mean sea level. m


	// Temperature parameters
	NSMMetParameter  airTemp;		// Measured air temperature. Units: C
	
	NSMMetParameter  dewPointTemp;	// Measured dew point temp. Units: C
	
	NSMMetParameter  averageAirTemp;		// Daily average air temperature. Units: C
	
	NSMMetParameter  temperatureMeasurementHeight; // Distance above soil surface for temp. measurement. m.


	// Pressure
	NSMMetParameter  barometricPressure;		// kPa  (1013.2 mbar == 101.32 kPa)


	// Wind speed
	NSMMetParameter  windSpeed;				// m/s
	
	NSMMetParameter  windSpeedMeasurementHeight;	// Distance above soil surface for wind measurement. m.


	// Rainfall parameters
	NSMMetParameter  rainfallRate;			// Rainfall rate. Units: meters/day
	
	NSMMetParameter  rainfallNitrogenConc;	// Conc of nitrogen in rainfall. Units: (mg N)/Liter


	// Solar radiation parameters
	NSMMetParameter  solarRadiation;		// Measured radiation Units: MJ/(m^2.h)

	// To convert from Watts/m^2 to MJ/(m^2.h), multiply by 0.0036


	NSMMetParameter  albedo;

	NSMMetParameter  cloudCoverFraction;		// Fraction of sky cover by clouds


	// Snow parameters
	NSMMetParameter  snowCoverWaterContent;	// Water content of snow cover. Units: millimeters

	struct NSMMetEnvironment* metEnvironment;
};

typedef struct NSMMetEnvironment  NSMMetEnvironment;



//NSMMetEnvironment   P R O T O T Y P E S


EXTERN IMPORT  NSMMetEnvironment*  NSMMetEnvironment_Create();

EXTERN IMPORT  NSMMetEnvironment* NSMMetEnvironment_Delete(NSMMetEnvironment* me);

EXTERN IMPORT  void  NSMMetEnvironment_Update(NSMMetEnvironment* me);


EXTERN IMPORT  double  NSMMetEnvironment_CalculateWaterTemp(NSMMetEnvironment* me);









// L A N D   U S E   P R O P E R T I E S 

struct NSMLandUseProperties
{

	double soilDetachability;		// Soil Detachability parameter. Units: [kg/m^3]

};	// TKG:3-2009:  NSM_RV6.0

typedef  struct NSMLandUseProperties  NSMLandUseProperties;



// NSMLandUseProperties prototypes


EXTERN IMPORT  NSMLandUseProperties*  NSMLandUseProperties_Create();


EXTERN IMPORT  NSMLandUseProperties*  NSMLandUseProperties_Delete(NSMLandUseProperties* lup);


EXTERN IMPORT  void NSMLandUSeProperties_Init(NSMLandUseProperties* lup);










// N U T R I E N T   P R O P E R T I E S

typedef struct 
{
	// Dissolved Transfer Properties
    double NO2_Km;
	double NO3_Km;
	double NH4_Km;
	double ON_Km;
	double OP_Km;
	double DP_Km;

	// Sorption Properties
	double NO2_KOC;
	double NO3_KOC;
	double NH4_KOC;
	double ON_KOC;
	double OP_KOC;
	double DP_KOC;

} NSMNutrientProperties;			//TKG:3-2008:  NSM_RV3.2
	


// NSMNutrientProperties  Prototypes


EXTERN IMPORT  NSMNutrientProperties*  NSMNutrientProperties_Create();


EXTERN IMPORT  NSMNutrientProperties*  NSMNutrientProperties_Delete(NSMNutrientProperties*);


EXTERN IMPORT  void  NSMNutrientProperties_Init(NSMNutrientProperties* se);














// A Q U A T I C   E N V I R O N M E N T 

typedef struct
	{
	double rc;      			/* reaction_coefficient */
	double etcc;				/* Temperature Correction Coefficient (theta) */
	} NSM_COEFFICIENT;





typedef struct
	{
	NSM_COEFFICIENT alpha0;		//TC: alpha0 := Fraction of chlorophil-A to Algal biomass.
								//Units: (ug-Chl-A/mg-A). Range: {10-100} 

	NSM_COEFFICIENT alpha1;		//TC: alpha1 := Fraction of algal biomass that is nitrogen.
								//Units: (mg-N/mg-A). Range: {0.07-0.09} 

	NSM_COEFFICIENT alpha2;		//TC: alpha2 := Fraction of algal biomass that is Phosphorus.
								//Units: (mg-N/mg-A). Range: {0.01-0.02} 
	 
	NSM_COEFFICIENT alpha3;		//TC: alpha3 := O2 production per unit of algal growth.
								//Units: (mg-O/mg-A). Range: {1.4-1.8} 

	NSM_COEFFICIENT alpha4;		//TC: alpha4 := O2 uptake per unit algae respired.
								//Units: (mg-O/mg-A). Range: {1.6-2.3} 

	NSM_COEFFICIENT alpha5;		//TC: alpha5 := O2 uptake per unit NH3 oxidized.
								//Units: (mg-O/mg-NH3). Range: {3.0-4.0} 

	NSM_COEFFICIENT alpha6;		//TC: alpha6 := O2 uptake per unit NO2 oxidized.
								//Units: (mg-O/mg-NO2). Range: {1.0-1.14} 



	NSM_COEFFICIENT beta1;		//TVSV: beta1 := Ammonia decay (rate constant for oxidation of NH3 (ammonia) to NO2): 
								//Units: (1/day). Range: {0.1-1.0} 

	NSM_COEFFICIENT beta2;		//TVSV: beta2 := Nitrite decay (rate constant for the biological oxidation of NO2 to NO3:
								//Units: (1/day). Range: {0.2-2.0} 

	NSM_COEFFICIENT beta3;		//TVSV: beta3 := Nitrite decay (rate constant for the hydrolysis of organic N to ammonia). 
								//Units: (1/day). Range: {0.02-0.4} 

	NSM_COEFFICIENT beta4;		//TVSV: beta4 := Organic Phosphorus decay (rate constant for the decay of organic-P to dissolved P).
								//Units: (1/day). Range: {0.01-0.7} 



	NSM_COEFFICIENT K1;			//TVSV: K1 := Carbonaceous Deoxygenation rate constant.
								//Units: (1/day). Range: {0.02-3.4} 

	NSM_COEFFICIENT K2;			//TVSV: K2 := Reaeration Rate Coefficient.
								//Units: (1/day). Range: {0.0-100.0} 

	NSM_COEFFICIENT K3;			//TVSV: K3 := CBOD settling rate.
								//Units: (1/day). Range: {-0.36-0.36} 

	NSM_COEFFICIENT K4;			//TVSV: K4 := SOD (Benthic Oxygen) uptake. 
								//Units: (mg-O/ft^2-day). Range: {variable}



	NSM_COEFFICIENT mu;			//TV: mu := Algal growth rate. 
								//Units: (1/day). Range: {1.0-3.0} 

	NSM_COEFFICIENT rho;		//TV: rho := Algal respiration rate. 
								//Units (1/day). Range: {0.05,0.5} 



	NSM_COEFFICIENT sigma1;		//TVSV: sigma3 := Algal settling rate. 
								//Units (ft/day).. Range: {variable} 

	NSM_COEFFICIENT sigma2;		//TVSV: sigma2 := Benthos source rate for dissolved Phosphorus. 
								//Units: (mg/ft^2-day). Range: {variable} 

	NSM_COEFFICIENT sigma3;		//TVSV: sigma3 := Benthos source rate for ammonia Nitrogen. 
								//Units: (mg/ft^2-day). Range: {variable} 

	NSM_COEFFICIENT sigma4;		//TVSV: sigma4 := Organic Nitrogen settling rate: 
								//Units: (1/day). Range: {0.001-0.1} 

	NSM_COEFFICIENT sigma5;		//TVSV: sigma5 := Organic Phosphorus settling rate: 
								//Units: (1/day). Range: {0.001-0.1} 


	double Kl;					//Michaelis-Menton half-saturation constant for light. 
								//Range: {0.02-0.10}. Units: (Btu/(ft^2-min)).

	double Kn;					//Michaelis-Menton half-saturation constant for nitrogen.
								//Range: {0.01-0.02}. Units: (mg-N/L).	

	double Kp;					//Michaelis-Menton half-saturation constant for phosphorus. 
								//Range: {0.01-0.20}. Units: (mg-P/L).    


	double Pn;					//Algal perference factor for ammonia. 
								//Range: {0.0-1.0}. Units: nondimensional.	


	double lambda0;				//Non-algal self-shading coefficient. 
								//Range: {variable}. Units: (1/ft). 

	double lambda1;				//Linear algal self-shading coefficient. 
								//Range: {0.002-0.02}. Units: ( (1/ft) / (ug Chl-A/L) ) 

	double lambda2;				//Nonlinear algal self-shading coefficient. 
								//Range: {0.0165}. Units: ( (1/ft) / (ug Chl-A/L)^(2/3) )

	double min_do_conc;			//Min [DO] for kinetics calculations 

	
} NSM_AQUATIC_ENVIRONMENT;



	
// Aquatic Environment prototypes


EXTERN IMPORT  NSM_AQUATIC_ENVIRONMENT*  NSM_AQUATIC_ENVIRONMENT_Create();


EXTERN IMPORT  NSM_AQUATIC_ENVIRONMENT*  NSM_AQUATIC_ENVIRONMENT_Init(NSM_AQUATIC_ENVIRONMENT*);


EXTERN IMPORT  NSM_AQUATIC_ENVIRONMENT*  NSM_AQUATIC_ENVIRONMENT_Delete(NSM_AQUATIC_ENVIRONMENT*);















// C H A N N E L  & O V E R L A N D  S P E C F I C   N U T R I E N T   I N F O



enum  NSM_PARTITIONABLE_ENUM
{
    NSM_IS_PARTITIONABLE_Enum,      // Use if: Want partitioning between dissolved, bound and sediment-sorbed
    NSM_IS_FULLY_DISSOLVED_Enum,    // Use if: Want only dissolved form
    NSM_IS_FULLY_SORBED_Enum        // Use if: Want only sediment-sorbed form
};



enum  NSM_SPECIES_ENUM  
{
	NSM_NO2_Enum = 0,
	NSM_NO3_Enum, 
	NSM_NH4_Enum, 
	NSM_ORGANIC_NITROGEN_Enum, 
	NSM_ORGANIC_PHOSPHORUS_Enum,
	NSM_DISSOLVED_PHOSPHORUS_Enum,
	NSM_ALGAE_Enum,
	NSM_CBOD_Enum,
	NSM_DISSOLVED_OXYGEN_Enum,
	NSM_Soil_N_organicActive_Enum = 100,
	NSM_Soil_N_organicStable_Enum,
	NSM_Soil_P_organicActive_Enum,
	NSM_Soil_P_organicStable_Enum,
	NSM_Soil_P_mineralActive_Enum,
	NSM_Soil_P_mineralStable_Enum,
    NSM_NULL_Enum
};



//  NSM_SPECIES_ENUM   F U N C T I O N   P R O T O T Y P E S  

EXTERN IMPORT  int  NSM_SPECIES_ENUM_IsChannelSpecie(enum NSM_SPECIES_ENUM specieEnum);

EXTERN IMPORT  int  NSM_SPECIES_ENUM_IsOverlandSpecie(enum NSM_SPECIES_ENUM specieEnum);

EXTERN IMPORT  int  NSM_SPECIES_ENUM_IsSoilSpecie(enum NSM_SPECIES_ENUM specieEnum);

EXTERN IMPORT  enum  NSM_PARTITIONABLE_ENUM  NSM_SPECIES_ENUM_IsPartitionable(enum NSM_SPECIES_ENUM specieEnum);






typedef struct
{
	double  conc;
	double  mass;
	double  dmdt;

    enum NSM_PARTITIONABLE_ENUM   partitionable; // New in 4.6 & later

	double  fractionDissolved;
	double  fractionBound;

	unsigned long  numberSorbingParticles;
    double*  fractionPerSorbingParticle;
	double  fsum; 	// Sorption validation check. Should be very near 1.0.
}  NSM_AQUATIC_SPECIES;





#define NSM_NUMBER_AQUATIC_SPECIES 9

typedef enum NSM_AQUATIC_SPECIES_ID 
{
    NSM_AQUATIC_NO2, 
    NSM_AQUATIC_NO3, 
    NSM_AQUATIC_NH4, 
    NSM_AQUATIC_ON,  
    NSM_AQUATIC_OP,  
    NSM_AQUATIC_DP, 
    NSM_AQUATIC_Alg, 
    NSM_AQUATIC_CBOD,
    NSM_AQUATIC_DO,
    NSM_AQUATIC_DOC,            // Per request
    NSM_AQUATIC_DOC_FOC,        // Per request
    NSM_AQUATIC_DISPERSION      // Per request
} NSM_AQUATIC_SPECIES_ID;       // Not used by NSM - juust for interfacing




typedef struct 
{
  long numberSpecies;
  NSM_AQUATIC_SPECIES aquatic_NO2;
  NSM_AQUATIC_SPECIES aquatic_NO3;
  NSM_AQUATIC_SPECIES aquatic_NH4;
  NSM_AQUATIC_SPECIES aquatic_ON;
  NSM_AQUATIC_SPECIES aquatic_OP;
  NSM_AQUATIC_SPECIES aquatic_DP;
  NSM_AQUATIC_SPECIES aquatic_Alg;
  NSM_AQUATIC_SPECIES aquatic_CBOD;
  NSM_AQUATIC_SPECIES aquatic_DO;
  NSM_AQUATIC_SPECIES* list[NSM_NUMBER_AQUATIC_SPECIES];
} NSM_AQUATIC_SPECIES_NET;














// NSMChannelCell  T Y P E D E F S

typedef struct 
{
	double waterTemp;		// Temperature of water in cell.  Units: C

	double depth;			// Depth of water in cell.  Units: m

    double minimumDepth;    // Minimum depth to run kinetics

	double vol;				// Volume of water in cell.  Units: m^3

	NSM_AQUATIC_SPECIES_NET aqSpecies;

	double dissolvedOrganicCarbon;		// Conc (M/L^3)

	double dissolvedOrganicCarbon_FOC; 	// Logically = 1.0.

} NSMChannelCell;




// NSMChannelCell Prototypes


EXTERN IMPORT  NSMChannelCell*  NSMChannelCell_Create();


EXTERN IMPORT  void  NSMChannelCell_Init(NSMChannelCell* pCell);


EXTERN IMPORT  NSMChannelCell*  NSMChannelCell_Delete(NSMChannelCell* pCell);


EXTERN IMPORT  void  NSMChannelCell_EnableSystem(int enable);


EXTERN IMPORT  void  NSMChannelCell_Kinetics(
                    NSMChannelCell* pCell, 
                    NSMMetEnvironment* me, 
                    NSM_AQUATIC_ENVIRONMENT* ae);


EXTERN IMPORT double NSMChannelCell_CalculateMaximumDissolvedOxygenConc(
                    double waterTempC);


EXTERN IMPORT  void  NSMChannelCell_SetMinimumDepth(
                    NSMChannelCell* pCell, 
                    double depth);



// Two techniques for calculating specie partitioning

EXTERN IMPORT  void NSMChannelCell_CalculateSedimentPartitioning(
                    NSMChannelCell* pCell, 
                    NSMNutrientProperties* np, 
                    NSMSedimentEnvironment* se);        // If using NSMSedimentEnvironment


EXTERN IMPORT  void  NSMChannelCell_CalculateSpecieSedimentPartitioning(
                    NSMChannelCell* pCell, 
                    enum NSM_SPECIES_ENUM specie, 
                    NSMNutrientProperties* np, 
                    NSMSedimentEnvironment* se);        // If using NSMSedimentEnvironment




EXTERN IMPORT  void  NSMChannelCell_CalculateSpeciePartitioning(
                    NSMChannelCell* pCell, 
                    enum NSM_SPECIES_ENUM specie, 
                    NSMNutrientProperties* np, 
                    NSMCSchemaParticleDistribution* pd);  // If using NSMCSchemaParticleDistribution












// NSMOverlandCell  T Y P E D E F S

typedef struct 
{
	double waterTemp;		// Temperature of water in cell.  Units: C

	double depth;			// Depth of water in cell.  Units: m

    double minimumDepth;    // Minimum depth to run kinetics

	double vol;				// Volume of water in cell.  Units: m^3

	NSM_AQUATIC_SPECIES_NET aqSpecies;

	double dissolvedOrganicCarbon;

	double dissolvedOrganicCarbon_FOC; 	// Logically = 1.0.

} NSMOverlandCell;








// NSMOverlandCell Prototypes


EXTERN IMPORT  NSMOverlandCell*  NSMOverlandCell_Create();


EXTERN IMPORT  void  NSMOverlandCell_Init(NSMOverlandCell* pCell);


EXTERN IMPORT  NSMOverlandCell*  NSMOverlandCell_Delete(NSMOverlandCell* pCell);


EXTERN IMPORT  void NSMOverlandCell_EnableSystem(int enable);


EXTERN IMPORT  void NSMOverlandCell_Kinetics(
                    NSMOverlandCell* pCell, 
                    NSMMetEnvironment* me, 
                    NSM_AQUATIC_ENVIRONMENT* ae);


EXTERN IMPORT double NSMOverlandCell_CalculateMaximumDissolvedOxygenConc(
                    double waterTempC);


EXTERN IMPORT  void NSMOverlandCell_SetMinimumDepth(
                    NSMOverlandCell* pCell, 
                    double depth);



// Two techniques for calculating specie partitioning

EXTERN IMPORT  void NSMOverlandCell_CalculateSedimentPartitioning(
                    NSMOverlandCell* pCell, 
                    NSMNutrientProperties* np, 
                    NSMSedimentEnvironment* se);            // If using NSMSedimentEnvironment


EXTERN IMPORT  void NSMOverlandCell_CalculateSpecieSedimentPartitioning(
                    NSMOverlandCell* pCell, 
                    enum NSM_SPECIES_ENUM specie, 
                    NSMNutrientProperties* np,
                    NSMSedimentEnvironment* se);            // If using NSMSedimentEnvironment



EXTERN IMPORT  void NSMOverlandCell_CalculateSpeciePartitioning(
                    NSMOverlandCell* pCell, 
                    enum NSM_SPECIES_ENUM specie, 
                    NSMNutrientProperties* np, 
                    NSMCSchemaParticleDistribution* pd);    // If using NSMCSchemaParticleDistribution









/*  F U N C T I O N   P R O T O T Y P E S    */

EXTERN IMPORT  void  NSM_init();













// P L A N T - S O I L     D E S C R I P T I O N 


#define MAX_NUMBER_SOIL_PARTICLES 25		//TKG:NOTE: Come back when I have more time and eliminate this.


#define NSM_NUMBER_SOIL_SPECIES 11



struct NSMPSPlantEnvironment
{
	char name[200];

	double coverFraction;

    struct NSMPSPlantEnvironment*  next;
};

typedef struct NSMPSPlantEnvironment NSMPSPlantEnvironment;







struct NSMPSCSchemaFlora
{
	unsigned long numberPlants;

	NSMListNode* plantList;
};

typedef struct NSMPSCSchemaFlora NSMPSCSchemaFlora;







struct NSMPSPool
{
	double  mass;

	double  fractionDissolved;

	double  fractionSolid;

	double  fsum;
};

typedef struct NSMPSPool NSMPSPool;







struct NSMPSCarbonEnvironment
{
	double fractionCarbonInSoil;

	NSMPSPool organicActive;

	NSMPSPool organicStable;
};

typedef struct NSMPSCarbonEnvironment NSMPSCarbonEnvironment;








struct NSMPSNitrogenEnvironment
{
	NSMPSPool NH4;

	NSMPSPool NO3;

	NSMPSPool organicFresh;

	NSMPSPool organicActive;

	NSMPSPool organicStable;
};

typedef struct NSMPSNitrogenEnvironment NSMPSNitrogenEnvironment;








struct NSMPSPhosphorusEnvironment
{
	NSMPSPool organicFresh;

	NSMPSPool organicActive;

	NSMPSPool organicStable;

	NSMPSPool mineralActive;

	NSMPSPool mineralStable;

	NSMPSPool mineralSolution;

	double phosphorusAvailabilityIndex;			// Fraction fertilizer P in solution 
												// after period of fast reaction
};

typedef struct NSMPSPhosphorusEnvironment NSMPSPhosphorusEnvironment;














struct NSMPSSoilLayerEnvironment
{
	long  index;					// Index of layer in soil layer stack

	double thickness;				// Soil layer thickness.  Units:m

	double depth;					// Depth to layer bottom. Units:m

	double area;					// Area of layer. Units:m^2 

	double density;					// Soil density. Units:kg/m^3

	double porosity;				// Porosity. Dimensionless.


	// 7-2006: Changing unit dimension from meters water to [water volume]/[soil volume]
	//double waterContent;			// WaterContent. Units:m
	//double fieldCapacity;			// Field capacity of soil layer. Units:m
	//double wiltingPoint;			// Wilting point of soil layer. Units:m

	double volumetricWaterContent;	// WaterContent. Units:dimensionless
	double volumetricFieldCapacity;	// Field capacity of soil layer. Units:dimensionless
	double volumetricWiltingPoint;	// Wilting point of soil layer. Units:dimensionless



	double floraWaterDemand;		// Water demanded by flora community. Units: m


	double pH;						// Units: dimensionless

	double soilTemp;				// Soil layer temperature. Units:C


	double fractionCarbonInSoil;	// Soil Organic Matter fraction

    double dissolvedOrganicCarbon;	// Pore water conc. (Till carbon cycle is full-featured)

	double residue;					// Mass of plant residue in layer. Units:Kg


	NSMPSCarbonEnvironment  carbonEnvironment;

	NSMPSNitrogenEnvironment  nitrogenEnvironment;

	NSMPSPhosphorusEnvironment  phosphorusEnvironment;
};

typedef struct NSMPSSoilLayerEnvironment NSMPSSoilLayerEnvironment;
	





struct NSMPSCSchemaSoil
{
	unsigned long numberSoilLayers;

	NSMListNode* soilLayerList;
};

typedef struct NSMPSCSchemaSoil NSMPSCSchemaSoil;








struct NSMPSCellEnvironment
{
	// Plant environment
	unsigned long  numberPlants;

	NSMPSPlantEnvironment**  plantEnvironments;


	// Soil environment
	unsigned long  numberSoilLayers;

	NSMPSSoilLayerEnvironment**  soilLayerEnvironments;
};

typedef struct NSMPSCellEnvironment NSMPSCellEnvironment;









struct NSMCPSCell
{
	NSMPSCellEnvironment*  cellEnvironment;	// Each cell has it's own NSMPSCellEnvironment

	struct NSMCPSCell*  cell;
};

typedef struct NSMCPSCell NSMCPSCell;








// E X P O R T E D    P R O T O T Y P E S



// NSMPSPlantEnvironment prototypes

EXTERN IMPORT  NSMPSPlantEnvironment*  NSMPSPlantEnvironment_Create(
                    char* name, 
                    double coverFraction);


EXTERN IMPORT  NSMPSPlantEnvironment*  NSMPSPlantEnvironment_Delete(NSMPSPlantEnvironment* pe);


EXTERN IMPORT  void  NSMPSPlantEnvironment_Init(
                    NSMPSPlantEnvironment* pe, 
                    char* name, 
                    double coverFraction);


EXTERN IMPORT  NSMPSPlantEnvironment*  NSMPSPlantEnvironment_Clone(NSMPSPlantEnvironment* pe);






// NSMPSCSchemaFlora  prototypes


EXTERN IMPORT  NSMPSCSchemaFlora*  NSMPSCSchemaFlora_Create();


EXTERN IMPORT  void  NSMPSCSchemaFlora_Init(NSMPSCSchemaFlora* fe);


EXTERN IMPORT  unsigned long  NSMPSCSchemaFlora_Add(
                    NSMPSCSchemaFlora* fe, 
                    NSMPSPlantEnvironment* pe);


EXTERN IMPORT  unsigned long  NSMPSCSchemaFlora_AddList(
                    NSMPSCSchemaFlora* fe, 
                    int numberPlants, 
                    NSMPSPlantEnvironment** plantEnvs);


EXTERN IMPORT  unsigned long  NSMPSCSchemaFlora_AddPlantEnvs(
                    NSMPSCSchemaFlora* fe, 
                    int numberPlants, 
                    ...);


EXTERN IMPORT  NSMPSCSchemaFlora*  NSMPSCSchemaFlora_Delete(NSMPSCSchemaFlora* fe);


EXTERN IMPORT  unsigned long  NSMPSCSchemaFlora_GetCount(NSMPSCSchemaFlora* fe);


EXTERN IMPORT  NSMPSPlantEnvironment*  NSMPSCSchemaFlora_GetNext(
                    NSMPSCSchemaFlora* fe, 
                    NSMPSPlantEnvironment* last);


EXTERN IMPORT  unsigned long  NSMPSCSchemaFlora_GetList(
                    NSMPSCSchemaFlora* fe, 
                    unsigned long numberPlantsRequested, 
                    NSMPSPlantEnvironment** list);










// NSMPSSoilLayerEnvironment prototypes


EXTERN IMPORT  NSMPSSoilLayerEnvironment*  NSMPSSoilLayerEnvironment_Create();


EXTERN IMPORT  NSMPSSoilLayerEnvironment*  NSMPSSoilLayerEnvironment_Delete(NSMPSSoilLayerEnvironment* p);


EXTERN IMPORT  void  NSMPSSoilLayerEnvironment_Init(NSMPSSoilLayerEnvironment* sle);


EXTERN IMPORT  NSMPSSoilLayerEnvironment*  NSMPSSoilLayerEnvironment_Clone(NSMPSSoilLayerEnvironment* sle);








// NSMPSCSchemaSoil prototypes


EXTERN IMPORT  NSMPSCSchemaSoil*  NSMPSCSchemaSoil_Create();


EXTERN IMPORT  void  NSMPSCSchemaSoil_Init(NSMPSCSchemaSoil* env);


EXTERN IMPORT  unsigned long  NSMPSCSchemaSoil_Add(
                    NSMPSCSchemaSoil* se, 
                    NSMPSSoilLayerEnvironment* sle);


EXTERN IMPORT  unsigned long  NSMPSCSchemaSoil_AddList(
                    NSMPSCSchemaSoil* sc, 
                    int numberSoilLayers, 
                    NSMPSSoilLayerEnvironment** soilLayerEnvs);


EXTERN IMPORT  unsigned long  NSMPSCSchemaSoil_AddSoilLayerEnvs(
                    NSMPSCSchemaSoil* fe, int numberSoilLayers, 
                    ...);


EXTERN IMPORT  NSMPSCSchemaSoil*  NSMPSCSchemaSoil_Delete(NSMPSCSchemaSoil* sse);


EXTERN IMPORT  unsigned long  NSMPSCSchemaSoil_GetCount(NSMPSCSchemaSoil* sse);


EXTERN IMPORT  NSMPSSoilLayerEnvironment*  NSMPSCSchemaSoil_GetNext(
                    NSMPSCSchemaSoil* sse, 
                    NSMPSSoilLayerEnvironment* last);


EXTERN IMPORT  unsigned long  NSMPSCSchemaSoil_GetList(
                    NSMPSCSchemaSoil* sc, 
                    unsigned long numberSoilLayersRequested,
                    NSMPSSoilLayerEnvironment** list);








// NSMCellEnvironment prototypes


EXTERN IMPORT  NSMPSCellEnvironment*  NSMPSCellEnvironment_Create(long numberPlants, long numberSoilLayers);


EXTERN IMPORT  NSMPSCellEnvironment*  NSMPSCellEnvironment_CreateViaSchemas(
                    NSMPSCSchemaFlora* fse, 
                    NSMPSCSchemaSoil* sse);


EXTERN IMPORT  NSMPSCellEnvironment*  NSMPSCellEnvironment_Delete(NSMPSCellEnvironment* ce);


EXTERN IMPORT  NSMPSCellEnvironment*  NSMPSCellEnvironment_Copy(NSMPSCellEnvironment* cea);






// NSMCell prototypes


EXTERN IMPORT  NSMCPSCell*  NSMCPSCell_Create(NSMPSCellEnvironment* cellEnv);


EXTERN IMPORT  NSMCPSCell*  NSMCPSCell_Delete(NSMCPSCell*);


EXTERN IMPORT  NSMPSCellEnvironment*  NSMCPSCell_Advance(
                    NSMMetEnvironment* me, 
                    NSMCPSCell* p);


EXTERN IMPORT  void  NSMCPSCell_Update(
                    NSMCPSCell* p);


EXTERN IMPORT  void  NSMCPSCell_UpdateAllSpecies(
                    NSMCPSCell* p);


EXTERN IMPORT  void  NSMCPSCell_UpdateAllSpeciesByLayer(
                    NSMCPSCell* p, 
                    int layer);


EXTERN IMPORT  void  NSMCPSCell_UpdateSpecieByLayer(
                    NSMCPSCell* p, 
                    enum NSM_SPECIES_ENUM specie, 
                    int layer);


EXTERN IMPORT  unsigned long  NSMCPSCell_GetNumberOfSoilLayers(
                    NSMCPSCell* p);


EXTERN IMPORT  unsigned long  NSMCPSCell_GetSoilWaterDemand(
                    NSMCPSCell* p, 
                    double* buffer);


EXTERN IMPORT  double NSMCPSCell_GetSoilNutrientMass(
                    NSMCPSCell* pCell, 
                    enum NSM_SPECIES_ENUM id, 
                    unsigned long layer);


EXTERN IMPORT  unsigned long  NSMCPSCell_GetSoilNutrientConc(
                    NSMCPSCell* pCell, 
                    enum NSM_SPECIES_ENUM id, 
                    double* buffer);


EXTERN IMPORT  double  NSMCPSCell_ApplyMassFluxToSoilSurface(
                    NSMCPSCell* pCell, 
                    enum NSM_SPECIES_ENUM id, 
                    double flux,
                    double time_step);






// Q2<-->PS   I N T E R F A C E    P R O T O T Y P E S


// Allow user to disable all overland/soil mass transfer

EXTERN IMPORT  void  NSMOverlandCell_EnableSoilExchange(int enable);



// New version 11-2008

EXTERN IMPORT  double NSMGetSedimentWaterMassFlux(
                    NSMOverlandCell* ovCell, 
                    NSMCPSCell* psCell,
                    enum NSM_SPECIES_ENUM specieEnum,
                    NSMNutrientProperties* np,
                    NSMLandUseProperties* lup, 
                    NSMMetEnvironment* me);


// 3-2009: Added to ease computational burden associated with using NSMMetEnvironmet struct!

EXTERN IMPORT  double NSMGetOverlandSoilMassFlux(
                    NSMOverlandCell* ovCell, 
                    NSMCPSCell* psCell,
                    enum NSM_SPECIES_ENUM specieEnum,
                    NSMNutrientProperties* np,
                    NSMLandUseProperties* lup, 
                    double rainfallRate);







//                        I N F O   &   N O T E S
//
//        INTER-COMPARTMENT SPECIES MAPPING
//
//        OVERLAND               SOIL                         INDEX
//        oc.NH4 <-------------> ne.NH4                          0
//        oc.NO3 <-------------> ne.NO3                          1
//        oc.ON  <-------------> ne.OrganicFresh                 2
//        oc.DP  <-------------> pe.MineralSolution              3
//        oc.OP  <-------------> pe.OrganicActive                4


















// NSMK    A D D I T I O N S    /3/2010/


    // <NSMKMetProperties> 

    struct NSMKMetParameter
    {
        int specified;      // [0|1]: 0:value not specified, 1: value specified
        double value;
    };
    typedef struct NSMKMetParameter NSMKMetParameter;


    struct NSMKMetProperties
    {
        // Time parameters;
        long  year;

        long  month;

        long  day;

        double  julianDay;                      // Current time in days

        double  timeStepInDays;	                // Current time step in days




        // Site specific
        double  latitudeDegrees;                // degrees - no minutes or seconds
    	
        char  latitudeDirection;                // 'N' for north or 'S' for south

        double  longitudeDegrees;               // degrees - no minutes or seconds

        char  longitudeDirection;               // 'E' for east or 'W' for west

        double  elevationMSL;                   // Elevation above mean sea level. m




        // Temperature parameters
        NSMKMetParameter  airTemp;               // Measured air temperature. Units: C
    	
        NSMKMetParameter  dewPointTemp;	        // Measured dew point temp. Units: C
    	
        NSMKMetParameter  averageAirTemp;        // Daily average air temperature. Units: C
    	
        NSMKMetParameter  temperatureMeasurementHeight; // Distance above soil surface for temp. measurement. m.




        // Pressure
        NSMKMetParameter  barometricPressure;    // kPa  (1013.2 mbar == 101.32 kPa)




        // Wind speed
        NSMKMetParameter  windSpeed;             // m/s
    	
        NSMKMetParameter  windSpeedMeasurementHeight;   // Distance above soil surface for wind measurement. m.




        // Rainfall parameters
        NSMKMetParameter  rainfallRate;	        // Rainfall rate. Units: meters/day
    	
        NSMKMetParameter  rainfallNitrogenConc;  // Conc of nitrogen in rainfall. Units: (mg N)/Liter




        // Solar radiation parameters
        NSMKMetParameter  solarRadiation;        // Measured radiation Units: MJ/(m^2.h)

        // Watts/m^2 * 0.0036 == 1 MJ/(m^2.h)
        // langleys/min. x 697.32 == watts/meter2 
        // langleys/min. x 4.1855 == joules/centimeter2 min. 
        // langleys/min. x 3.692 == BTU/foot2 min
        // one langley == 4.184E+04 joules per square meter


        NSMKMetParameter  albedo;

        NSMKMetParameter  cloudCoverFraction;   // Fraction of sky cover by clouds



        // Snow parameters
        NSMKMetParameter snowCoverWaterContent; // Water content of snow cover. Units: millimeters

        //struct NSMKMetProperties* metProperties;



        double CO2PartialPressure;	        // Units: [atm]
    };
    typedef struct NSMKMetProperties NSMKMetProperties;


    EXTERN IMPORT  NSMKMetProperties*  NSMKMetProperties_Create();

    EXTERN IMPORT  NSMKMetProperties*  NSMKMetProperties_Delete(NSMKMetProperties* me);

    // </NSMKMetProperties> 









    // <LandUseProperties>

    struct NSMKLandUseProperties
    {
    	double soilDetachability;		// Soil Detachability parameter. Units: [kg/m^3]
    };
    typedef struct  NSMKLandUseProperties  NSMKLandUseProperties;

    EXTERN IMPORT  NSMKLandUseProperties*  NSMKLandUseProperties_Create();

    EXTERN IMPORT  NSMKLandUseProperties*  NSMKLandUseProperties_Delete(NSMKLandUseProperties* lup);

    // </LandUseProperties>










    // <NSMKPoolParameters>

    enum NSMKEnumBottomAlgaePhotolysisOptions
    {
        ZeroOrderOption,
        ZeroFirstOption,
    };


    enum NSMKEnumOxygenHydraulicReaerationOptions
    {
        OConner_Dobbins,
        Churchill,
        Owens_Gibbs,
        Thackston_Dawson,
        Covar
    };


    enum NSMKEnumOxygenWindReaerationOptions
    {
        None,
        Banks_Herrera,
        Wanninkhof
    };


    enum NSMKEnumOxidationAttenuationOptions
    {
        Exponential_attenuation,
        DO_half_saturation_attenuation
    };


    enum NSMKEnumPhotosynthesisModelOptions  
    {
        MichaelisMenton,
        SmithsModel,
        SteelesModel
    };

    
    struct NSMKPoolOptions
    {
        enum NSMKEnumBottomAlgaePhotolysisOptions  bottomAlgaePhotolysisEnum;

        enum NSMKEnumOxygenHydraulicReaerationOptions  oxygenHydraulicReaerationEnum;

        enum NSMKEnumOxygenWindReaerationOptions  oxygenWindReaerationEnum;

        enum NSMKEnumPhotosynthesisModelOptions  photosynthesisModelEnum;

        enum NSMKEnumOxidationAttenuationOptions  oxidationAttenuationEnum;
    };
    typedef struct NSMKPoolOptions NSMKPoolOptions;




    struct TD_COEFFICIENT
    {
        double rc;              /* reaction_coefficient */
        double etcc;            /* Temperature Correction Coefficient (theta) NEVER ZERO!*/
    };
    typedef struct TD_COEFFICIENT TD_COEFFICIENT;



    struct  NSMKPoolParameters
    {
        TD_COEFFICIENT Kdp;     // Kdp = Temperature-dependent Phytoplankton death rate [1/day]
        TD_COEFFICIENT Kdb;     // Kdb = Temperature-dependent bottom Algae death rate [1/day]

        TD_COEFFICIENT Keb;     // Keb = Phytoplankton light extinction coefficient [1/m]

        TD_COEFFICIENT Khn;     // khn = temperature-dependent PON hydrolysis rate coefficient [1/day]
        TD_COEFFICIENT Kmn;     // kmn = temperature-dependent DON mineralization rate coefficient [1/day]
        
        TD_COEFFICIENT Krp;     // krp = temperature-dependent phytoplankton respiration rate [1/day]
        TD_COEFFICIENT Knit;    // Knit = temperature-dependent nitrification rate coefficient [1/day]
        TD_COEFFICIENT Kdint;   // kdnit = temperature-dependent denitrification rate coefficient [1/day]

        double Khnxp;           // Khnxp = preference coefficient of phytoplankton for NH4+
        double Khnxb;           // Khnxb = preference coefficient of bottom algae for NH4+

        double Up;              // Up = phytoplankton photosynthesis rate [1/day]
        double Ub;              // Ub = bottom algae photosynthesis rate [1/day]

        double Kscdn;           // Kscdn = DOC half-saturation constant for denitrification [gC/m3] [M/L3]
      
        double Ksnp;            // Ksnp = Nitrogen half-saturation constant [ugN/L]
        double Kspp;            // Kspp = Phosphorus half-saturation constant [ugP/L]
        double Kscp;            // Kscp = Carbon half-saturation constant [mole/L]

        double CO2;	            // CO2 = Dissolved carbon dioxide concentration - [mole/L]
        double CO3;	            // CO3 = Dissolved bicarbonate concentration - [mole/L]

        double alpha_o;         // alpha_o = effect of particulate organic matter on light attenuation: [Langleys/mgD/m]
        double alpha_i;         // alpha_i = effect of inorganic suspended solids on light attenuation: [Langleys/mgD/m]
        double alpha_p;	        // alpha_p = Linear effect of chlorophyll on light attenuation: [Langleys/ugD/m]
        double alpha_pn;        // alpha_pn = non-linear effect of chlorphyll on light attenuation: [(Langleys/mgD)^2/3/m]

        double Klp;	            // Phytplankton light parameter
        double Kgp;	            // Maximum photosynthesis rate at temperature T

        double qn;              // Nitrogen cell quota [mgN/mgA]
        double q0n;             // Minimum Nitrogen cell quota [mgN/mgA]
        double qp;              // Phosphorus cell quota [mgP/mgA]
        double q0p;             // Minimum Phosphorus cell quota [mgP/mgA]

        double rho_mN;          // Maximum Nitrogen uptake rate [mgN/mgA/day]
        double Kqn;	            // Half-saturation constant for intracellular Nitrogen [mgN/mgA];

        double rho_mP;          // Maximum Phosphorus uptake rate [mgN/mgA/day]
        double Kqp;	            // Half-saturation constant for intracellular Phosphorus [mgP/mgA];
        double Kspb;            // Half-saturation constant for external Phosphorus [upP/L];

        TD_COEFFICIENT Krb;	    // Temperature-dependent Bottom Algae respiration rate [1/day]
        TD_COEFFICIENT Krd;	    // Temperature-dependent Bottom Algae death rate [1/day]
      
        TD_COEFFICIENT Kdt;	    // Temperature-dependent Detritus dissolution rate [1/day]

        TD_COEFFICIENT Khc;     // Temperature-dependent slow CBOD hydrolysis rate [1/day]

        TD_COEFFICIENT Kdcs;    // Temperature-dependent slow CBOD oxidation rate [1/day]

        TD_COEFFICIENT Kdc;	    // Temperature-dependent fast CBOD oxidation rate [1/day]

        TD_COEFFICIENT Kn;      // Temperature-dependent nitrification rate [1/day]

        TD_COEFFICIENT Khp;     // Temperature-dependent Organic Phosphorus hydrolysis rate [1/day]

        TD_COEFFICIENT Kmp;     // Temperature-dependent DOP mineralization rate [1/day]

        TD_COEFFICIENT Ka;      // Temperature-dependent DOX reaeration coefficient [ 1/day] !!!DYNAMIC!!!
                                // (.rc updated with unique value for every cell and time step)

        double Ff;              // Fraction of Detrital dissolution consumed via fast reacting CBOD [dimensionless]


        TD_COEFFICIENT Kac;     // Temperature-dependent CO2 reaeration coefficient [1/day]

        double Fic_co2;         // Fraction of total inorganic carbon in CO2

        TD_COEFFICIENT Kmc;     // Temperature-dependent DOC mineralization rate [1/day]

        TD_COEFFICIENT Khyc;    // Temperature-dependent POC hydrolysis rate [1/day]


        NSMKPoolOptions options;
    };
    typedef struct NSMKPoolParameters NSMKPoolParameters;


    EXTERN IMPORT  NSMKPoolParameters*  NSMKPoolParameters_Create();

    EXTERN IMPORT  NSMKPoolParameters*  NSMKPoolParameters_Delete(NSMKPoolParameters*);

    // </NSMKPoolParameters>





    // <NSMKParticleProperties>

    struct NSMKParticleProperties
    {
    	double  conc;                   // [g/m^3] 

     	double  FOC;			        // Fraction Organic Carbon

        double  velocity;               // Velocity [m/s]
                                        // Particle settling velocity if aquatic 
                                        // Particle resuspension velocity if benthic 
    };
    typedef struct NSMKParticleProperties NSMKParticleProperties;

    // </NSMKParticleProperties>





    // <NSMKSedProperties>

    struct NSMKSedProperties
    {
        int numberParticles;

        NSMKParticleProperties** particleProperties; // "C" type array of size: numberSediments

        double  totalParticleConc;                  // Sum of sediment particle concs

    	double  bindingEffectivenessCoefficient;    // De: DOC binding

    	double  particleInteractionParameter;       // DiToro's Vx parameter 
    };
    typedef struct  NSMKSedProperties  NSMKSedProperties;

    // </NSMKSedProperties>






    const int NSMK_POOL_COUNT = 18;


    enum NSMKPoolEnum
    {
        Enum_AB1,           // Algae - Bottom substrate group 1

        Enum_AB2,           // Algae - Bottom substrate group 2

        Enum_AP1,           // Algae - Phytoplankton group 1

        Enum_AP2,           // Algae - Phytoplankton group 2

        Enum_CBOD_fast,     // CBOD fast

        Enum_CBOD_slow,	    // CBOD slow

        Enum_DIC,           // Dissolved Inorganic Carbon  

        Enum_DIP,           // Dissolved Inorganic Phosphorus

        Enum_DOC,           // Dissolved Organic Carbon  

        Enum_DON,           // Dissolved Organic Nitrogen

        Enum_DOP,           // Dissolved Organic Phosphorus

        Enum_DOX,           // Dissolved Oxygen

        Enum_NH4,           // Ammonium Nitrogen

        Enum_NO3,           // Nitrate Nitrogen

        Enum_POC,           // Particulate Organic Carbon  

        Enum_POM,           // Particulate Organic Matter - Detritus

        Enum_PON,           // Particulate Organic Nitrogen

        Enum_POP            // Particulate Organic Phosphorus
    };





    // <PoolProperties>

    struct NSMKPoolAttribute
    {
        enum NSMKPoolEnum myEnum;


        // Transfer properties

        int transfers;                  // [0|1]: 0:pool does not transfer, 1:does transfer

        double Km;                      // Diffusive Transfer Velocity: [meter/sec]


        // Partitioning properties

        int partitions;                 // [0|1]: 0:pool does not partition, 1:does partition

        double KOC;                     // Partition coefficient


        // Settling properties

        int settles;                    // [0|1]: 0:pool does not settle, 1: does settle

        double settling_velocity;       // Pool settling velocity: [meter/sec]

        
        // Resuspension properties

        int resuspends;                 // [0|1]: 0:pool does not resuspend, 1: does resuspend

        double resuspension_velocity;   // Pool resuspension velocity: [meter/sec]
    };	
    typedef struct NSMKPoolAttribute NSMKPoolAttribute;
    	


    struct NSMKPoolAttributeGroup
    {
        int numberPools;

        NSMKPoolAttribute* pools; // "C" type array of NSMKPoolAttribute of length "numberPools"
    };
    typedef struct NSMKPoolAttributeGroup NSMKPoolAttributeGroup;


    EXTERN IMPORT  NSMKPoolAttributeGroup*  NSMKPoolAttributeGroup_Create();

    EXTERN IMPORT  NSMKPoolAttributeGroup*  NSMKPoolAttributeGroup_Delete(NSMKPoolAttributeGroup*);

    // </PoolProperties>





     

    struct NSMKPartitionDistribution
    {
        int numberParticles;

        double pool;                    // That which isn't bound or sorbed

        double bound;                   // That which is bound to carbon

        double sorbed;                  // That which is sorbed to all particles

        double* by_particle;            // That which is sorbed to particle n

        double sum;                     // Total of above
    };
    typedef struct NSMKPartitionDistribution NSMKPartitionDistribution;



    struct NSMKPool
    {
        enum NSMKPoolEnum  myEnum;

        double mass;                    // [g]

        double conc;                    // [g/m^3]
        
        double kinetic_dmdt;            // [g/sec]

        double transfer_dmdt;           // [g/sec]


        // Partitioning Distributions

        NSMKPartitionDistribution* settling_dmdt;      // [g/sec]

        NSMKPartitionDistribution* resuspension_dmdt;  // [g/sec]

        NSMKPartitionDistribution* partition;          // Non-dimensional partition fractions
    };
    typedef struct NSMKPool NSMKPool;










    struct NSMKAquaticReactor
    {
        int numberPools;
        NSMKPool** pools;

        double temperature;         // [Celcius]
        double volume;              // [meters^3]
        double depth;               // [meters]
        double meanWaterSpeed;      // [meters/sec]

        NSMKSedProperties sedProperties;
    };
    typedef struct NSMKAquaticReactor NSMKAquaticReactor;



    // NSMKAquaticReactor   P R O T O T Y P E S

    EXTERN IMPORT  NSMKAquaticReactor* NSMKAquaticReactor_Create(NSMKPoolAttributeGroup*, double area, double slope,  double hydraulicRadius, int numberParticles);


    EXTERN IMPORT  NSMKPool* NSMKAquaticReactor_PartitionByPool(NSMKAquaticReactor*, enum NSMKPoolEnum);

    EXTERN IMPORT  NSMKPool** NSMKAquaticReactor_Partition(NSMKAquaticReactor*);


    EXTERN IMPORT  void NSMKAquaticReactor_SettlingDmDt(NSMKAquaticReactor*, NSMKPoolAttributeGroup*);

    EXTERN IMPORT  void NSMKAquaticReactor_KineticDmDt(NSMKAquaticReactor*, NSMKMetProperties*, NSMKPoolParameters*);


    EXTERN IMPORT  NSMKPool** NSMKAquaticReactor_GetPools(NSMKAquaticReactor*);

    EXTERN IMPORT  NSMKPool* NSMKAquaticReactor_GetPool(NSMKAquaticReactor*, enum NSMKPoolEnum);





    struct NSMKBenthicReactor
    {
        int numberPools;
        NSMKPool** pools;

        double temperature;                 // [Celcius]
        double porosity;                    // Average porosity
        double density;                     // sediment density [g/m^3]
        double depth;                       // [meters]
        double volume;                      // [meters^3]

        NSMKSedProperties sedProperties;
    };
    typedef struct NSMKBenthicReactor NSMKBenthicReactor;


    // NSMKBenthicReactor   P R O T O T Y P E S

    EXTERN IMPORT  NSMKBenthicReactor* NSMKBenthicReactor_Create(NSMKPoolAttributeGroup*, double area, int numberParticles);


    EXTERN IMPORT  NSMKPool* NSMKBenthicReactor_PartitionByPool(NSMKBenthicReactor*, enum NSMKPoolEnum);

    EXTERN IMPORT  NSMKPool** NSMKBenthicReactor_Partition(NSMKBenthicReactor*);


    EXTERN IMPORT  void NSMKBenthicReactor_DmDtResuspension(NSMKBenthicReactor*, NSMKPoolAttributeGroup*);

    EXTERN IMPORT  void NSMKBenthicReactor_KineticDmDt(NSMKBenthicReactor*, NSMKMetProperties*,  NSMKPoolParameters*);


    EXTERN IMPORT  NSMKPool** NSMKBenthicReactor_GetPools(NSMKBenthicReactor*);

    EXTERN IMPORT  NSMKPool* NSMKBenthicReactor_GetPool(NSMKBenthicReactor*, enum NSMKPoolEnum);






    struct NSMKBenthicExchanger
    {
        NSMKAquaticReactor* pNSMKAquaticReactor;
        NSMKBenthicReactor* pNSMKbenthicReactor;
    };
    typedef struct NSMKBenthicExchanger NSMKBenthicExchanger;


    // NSMKBenthicExchanger P R O T O T Y P E S

    EXTERN IMPORT  NSMKBenthicExchanger* NSMKBenthicExchanger_Create(NSMKAquaticReactor*, NSMKBenthicReactor*);

    EXTERN IMPORT  void NSMKBenthicExchanger_DoParticleKineticsExchanges(NSMKBenthicExchanger*, NSMKPoolAttributeGroup*);

    EXTERN IMPORT  void NSMKBenthicExchanger_DoPoolTransfers(NSMKBenthicExchanger*, NSMKPoolAttributeGroup*);







     struct NSMKChannelCell
     {
        NSMKAquaticReactor* ar;

        NSMKBenthicReactor* br;

        NSMKBenthicExchanger* be;
     };
     typedef struct NSMKChannelCell NSMKChannelCell;

     
    // NSMKChannelCell   P R O T O T Y P E S

     EXTERN IMPORT  NSMKChannelCell* NSMKChannelCell_Create(NSMKPoolAttributeGroup*, double area, double slope,  double hydraulicRadius, int numberParticles);


     EXTERN IMPORT  void NSMKChannelCell_PartitionByPool(NSMKChannelCell*, enum NSMKPoolEnum, NSMKPoolAttributeGroup*);

     EXTERN IMPORT  void NSMKChannelCell_Partition(NSMKChannelCell*, NSMKPoolAttributeGroup*);


     EXTERN IMPORT  void NSMKChannelCell_DoParticleKinetics(NSMKChannelCell*, NSMKPoolAttributeGroup*);

     EXTERN IMPORT  void NSMKChannelCell_DoPoolKinetics(NSMKChannelCell*, NSMKMetProperties*, NSMKPoolParameters*);

     EXTERN IMPORT  void NSMKChannelCell_DoPoolTransfers(NSMKChannelCell*, NSMKPoolAttributeGroup*);


     EXTERN IMPORT  NSMKPool** NSMKChannelCell_GetAquaticPools(NSMKChannelCell*);

     EXTERN IMPORT  NSMKPool* NSMKChannelCell_GetAquaticPool(NSMKChannelCell*, enum NSMKPoolEnum);

     EXTERN IMPORT  NSMKPool** NSMKChannelCell_GetBenthicPools(NSMKChannelCell*);

     EXTERN IMPORT  NSMKPool* NSMKChannelCell_GetBenthicPool(NSMKChannelCell*, enum NSMKPoolEnum);


     EXTERN IMPORT  NSMKSedProperties* NSMKChannelCell_GetBenthicSedProperties(NSMKChannelCell*);

     EXTERN IMPORT  NSMKSedProperties* NSMKChannelCell_GetAquaticSedProperties(NSMKChannelCell*);


#endif
