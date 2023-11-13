#pragma once


#ifdef __USRDLL
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif


#ifdef __cplusplus
extern "C"
{
#endif

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


    EXPORT  NSMKMetProperties*  NSMKMetProperties_Create();

    EXPORT  NSMKMetProperties*  NSMKMetProperties_Delete(NSMKMetProperties* me);

    // </NSMKMetProperties> 









    // <LandUseProperties>

    struct NSMKLandUseProperties
    {
    	double soilDetachability;		// Soil Detachability parameter. Units: [kg/m^3]
    };
    typedef struct  NSMKLandUseProperties  NSMKLandUseProperties;

    EXPORT  NSMKLandUseProperties*  NSMKLandUseProperties_Create();

    EXPORT  NSMKLandUseProperties*  NSMKLandUseProperties_Delete(NSMKLandUseProperties* lup);

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


    EXPORT  NSMKPoolParameters*  NSMKPoolParameters_Create();

    EXPORT  NSMKPoolParameters*  NSMKPoolParameters_Delete(NSMKPoolParameters*);

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

    EXPORT  NSMKPoolAttributeGroup*  NSMKPoolAttributeGroup_Create();

    EXPORT  NSMKPoolAttributeGroup*  NSMKPoolAttributeGroup_Delete(NSMKPoolAttributeGroup*);

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

    EXPORT  NSMKAquaticReactor* NSMKAquaticReactor_Create(NSMKPoolAttributeGroup*, double area, double slope,  double hydraulicRadius, int numberParticles);


    EXPORT  NSMKPool* NSMKAquaticReactor_PartitionByPool(NSMKAquaticReactor*, enum NSMKPoolEnum);

    EXPORT  NSMKPool** NSMKAquaticReactor_Partition(NSMKAquaticReactor*);


    EXPORT  void NSMKAquaticReactor_SettlingDmDt(NSMKAquaticReactor*, NSMKPoolAttributeGroup*);

    EXPORT  void NSMKAquaticReactor_KineticDmDt(NSMKAquaticReactor*, NSMKMetProperties*, NSMKPoolParameters*);


    EXPORT  NSMKPool** NSMKAquaticReactor_GetPools(NSMKAquaticReactor*);

    EXPORT  NSMKPool* NSMKAquaticReactor_GetPool(NSMKAquaticReactor*, enum NSMKPoolEnum);





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

    EXPORT  NSMKBenthicReactor* NSMKBenthicReactor_Create(NSMKPoolAttributeGroup*, double area, int numberParticles);


    EXPORT  NSMKPool* NSMKBenthicReactor_PartitionByPool(NSMKBenthicReactor*, enum NSMKPoolEnum);

    EXPORT  NSMKPool** NSMKBenthicReactor_Partition(NSMKBenthicReactor*);


    EXPORT  void NSMKBenthicReactor_DmDtResuspension(NSMKBenthicReactor*, NSMKPoolAttributeGroup*);

    EXPORT  void NSMKBenthicReactor_KineticDmDt(NSMKBenthicReactor*, NSMKMetProperties*,  NSMKPoolParameters*);


    EXPORT  NSMKPool** NSMKBenthicReactor_GetPools(NSMKBenthicReactor*);

    EXPORT  NSMKPool* NSMKBenthicReactor_GetPool(NSMKBenthicReactor*, enum NSMKPoolEnum);






    struct NSMKBenthicExchanger
    {
        NSMKAquaticReactor* pNSMKAquaticReactor;
        NSMKBenthicReactor* pNSMKbenthicReactor;
    };
    typedef struct NSMKBenthicExchanger NSMKBenthicExchanger;


    // NSMKBenthicExchanger P R O T O T Y P E S

    EXPORT  NSMKBenthicExchanger* NSMKBenthicExchanger_Create(NSMKAquaticReactor*, NSMKBenthicReactor*);

    EXPORT  void NSMKBenthicExchanger_DoParticleKineticsExchanges(NSMKBenthicExchanger*, NSMKPoolAttributeGroup*);

    EXPORT  void NSMKBenthicExchanger_DoPoolTransfers(NSMKBenthicExchanger*, NSMKPoolAttributeGroup*);







     struct NSMKChannelCell
     {
        NSMKAquaticReactor* ar;

        NSMKBenthicReactor* br;

        NSMKBenthicExchanger* be;
     };
     typedef struct NSMKChannelCell NSMKChannelCell;

     
    // NSMKChannelCell   P R O T O T Y P E S

     EXPORT  NSMKChannelCell* NSMKChannelCell_Create(NSMKPoolAttributeGroup*, double area, double slope,  double hydraulicRadius, int numberParticles);


     EXPORT  void NSMKChannelCell_PartitionByPool(NSMKChannelCell*, enum NSMKPoolEnum, NSMKPoolAttributeGroup*);

     EXPORT  void NSMKChannelCell_Partition(NSMKChannelCell*, NSMKPoolAttributeGroup*);


     EXPORT  void NSMKChannelCell_DoParticleKinetics(NSMKChannelCell*, NSMKPoolAttributeGroup*);

     EXPORT  void NSMKChannelCell_DoPoolKinetics(NSMKChannelCell*, NSMKMetProperties*, NSMKPoolParameters*);

     EXPORT  void NSMKChannelCell_DoPoolTransfers(NSMKChannelCell*, NSMKPoolAttributeGroup*);


     EXPORT  NSMKPool** NSMKChannelCell_GetAquaticPools(NSMKChannelCell*);

     EXPORT  NSMKPool* NSMKChannelCell_GetAquaticPool(NSMKChannelCell*, enum NSMKPoolEnum);

     EXPORT  NSMKPool** NSMKChannelCell_GetBenthicPools(NSMKChannelCell*);

     EXPORT  NSMKPool* NSMKChannelCell_GetBenthicPool(NSMKChannelCell*, enum NSMKPoolEnum);


     EXPORT  NSMKSedProperties* NSMKChannelCell_GetBenthicSedProperties(NSMKChannelCell*);

     EXPORT  NSMKSedProperties* NSMKChannelCell_GetAquaticSedProperties(NSMKChannelCell*);

#ifdef __cplusplus
 }
#endif