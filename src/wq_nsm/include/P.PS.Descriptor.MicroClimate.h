#pragma once


#include "C.PS.Environment.Cell.h"

#include "C.PS.Environment.Plant.h"

#include "C.PS.CSchema.Flora.h"

#include "C.PS.Environment.Nutrients.h"

#include "C.PS.Environment.Particle.h"

#include "C.PS.Environment.SoilLayer.h"

#include "C.PS.CSchema.Soil.h"





namespace NSMPS
{



	struct MicroClimateDescriptor
	{

		// C O N S T A N T S

		static const double PI;

		static const double SOLAR_CONSTANT;						// MJ/(m^2 minute)

		static const double STEFFAN_BOLTZMAN_CONSTANT;			// MJ/(m^2.h.K^4)

		static const double WATER_LATENT_HEAT_OF_VAPORIZATION;	// MJ/kg

		static const double DRY_ADIABATIC_LAPSE_RATE;			// Kelvin/m

		static const double GRAVITY;							// m/sec^2

		static const double AIR_MOLECULAR_WEIGHT;				// kg/mole

		static const double UNIVERSAL_GAS_CONSTANT;				// J/(mole.Kelvin)

		static const double EARTHS_ANGULAR_VELOCITY;			// Earths rotation rate:rad/hour



		// S P E C I F I E D    V A L U E S

		double latitudeDegrees;					// Degrees {North:+, South:-} 

		double longitudeDegrees;				// Degrees {East:+, West:-} 

		double elevationMSL;					// Elevation. Meters

		double timeStepInDays;					// Current timestep in DAYS.

		double julianDay;						// Day of the year (1-366)


		double temperature;						// Current temperature. C.

		double dewPointTemperature;				// Current dew point temp. C.

		double temperatureMeasurementHeight;	// Distance above soil surface for temp. measurement. m.


		double measuredWindSpeed;				// Wind speed in m/s.

		double windSpeedMeasurementHeight;		// Distance above soil surface for wind measurement. m.


		double barometricPressure;				// kPa

		double rainfallRate;					// m/day

		double rainfallNitrogenConcentration;	// Units: (mg N)/Liter

		double measuredSolarRadiation;			// MJ/(m^2.h)



		// C A L C U L A T E D   V A L U E S

		double latitude;						// Radians {North:+, South:-} 

		double longitude;						// Radians {East:+, West:-}

		double airDensity;						// Density of air. (kg/m^3);

		double extraterrestrialRadiation;		// MJ/(m^2.h)

		double clearSkyTotalGlobalSolarRadiation;	// MJ/(m^2.h)

		double saturationVaporPressure;			// kPa

		double vaporPressure;					// kPa

		// Slope of the saturation vapor pressure curve at mean air temp
		double saturationVaporPressureCurveSlope;	// kPa/Celsius

		double apparentClearSkyEmissivity;

		double cloudinessFactor;	


		double netShortWaveRadiation;			// MJ/(m^2.h)

		double netLongWaveRadiation;			// MJ/(m^2.h)

		double netRadiation;					// MJ/(m^2.h)	

		double soilHeatFluxDensity;				// MJ/(m^2.h)


		double psychrometricConstant;			// kPa/Celseus

		double modifiedPsychrometricConstant;	// kPa/Celseus

		double aerodynamicResistance;			// sec/m



		bool dayTime;		// TRUE if current time is during daylight
							// FALSE if current time is during nighttime


		double dayLengthHours;	// Length of day in hours
		
        private: MicroClimateDescriptor();

        public: MicroClimateDescriptor(NSMMetEnvironment* me);

        public: void update(NSMMetEnvironment* me);

        private: void Clear();

        private: void UpdateMetEnvironment(NSMMetEnvironment* me);
    };


}