#pragma once
#ifndef MET_ENVIRONMENT_H

#define MET_ENVIRONMENT_H


#include "LINKAGE.h"





typedef long boolean;


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









// E X P O R T E D    P R O T O T Y P E S


// Met functions

EXTERN EXPORT  NSMMetEnvironment* NSMMetEnvironment_Create();

EXTERN EXPORT  NSMMetEnvironment* NSMMetEnvironment_Delete(NSMMetEnvironment* me);

EXTERN EXPORT  void NSMMetEnvironment_Update(NSMMetEnvironment* me);





// Quasi Met functions

EXTERN EXPORT double NSMMetEnvironment_CalculateWaterTemp(NSMMetEnvironment* me);

#endif
