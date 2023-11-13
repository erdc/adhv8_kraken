#pragma once
#ifndef C_PS_ENVIRONMENT_SOILLAYER_H
#define C_PS_ENVIRONMENT_SOILLAYER_H


#include "LINKAGE.h"

#include "C.METEnvironment.h"

#include "C.PS.Environment.Nutrients.h"

#include "C.PS.Environment.Particle.h"







struct NSMPSSoilLayerEnvironment
{
	long  index;					// Index of layer in soil layer stack

	double thickness;				// Soil layer thickness.  Units:m

	double depth;					// Depth to layer bottom.  Units:m

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
	





// S O I L   L A Y E R   E N V I R O N M E N T   P R O T O T Y P E S 


EXTERN EXPORT  NSMPSSoilLayerEnvironment* NSMPSSoilLayerEnvironment_Create();

EXTERN EXPORT  void  NSMPSSoilLayerEnvironment_Init(NSMPSSoilLayerEnvironment* sle);

EXTERN EXPORT  NSMPSSoilLayerEnvironment*  NSMPSSoilLayerEnvironment_Delete(NSMPSSoilLayerEnvironment* p);

EXTERN EXPORT  NSMPSSoilLayerEnvironment*  NSMPSSoilLayerEnvironment_Clone(NSMPSSoilLayerEnvironment* sle);
#endif
