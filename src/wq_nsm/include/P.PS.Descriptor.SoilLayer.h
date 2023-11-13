#pragma once

#include <iostream>

#include <fstream>

#include <sstream>

#include <map>

#include <string>

#include <vector>


#include "C.PS.Environment.Cell.h"

#include "C.PS.Environment.Plant.h"

#include "C.PS.CSchema.Flora.h"

#include "C.PS.Environment.Nutrients.h"

#include "C.PS.Environment.Particle.h"

#include "C.PS.Environment.SoilLayer.h"

#include "C.PS.CSchema.Soil.h"










namespace NSMPS
{


	struct SoilLayerDescriptor
	{   
		unsigned long layer;	// Layer id.

		double depth;			// Depth @ layer bottom. Units:m

		double thickness;		// Soil layer thickness.  Units:m

		double area;			// Area of soil layer. Units: m^2

		double density;			// Soil layer density. Units: kg/m^3

		double porosity;		// Soil layer porosity. Dimensionless.


		double fractionCarbonInSoil;	// Soil Organic Matter fraction

		double dissolvedOrganicCarbon;	// conc. (Till carbon cycle is full-featured)


		double waterContent;	// Water content of soil layer. Units: m

		double fieldCapacity;	// Field capacity of soil layer. Units: m

		double wiltingPoint;	// Wilting point of soil layer. Units: m

		double volumetricWaterContent;	// WaterContent. Units:dimensionless

		double volumetricFieldCapacity;	// Field capacity of soil layer. Units:dimensionless

		double volumetricWiltingPoint;	// Wilting point of soil layer. Units:dimensionless



		double floraWaterDemand;// Water demanded by flora community. Units: m

		double floraWaterSupply;// Water supplied to flora community. Units: m


		double pH;				// Units: dimensionless

		double soilTemp;		// Soil layer temperature. Units: C

		double residue;			// Mass of plant residue in layer. Units: Kg



		SoilLayerDescriptor(unsigned long layerId=0);

		SoilLayerDescriptor(NSMPSSoilLayerEnvironment*);


        // SCREAM
        double GetSoilWaterDemand();

        void SetSoilWaterSupply(double supply);
	};





}