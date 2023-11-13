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


#include "P.PS.Descriptor.Plant.h"

#include "P.PS.Descriptor.Pool.h"

#include "P.PS.Descriptor.SoilLayer.h"

#include "P.PS.Descriptor.SoilSchema.h"

#include "P.PS.Descriptor.SoilParticle.h"

#include "P.PS.Descriptor.SoilParticleDistribution.h"







namespace NSMPS
{



	struct CarbonCycleDescriptor
	{
		// F I E L D S 

		PoolDescriptor organicActive;

		PoolDescriptor organicStable;
		
		double fractionCarbonInSoil;	// Soil Organic Matter fraction



		// C O N S T R U C T O R S

		CarbonCycleDescriptor();

		CarbonCycleDescriptor(SoilLayerDescriptor*);

		CarbonCycleDescriptor(SoilSchemaDescriptor*);

		CarbonCycleDescriptor(NSMPSCarbonEnvironment*);



		// M E T H O D S

		CarbonCycleDescriptor& operator+=(const CarbonCycleDescriptor&);

		CarbonCycleDescriptor& operator<<(const NSMPSCarbonEnvironment&);

		CarbonCycleDescriptor& operator>>(NSMPSCarbonEnvironment&);
	};



}
