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


#include "P.PS.Descriptor.SoilLayer.h"


namespace NSMPS
{

	// Forward declarations

	struct  PlantDescriptor;				

	struct  SoilParticleDescriptor;

	struct  PoolDescriptor;

	struct  CarbonCycleDescriptor;

	struct  NitrogenCycleDescriptor;

	struct  PhosphorusCycleDescriptor;

	struct  SoilLayerDescriptor;	



	struct SoilSchemaDescriptor
	{
		SoilLayerDescriptor*  soilLayerDescriptor;

		CarbonCycleDescriptor*  carbonCycleDescriptor;

		NitrogenCycleDescriptor*  nitrogenCycleDescriptor;

		PhosphorusCycleDescriptor*  phosphorusCycleDescriptor;

		SoilSchemaDescriptor();

		SoilSchemaDescriptor(SoilSchemaDescriptor&);
	};

}
