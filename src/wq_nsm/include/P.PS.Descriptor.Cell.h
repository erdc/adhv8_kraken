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


#include "P.PS.Descriptor.CarbonCycle.h"

#include "P.PS.Descriptor.NitrogenCycle.h"

#include "P.PS.Descriptor.PhosphorusCycle.h"

#include "P.PS.Descriptor.Plant.h"

#include "P.PS.Descriptor.Pool.h"

#include "P.PS.Descriptor.SoilLayer.h"

#include "P.PS.Descriptor.SoilSchema.h"

#include "P.PS.Descriptor.SoilParticle.h"

#include "P.PS.Descriptor.SoilParticleDistribution.h"

#include "P.PS.Output.h"






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

	struct  SoilSchemaDescriptor;	





	struct CellDescriptor
	{

		// F I E L D S

		// List of flora within a cell

		std::vector<PlantDescriptor*>  plantDescriptors;


		// Vertical schema describing soil profile

		std::vector<SoilSchemaDescriptor*>  soilSchemaDescriptors;	




		// C O N S T R U C T O R S

		CellDescriptor();

		CellDescriptor(NSMPSCellEnvironment*);

		CellDescriptor(CellDescriptor&);



		// M E T H O D S
        
		CellDescriptor& operator<<(const NSMPSCellEnvironment&);

		CellDescriptor& operator>>(NSMPSCellEnvironment&);

	};




}