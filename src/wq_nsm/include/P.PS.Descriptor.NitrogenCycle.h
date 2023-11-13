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

#include "P.PS.Output.h"





namespace NSMPS
{



	struct NitrogenCycleDescriptor
	{

		// F I E L D S 

		PoolDescriptor organicFresh;

		double fractionHumicActive;		// Fraction of active pool in humic pool

		PoolDescriptor organicActive;

		PoolDescriptor organicStable;

		PoolDescriptor mineralNH4;

		PoolDescriptor mineralNO3;


		// C O N S T R U C T O R S

		NitrogenCycleDescriptor();

		// Can initialize from a user-supplied "NSMPSNitrogenEnvironment"

		NitrogenCycleDescriptor(NSMPSNitrogenEnvironment* ne);

		// Can initialize from a CarbonCycleDescriptor also.

		NitrogenCycleDescriptor(double depth, CarbonCycleDescriptor* ccd);

		NitrogenCycleDescriptor(SoilSchemaDescriptor* ssd);

		//NitrogenCycleDescriptor(NitrogenCycleDescriptor&);


		// M E T H O D S

		NitrogenCycleDescriptor& operator+=(const NitrogenCycleDescriptor&);


		NitrogenCycleDescriptor& operator<<(const NSMPSNitrogenEnvironment&);

		NitrogenCycleDescriptor& operator>>(NSMPSNitrogenEnvironment&);


        void  PullSpecieFrom(const NSM_SPECIES_ENUM specie, const NSMPSNitrogenEnvironment*);

        void  PushSpecieTo(NSM_SPECIES_ENUM specie, NSMPSNitrogenEnvironment*);


		void zero();

		void print(double time, OutputFile* of);

		void calculatePartitioning(SoilLayerDescriptor*);
	};




}
