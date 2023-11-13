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

#include "P.PS.Descriptor.NitrogenCycle.h"

#include "P.PS.Output.h"





namespace NSMPS
{

	





	struct PhosphorusCycleDescriptor
	{
		// F I E L D S 

		PoolDescriptor organicFresh;

		PoolDescriptor organicActive;

		PoolDescriptor organicStable;

		PoolDescriptor mineralActive;

		PoolDescriptor mineralStable;

		PoolDescriptor mineralSolution;


		// C O N S T R U C T O R S

		PhosphorusCycleDescriptor();

		// Can initialize from a user-supplied "NSMPSNitrogenEnvironment"

		PhosphorusCycleDescriptor(NSMPSPhosphorusEnvironment* pe);


		// M E T H O D S

		PhosphorusCycleDescriptor& operator+=(const PhosphorusCycleDescriptor&);


		PhosphorusCycleDescriptor& operator<<(const NSMPSPhosphorusEnvironment&);
		
		PhosphorusCycleDescriptor& operator>>(NSMPSPhosphorusEnvironment&);


        void  PullSpecieFrom(const NSM_SPECIES_ENUM specie, const NSMPSPhosphorusEnvironment*);

        void  PushSpecieTo(const NSM_SPECIES_ENUM specie, NSMPSPhosphorusEnvironment*);


		void zero();

		void print(double time, OutputFile* of);

		void calculatePartitioning(SoilLayerDescriptor*);
	};





}