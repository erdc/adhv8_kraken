#pragma once

#include <iostream>

#include <fstream>

#include <sstream>

#include <map>

#include <string>

#include <vector>


#include "C.Enumerations.h"

#include "C.PS.Environment.Cell.h"

#include "C.PS.Environment.Plant.h"

#include "C.PS.CSchema.Flora.h"

#include "C.PS.Environment.Nutrients.h"

#include "C.PS.Environment.Particle.h"

#include "C.PS.Environment.SoilLayer.h"

#include "C.PS.CSchema.Soil.h"


#include "P.PS.Descriptor.SoilLayer.h"

#include "P.PS.Descriptor.SoilParticle.h"

#include "P.PS.Descriptor.SoilParticleDistribution.h"




namespace NSMPS
{


	struct PoolDescriptor
	{
		double mass;

		double KOC;

		double fractionDissolved;

		double fractionSolid;

		double fsum;		// Validity check

        NSM_PARTITIONABLE_ENUM partitionAttr;

		PoolDescriptor();

		PoolDescriptor(PoolDescriptor&);

		void calculatePartitioning(SoilLayerDescriptor*);
	};



}