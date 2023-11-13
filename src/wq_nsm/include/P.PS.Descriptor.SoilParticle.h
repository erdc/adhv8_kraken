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

	struct SoilParticleDescriptor
	{
		public:double fraction;		// Fraction of particle in soil

		public:double FOC;			// Fraction organic carbon -- not percentage

		public:SoilParticleDescriptor();

		public:SoilParticleDescriptor(NSMSoilParticleEnvironment*);
	};

}