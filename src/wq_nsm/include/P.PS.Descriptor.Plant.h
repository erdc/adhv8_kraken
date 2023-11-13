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



	struct PlantDescriptor
	{
		std::string name;

		double coverFraction;

		PlantDescriptor();

		PlantDescriptor(NSMPSPlantEnvironment*);
	};


}


