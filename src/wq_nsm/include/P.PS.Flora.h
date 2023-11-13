#pragma once

#include <iostream>

#include <string>

#include <vector>


#include "C.PS.CSchema.Flora.h"


#include "P.PS.Descriptor.Plant.h"

#include "P.PS.Descriptor.SoilLayer.h"

#include "P.PS.Descriptor.Cell.h"

#include "P.PS.Descriptor.MicroClimate.h"


#include "P.PS.Plant.h"



namespace NSMPS
{



	struct Flora
	{
		std::vector<PlantDescriptor*> plantDescriptors;

		std::vector<Plant*> plants;


		public:
		Flora(CellDescriptor* cd);

		void  advance(const MicroClimateDescriptor* mcd);

		void distributeNitrogenUptakeAmongPlants(double totalNitrogenUptakeMass);

		void calculatePlantGrowthConstraints(const MicroClimateDescriptor* mcd);




		// UTILITY METHODS 

		double  getTotalBiomass();

		double  calculatePotentialWaterDemandRateByLayer(SoilLayerDescriptor* sld);

		double  calculatePotentialNitrogenDemandRateByLayer(SoilLayerDescriptor* sld);

		double  getTotalNitrogenMass();

	};


}