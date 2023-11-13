#pragma once

#include <iostream>

#include <fstream>


#include <string>

#include <vector>

#include "C.PS.CSchema.Flora.h"


#include "P.PS.Descriptor.Plant.h"

#include "P.PS.Descriptor.SoilLayer.h"

#include "P.PS.Descriptor.MicroClimate.h"





namespace NSMPS
{

	struct Plant;


	namespace FLORA
	{


		struct CanopyZone
		{
			double height;
			double maxCanopyHeight;

			double LAI;
			double previousLAI;

			double fr_PHU_sen;
			double LAI_max;

			double maximumStomatalConductance;

			double LAIShapeCoeff1;
			double LAIShapeCoeff2;

			double fractionMaximumLeafAreaIndex;
			double previousFractionMaximumLeafAreaIndex;

			CanopyZone(PlantDescriptor*);

			void calculateLeafArea(Plant*);
		};





		struct RootZone
		{
			double rootDepth;
			double dailyRootBiomassFraction;

			double maxRootDepth;
			double maxStomaticalConductance;
			double min_USLE_C_factor;

			RootZone(PlantDescriptor*);

			void calculateRootDepth(Plant*);
		};






		struct WaterDemand
		{
			double transpirationRate;

			WaterDemand();

			void calculateTranspirationRate(Plant* plant, const MicroClimateDescriptor* mcd);

			double calculatePotentialDemandRateByLayer(Plant* plant, SoilLayerDescriptor* sld);
		};






		struct PotentialGrowth
		{
			double daily_biomass_delta;
			double potentialDailyBiomassIncrease;
		};






		struct Nutrient
		{
			// Normal fraction of nutrient in plant biomass @ emergence
			double  normalFractionInPlantBiomassAtEmergence;

			// Normal fraction of nutrient in plant biomass @ maturity
			double  normalFractionInPlantBiomassAtMaturity;

			double  uptakeShapeCoefficient1;		// Unique per plant

			double  uptakeShapeCoefficient2;

			double  stressFactor;

			double  fractionOfPlantBiomass;

			
			double  demandMass;					// Nutrient demand mass this time-step

			double  demandRate;					// Nutrient demand rate this time-step

			double  uptakeMass;					// Actual nutrient uptake mass this time-step

			double  totalMass;					// Total mass of nutrient in plant @ time=n+1

			double  previousTotalMass;			// Total mass of nutrient in plant @ time=n


			Nutrient();

			virtual void  calculatePotentialDemandRate(Plant* plant, double timeStep);

			virtual double  calculatePotentialDemandRateByLayer(Plant* plant, SoilLayerDescriptor* sld);

			virtual double  calculateStress();
		};








		struct Nitrogen : Nutrient
		{
			Nitrogen(PlantDescriptor*);
		};






		struct Phosphorus : Nutrient
		{
			Phosphorus(PlantDescriptor*);
		};





		struct GrowthParameters
		{
			// Contains plant specific growth parameters

			double baseTemperature;				// Plant specific base temperature. Units:C

			double optimalTemperature;			// Plant specific optimal temperature. Units:C



			struct BiomassProduction
			{
				double RUE;				// Radiation use efficiency

				double deltaRUE_dcl;

				double RUE_hi;

				double CO2_hi;

			} biomassProduction;



			public:
			GrowthParameters(PlantDescriptor*);
		};




	}



	struct  Plant
	{
		PlantDescriptor*  plantDescriptor;		// Every plant has one.


		// Energy parameters

		double  heatUnits;

		double  potentialHeatUnits;

		double  fractionPotentialHeatUnits;



		// Growth parameters

		double  H_phosyn;			// Solar radiation intercepted by leaf area per timestep.
									// Note: timestep is in seconds. Units: MJ/m^2.
									// ENERGY NOT POWER.

		double  biomass;			// Plant biomass. Units:Kg.

		double  deltaBiomass;		// Change in plant biomass. Units:Kg.

		double  optimalBiomass;		// Optimal biomass. Units:Kg. (growth with no limiters)

		double  plantGrowthConstraintFactor;	// Fraction of potential growth acheivable



		// Composites

		FLORA::CanopyZone*  canopyZone;

		FLORA::RootZone*  rootZone;

		FLORA::GrowthParameters*  growthParameters;




		// Resources

		FLORA::WaterDemand*  waterDemand;

		FLORA::Nitrogen*  nitrogen;




		public:
		Plant(PlantDescriptor* pd);



		// Growth functions

		void  calculateHU(const MicroClimateDescriptor* mcd);

		void  calculateSolarRadiationReceived(const MicroClimateDescriptor* mcd);

		double  calculatePlantGrowthConstraintFactor(const MicroClimateDescriptor* mcd);

		void  grow(const MicroClimateDescriptor* mcd);

		double calculateTemperatureStress(const MicroClimateDescriptor* mcd);


	};

	
	

}