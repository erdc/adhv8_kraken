#pragma once

#include <iostream>

#include <map>

#include <string>

#include <vector>

#include "P.PS.DiGraph.h"

#include "P.PS.Soil.h"




namespace NSMPS
{


	namespace NITROGEN
	{



		struct ZeroOrderPathway : Pathway
		{
			private:double demandRate;		// Kg/day

			public:ZeroOrderPathway(std::string name, Pool* donor, Pool* acceptor);

			public:double zeroOrderKinetics();

			public:void setDemandRate(double demandRate);

			public:void dmdt(ProcessContext* processContext);
		};





		struct FirstOrderPathway : Pathway
		{
			double beta;

			FirstOrderPathway(std::string name, Pool* donor, Pool* acceptor, double beta=0.0);

			double firstOrderKinetics(double mass);

			void dmdt(ProcessContext* processContext);
		};





		struct DemandPathway : Pathway
		{
			double demandRate;

			DemandPathway(std::string name, Pool* donor, Pool* acceptor);

			void dmdt(ProcessContext* processContext);
		};





		struct MineralizationPathway : FirstOrderPathway
		{
			MineralizationPathway(std::string name, Pool* donor, Pool* acceptor);

			void dmdt(ProcessContext* processContext);
		};





		struct ResidueMineralizationPathway : FirstOrderPathway
		{
			ResidueMineralizationPathway(std::string name, Pool* donor, Pool* acceptor);

			void dmdt(ProcessContext* processContext);
		};





		struct DenitrificationPathway : FirstOrderPathway
		{
			DenitrificationPathway(std::string name, Pool* donor, Pool* acceptor);

			void dmdt(ProcessContext* processContext);
		};




		
		struct VolatilizationPathway : FirstOrderPathway
		{
			VolatilizationPathway(std::string name, Pool* donor, Pool* acceptor);

			void dmdt(ProcessContext* processContext);
		};




		
		struct NitrificationPathway : FirstOrderPathway
		{
			NitrificationPathway(std::string name, Pool* donor, Pool* acceptor);

			void dmdt(ProcessContext* processContext);
		};





		struct RainfallPathway : Pathway
		{
			RainfallPathway(std::string name, Pool* donor, Pool* acceptor);

			void dmdt(ProcessContext* processContext);
		};





		struct NitrogenDiGraph
		{

			// P O O L    S P E C I F I C A T I O N S

			// Organic Pools

			Pool*  organicFreshPool;

			Pool*  organicActivePool;

			Pool*  organicStablePool;



			// Mineral Pools

			Pool*  mineralNH4Pool;

			Pool*  mineralNO3Pool;



			// Applied Pools

			Pool*  appliedOrganicPool;

			Pool*  appliedNH4Pool;

			Pool*  appliedNO3Pool;



			// Atmospheric Pools

			Pool*  denitrifiedPool;

			Pool*  volatilizedPool;

			Pool*  rainfallNO3Pool;



			// Flora Pools

			Pool*  floraNH4Pool;

			Pool*  floraNO3Pool;


			std::vector<Pool*> pools;		// The set of all pools



			// P A T H W A Y    S P E C I F I C A T I O N S


			// Organic Pathways

			FirstOrderPathway*  decayPathway_OrganicFresh_OrganicActive;

			FirstOrderPathway*  decayPathway_OrganicFresh_OrganicStable;

			FirstOrderPathway*  decayPathway_OrganicStable_OrganicActive;


			// Organic-to-Mineral Pathways

			MineralizationPathway*  mineralizationPathway_OrganicActive_MineralNO3;

			ResidueMineralizationPathway*  residueMineralizationPathway_OrganicFresh_MineralNO3;


			// Mineral Pathways

			NitrificationPathway*  NitrificationPathway_MineralNH4_MineralNO3;


			// Air Pathways

			DenitrificationPathway*  denitrificationPathway_MineralNO3_Denitrified;

			VolatilizationPathway*  volatilizationPathway_MineralNH4_Volatilized;

			RainfallPathway*  rainfallPathway_RainfallNO3_Mineral_NO3;


			// Applied Pathways

			FirstOrderPathway*  sourcePathway_OrganicApplied_OrganicActive;

			FirstOrderPathway*  sourcePathway_AppliedNH4_MineralNH4;

			FirstOrderPathway*  sourcePathway_AppliedNO3_MineralNO3;


			// Mineral-to-Flora Pathways

			ZeroOrderPathway*  uptakePathway_MineralNH4_Flora;

			ZeroOrderPathway*  uptakePathway_MineralNO3_Flora;
		


			std::vector<Pathway*>  pathways;	// The set of all pathways


			// S U P P O R T   F I E L D S

			bool diGraphInSurfaceLayer;			// True if digraph in surface soil layer


			// DiGraph state: Internal use only

			bool enableOutput;					// Set to true if you want output

			OutputFile*  outputFile;	// For debugging graph




			// M E T H O D S

			NitrogenDiGraph(SoilLayer* soilLayer);

			public:bool  hasMass();

			public:void  setFloraNitrogenDemand(double demandRate);

			private:void  setFloraNH4PathwayDemand(double demandRate);	// OBSOLETE??

			private:void  setFloraNO3PathwayDemand(double demandRate);	// OBSOLETE??

			public:double  getFloraNitrogenMassUptake();


            public:
			void PushSpeciesTo(NitrogenCycleDescriptor*);

			void PushSpecieTo(const NSM_SPECIES_ENUM specie, NitrogenCycleDescriptor*);

		    NitrogenDiGraph& operator<<(const NitrogenCycleDescriptor&);


			void PullSpeciesFrom(const NitrogenCycleDescriptor*);

			void PullSpecieFrom(const NSM_SPECIES_ENUM specie, const NitrogenCycleDescriptor*);

		    NitrogenDiGraph& operator>>(NitrogenCycleDescriptor&);


			public:
			void  advance(ProcessContext*);

			void  printPools();
		};


	}

}

