#pragma once

#include <iostream>

#include <map>

#include <string>

#include <vector>


#include "P.PS.DiGraph.h"

#include "P.PS.NitrogenDiGraph.h"	// Some Phosphorus pathways are functions of nitrogen conc.

#include "P.PS.Soil.h"




namespace NSMPS
{


	namespace PHOSPHORUS
	{


	

		struct ZeroOrderPathway : Pathway
		{
			private:double  demandRate;		// Kg/day

			public:ZeroOrderPathway(std::string name, Pool* donor, Pool* acceptor);

			public:double  zeroOrderKinetics();

			public:void  setDemandRate(double demandRate);

			public:void  dmdt(ProcessContext* processContext);
		};





		struct FirstOrderPathway : Pathway
		{
			double  beta;

			FirstOrderPathway(std::string name, Pool* donor, Pool* acceptor, double beta=0.0);

			double  firstOrderKinetics(double mass);

			void  dmdt(ProcessContext* processContext);
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





		struct PhosphorusDiGraph
		{

			// P O O L    S P E C I F I C A T I O N S

			// Organic Pools

			Pool*  organicFreshPool;

			Pool*  organicActivePool;

			Pool*  organicStablePool;



			// Mineral Pools

			Pool*  mineralActivePool;

			Pool*  mineralStablePool;

			Pool*  mineralSolutionPool;



			// Applied Pools

			Pool*  appliedOrganicPool;

			Pool*  appliedMineralPool;



			// Flora Pools

			Pool*  floraMineralPool;



			std::vector<Pool*> pools;		// The set of all pools



			// P A T H W A Y    S P E C I F I C A T I O N S


			// Organic Pathways

			FirstOrderPathway*  decayPathway_OrganicFresh_OrganicActive;

			FirstOrderPathway*  decayPathway_OrganicFresh_OrganicStable;

			ZeroOrderPathway*  equilibriumPathway_OrganicStable_OrganicActive;

			ZeroOrderPathway*  equilibriumPathway_OrganicActive_OrganicStable;


			// Organic-to-Mineral Pathways

			MineralizationPathway*  mineralizationPathway_OrganicActive_MineralSolution;

			ResidueMineralizationPathway*  residueMineralizationPathway_OrganicFresh_MineralSolution;


			// Mineral Pathways

			FirstOrderPathway*  decayPathway_MineralStable_MineralActive;

			FirstOrderPathway*  decayPathway_MineralActive_MineralSolution;


			// Applied Pathways

			FirstOrderPathway*  sourcePathway_OrganicApplied_OrganicActive;

			FirstOrderPathway*  sourcePathway_MineralApplied_MineralSolution;


			// Mineral-to-Flora Pathways

			ZeroOrderPathway*  uptakePathway_MineralSolution_Flora;



			std::vector<Pathway*>  pathways;		// The set of all pathways




			// S U P P O R T   F I E L D S

			bool diGraphInSurfaceLayer;				// True if digraph associated 
													// with surface soil layer


			// DiGraph state: Internal use only

			bool enableOutput;						// Set to true if you want output

			OutputFile*  outputFile;		// For examining di-graph state





			// M E T H O D S

			PhosphorusDiGraph(SoilLayer* soilLayer);


			void calculateHumicEqilibriumTransferRate(ProcessContext* processContext);


			public:void setFloraPhosphorusDemand(double demandRate);

			public:double getFloraPhosphorusMassUptake();


            public:
			void PushSpeciesTo(PhosphorusCycleDescriptor*);

			void PushSpecieTo(const NSM_SPECIES_ENUM specie, PhosphorusCycleDescriptor*);

            PhosphorusDiGraph& operator>>(PhosphorusCycleDescriptor&);


			void PullSpeciesFrom(const PhosphorusCycleDescriptor*);

			void PullSpecieFrom(const NSM_SPECIES_ENUM specie, const PhosphorusCycleDescriptor*);

		    PhosphorusDiGraph& operator<<(const PhosphorusCycleDescriptor&);




			public:
			void advance(ProcessContext*);

			void printPools();
		};


	}

}



