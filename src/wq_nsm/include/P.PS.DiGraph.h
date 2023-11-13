#pragma once

#include <iostream>

#include <fstream>

#include <string>

#include <map>

#include <vector>


#include "P.PS.Descriptor.CarbonCycle.h"

#include "P.PS.Descriptor.NitrogenCycle.h"

#include "P.PS.Descriptor.PhosphorusCycle.h"

#include "P.PS.Descriptor.Plant.h"

#include "P.PS.Descriptor.Pool.h"

#include "P.PS.Descriptor.SoilLayer.h"

#include "P.PS.Descriptor.SoilSchema.h"

#include "P.PS.Descriptor.SoilParticle.h"

#include "P.PS.Descriptor.SoilParticleDistribution.h"

#include "P.PS.Descriptor.MicroClimate.h"

#include "P.PS.Output.h"



namespace NSMPS
{

	struct Phase;

	struct Pool;

	struct Pathway;



	struct Phase
	{
		double dissolved;

		double sorbed;
	};



	struct  Pool
	{
		std::string name;

		double mass;			// Actual mass resulting from all pathways with limits

		double potentialMass;	// Potential mass resulting from all pathways - no limits

		double k[5];

		std::vector<Pathway*> sinkPathways;		// Set of all "sourcing" pathways

		std::vector<Pathway*> sourcePathways;	// Set of all "sinking" pathways



		public:
		Pool(std::string name);

		void reset();


		void calculatePotentialMass(double timeStep);

		void updatePotentialMass();

		bool redistributeMass();

		bool redistributeSinkPathways();

		void updateActualMass();
	};






	struct PhasePool : Pool
	{
		Phase phase;

		PhasePool(std::string name);
	};






	struct SoilLayer;

	struct ProcessContext
	{
		MicroClimateDescriptor*  microClimateDescriptor;

		NitrogenCycleDescriptor*  nitrogenCycleDescriptor;

		PhosphorusCycleDescriptor*  phosphorusCycleDescriptor;


		long rkLevel;				// Runge-Kutta evaluation level

		SoilLayer* soilLayer;		// Pointer to enclosing soil layer

		ProcessContext(SoilLayer* soilLayer);
	};










	struct  Pathway
	{
		std::string name;

		double massIncrement;

		double massAccumulated;

		double potentialMassQuanta;

		double k[5];

		Pool* donorPool;

		Pool* acceptorPool;


		public:
		Pathway(std::string, Pool* donor, Pool* acceptor);

		void reset();

		virtual void dmdt(ProcessContext* processContext);

		void calculatePotentialMassQuanta(double timeStep);

		void updateActualMassQuanta();
	};



}

