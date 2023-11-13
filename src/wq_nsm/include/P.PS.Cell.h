#pragma once

#include <map>

#include <string>

#include <vector>

#include "C.NutrientProperties.h"

#include "C.PS.Environment.Nutrients.h"

#include "P.PS.Flora.h"

#include "P.PS.Soil.h"

#include "P.PS.DiGraph.h"

#include "P.PS.NitrogenDiGraph.h"

#include "P.PS.PhosphorusDiGraph.h"






// D E C L A R A T I O N S 

namespace NSMPS
{
	class Cell;
}





// D E F I N I T I O N S 

namespace NSMPS
{

	class  Cell
		{
		public: CellDescriptor cellDescriptor;



		public:

		// F L O R A   C O M P O N E N T S

		Flora* flora;



		// S O I L   C O M P O N E N T S

		std::vector<SoilLayerDescriptor*> soilLayerDescriptors;

		std::vector<SoilLayer*> soilLayers;




		// N U T R I E N T   C O M P O N E N T S

		std::vector<NSMPS::ProcessContext*> processContexts;


		std::vector<NSMPS::NitrogenCycleDescriptor*> nitrogenCycleDescriptors;

		std::vector<NSMPS::NITROGEN::NitrogenDiGraph*> nitrogenDiGraphs;


		std::vector<NSMPS::PhosphorusCycleDescriptor*> phosphorusCycleDescriptors;

		std::vector<NSMPS::PHOSPHORUS::PhosphorusDiGraph*> phosphorusDiGraphs;


		NSMPS::OutputFile* logFile;



		// S E L F 

		Cell(CellDescriptor* cd);

		void advance(NSMPS::MicroClimateDescriptor* mcd);

		void writeToLog(NSMPS::MicroClimateDescriptor* mcd);



		// U T I L I T Y   F U N C T I O N S

		double calculateFloraWaterDemand(double timeStep);

		double registerFloraNitrogenPotentialDemandRate();

		double retrieveFloraNitrogenActualDemandMass();

		void LoadNutrientProperties(NSMNutrientProperties* np);



        // T R A N S F E R

        void PullSpeciesFrom(const NSMPSCellEnvironment*);

        void PullSpeciesFrom(const NSMPSCellEnvironment* ce, const int layer);
       
        void PullSpecieFrom(const NSMPSCellEnvironment* ce, const NSM_SPECIES_ENUM specie, const int layer);


        void PushSpeciesTo(NSMPSCellEnvironment*);

        void PushSpeciesTo(NSMPSCellEnvironment* ce, const int layer);

        void PushSpecieTo(NSMPSCellEnvironment* ce, NSM_SPECIES_ENUM specie, const int layer);

		};

}


