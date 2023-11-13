#pragma once

#include <iterator>
#include <string>
#include <map>
#include <vector>


namespace NSMPS
{

	namespace PLANT_DATABASE
	{

		struct TemperatureRecord
		{
			double baseTemperature;

			double optimalTemperature;

			TemperatureRecord(double baseTemperature, double optimalTemperature);
		};




		struct BiomassProductionRecord
		{
			double RUE;

			double deltaRUE_dcl;

			double RUE_hi;

			double CO2_hi;

			BiomassProductionRecord(double RUE,
									double deltaRUE_dcl,
									double RUE_hi,
									double CO2_hi);
		};




		struct LeafAreaIndexRecord
		{
			double  LAI_max;

			double  fr_PHU_1;

			double  fr_LAI_1;

			double  fr_PHU_2;

			double  fr_LAI_2;

			double fr_PHU_sen;


			LeafAreaIndexRecord(double LAI_max,
								double  fr_PHU_1, 
								double  fr_LAI_1,
								double  fr_PHU_2,
								double  fr_LAI_2,
								double  fr_PHU_sen);
		};


		struct CanopyRootRecord
		{
			double maxStomaticalConductance;

			double maxCanopyHeight;

			double maxRootDepth;

			double min_USLE_C_factor;

			CanopyRootRecord(double maxStomaticalConductance,
							 double maxCanopyHeight,
							 double maxRootDepth,
							 double min_USLE_C_factor);
		};



		struct NutrientParameterRecord
		{
			double fr_N1;	// Normal fraction of N in plant biomass @ emergence

			double fr_N2;	// Normal fraction of N in plant biomass @ 50% maturity

			double fr_N3;	// Normal fraction of N in plant biomass @ maturity


			double fr_P1;	// Normal fraction of P in plant biomass @ emergence

			double fr_P2;	// Normal fraction of P in plant biomass @ 50% maturity

			double fr_P3;	// Normal fraction of P in plant biomass @ maturity

			NutrientParameterRecord(double fr_N1,
									double fr_N2,	
									double fr_N3,	
									double fr_P1,	
									double fr_P2,	
									double fr_P3);

		};



		class Database
		{
			private:

			static Database* pInstance;

			std::map<std::string, TemperatureRecord*> temperatureRecords;

			std::map<std::string, BiomassProductionRecord*> biomassProductionRecords;

			std::map<std::string, LeafAreaIndexRecord*> leafAreaIndexRecords;

			std::map<std::string, CanopyRootRecord*> canopyRootRecords;

			std::map<std::string, NutrientParameterRecord*> nutrientParameterRecords;


			public:

			static Database* getInstance();

			TemperatureRecord* getTemperatureRecord(std::string name);

			BiomassProductionRecord* getBiomassProductionRecord(std::string name);

			LeafAreaIndexRecord* getLeafAreaIndexRecord(std::string name);

			CanopyRootRecord* getCanopyRootRecord(std::string name);

			NutrientParameterRecord* getNutrientParameterRecord(std::string name);



			protected:

			Database();


			private:

			void buildTemperatureRecords();

			void buildBiomassProductionRecords();

			void buildLeafAreaIndexRecords();

			void buildCanopyRootRecords();

			void buildNutrientParameterRecords();
		};

	}

}
