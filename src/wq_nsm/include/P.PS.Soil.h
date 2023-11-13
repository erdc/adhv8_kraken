#pragma once

#include <vector>

#include "P.PS.Descriptor.SoilLayer.h"



namespace NSMPS
{

	struct SoilLayerComputedProperties
	{
		public:
		double  nutrientCyclingWaterFactor;

		double  nutrientCyclingTemperatureFactor;

		double  nutrientCyclingFactor;

		double  fractionOrganicCarbon;

		public:
		SoilLayerComputedProperties();

		void  update(SoilLayerDescriptor*  sld);
	};




	struct SoilLayer
	{
		SoilLayerDescriptor*  soilLayerDescriptor;

		SoilLayerComputedProperties  soilLayerComputedProperties;

		SoilLayer(SoilLayerDescriptor*  sld);

		void  update();
	};

}
