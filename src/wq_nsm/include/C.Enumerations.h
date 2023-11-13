#pragma once
#ifndef C_ENUMERATIONS_H
#define C_ENUMERATIONS_H

#include "LINKAGE.h"
#ifndef NSM_SEDIMENT_SORBING_ENUM
 enum  NSM_SEDIMENT_SORBING_ENUM
{
    NSM_SEDIMENT_IS_SORBING_Enum,
    NSM_SEDIMENT_IS_NOT_SORBING_Enum
}; 
#endif


 enum  NSM_PARTITIONABLE_ENUM
{
    NSM_IS_PARTITIONABLE_Enum,
    NSM_IS_FULLY_DISSOLVED_Enum,
    NSM_IS_FULLY_SORBED_Enum
}; 



enum  NSM_SPECIES_ENUM
{
	NSM_NO2_Enum = 0,
	NSM_NO3_Enum,
	NSM_NH4_Enum,
	NSM_ORGANIC_NITROGEN_Enum,
	NSM_ORGANIC_PHOSPHORUS_Enum,
	NSM_DISSOLVED_PHOSPHORUS_Enum,
	NSM_ALGAE_Enum,
	NSM_CBOD_Enum,
	NSM_DISSOLVED_OXYGEN_Enum,
	NSM_Soil_N_organicActive_Enum = 100,
	NSM_Soil_N_organicStable_Enum,
	NSM_Soil_P_organicActive_Enum,
	NSM_Soil_P_organicStable_Enum,
	NSM_Soil_P_mineralActive_Enum,
	NSM_Soil_P_mineralStable_Enum
}; 




//  F U N C T I O N   P R O T O T Y P E S

EXPORT  int NSM_SPECIES_ENUM_IsChannelSpecie(enum NSM_SPECIES_ENUM specieEnum);

EXPORT  int NSM_SPECIES_ENUM_IsOverlandSpecie(enum NSM_SPECIES_ENUM specieEnum);

EXPORT  int NSM_SPECIES_ENUM_IsSoilSpecie(enum NSM_SPECIES_ENUM specieEnum);

EXPORT  enum NSM_PARTITIONABLE_ENUM  NSM_SPECIES_ENUM_IsPartitionable(enum NSM_SPECIES_ENUM specieEnum);
#endif
