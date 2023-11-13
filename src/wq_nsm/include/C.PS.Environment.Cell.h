#pragma once
#ifndef C_PS_ENVIRONMENT_CELL_H
#define C_PS_ENVIRONMENT_CELL_H

#include "LINKAGE.h"

#include "C.PS.Environment.Plant.h"

#include "C.PS.CSchema.Flora.h"

#include "C.PS.Environment.Particle.h"

#include "C.PS.Environment.SoilLayer.h"

#include "C.PS.CSchema.Soil.h"









struct NSMPSCellEnvironment
{
	// Plant environment
	unsigned long  numberPlants;

	NSMPSPlantEnvironment**  plantEnvironments;


	// Soil environment
	unsigned long  numberSoilLayers;

	NSMPSSoilLayerEnvironment**  soilLayerEnvironments;
};

typedef struct NSMPSCellEnvironment NSMPSCellEnvironment;










// C E L L   E N V I R O N M E N T   P R O T O T Y P E S 


EXTERN EXPORT  NSMPSCellEnvironment*  NSMPSCellEnvironment_Create(long numberPlants, long numberSoilLayers);


EXTERN EXPORT  NSMPSCellEnvironment*  NSMPSCellEnvironment_CreateViaSchemas(NSMPSCSchemaFlora* fse,
								NSMPSCSchemaSoil* sse,  NSMSoilParticleDistributionEnvironment* spe);


EXTERN EXPORT  NSMPSCellEnvironment*  NSMPSCellEnvironment_Delete(NSMPSCellEnvironment* ce);


EXTERN EXPORT  NSMPSCellEnvironment*  NSMPSCellEnvironment_Copy(NSMPSCellEnvironment* cea);
#endif
