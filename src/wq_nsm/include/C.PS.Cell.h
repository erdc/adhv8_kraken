#pragma once
#ifndef C_PS_CELL_H
#define C_PS_CELL_H

#include "LINKAGE.h"

#include "C.Enumerations.h"

#include "C.METEnvironment.h"

#include "C.NutrientProperties.h"

#include "C.PS.Environment.Plant.h"

#include "C.PS.CSchema.Flora.h"

#include "C.PS.Environment.Nutrients.h"

#include "C.PS.Environment.Particle.h"

#include "C.PS.Environment.SoilLayer.h"

#include "C.PS.CSchema.Soil.h"




struct NSMCPSCell
{
	NSMPSCellEnvironment*  cellEnvironment;	// Each cell has it's own NSMPSCellEnvironment

	struct NSMCPSCell*  cell;
};

typedef struct NSMCPSCell NSMCPSCell;







//// NSMCPSCELL    P R O T O T Y P E S 


EXTERN EXPORT  NSMCPSCell*  NSMCPSCell_Create(NSMPSCellEnvironment* cellEnv);


EXTERN EXPORT  NSMCPSCell*  NSMCPSCell_Delete(NSMCPSCell*);


EXTERN EXPORT  NSMPSCellEnvironment*  NSMCPSCell_Advance(NSMMetEnvironment* me, NSMCPSCell* p);



EXTERN EXPORT  void  NSMCPSCell_Update(NSMCPSCell* p);


EXTERN EXPORT  void  NSMCPSCell_UpdateAllSpecies(NSMCPSCell* p);


EXTERN EXPORT  void  NSMCPSCell_UpdateAllSpeciesByLayer(NSMCPSCell* p, int layer);


EXTERN EXPORT  void  NSMCPSCell_UpdateSpecieByLayer(NSMCPSCell* p, enum NSM_SPECIES_ENUM specie, int layer);



EXTERN EXPORT  unsigned long  NSMCPSCell_GetNumberOfSoilLayers(NSMCPSCell* p);


EXTERN EXPORT  unsigned long  NSMCPSCell_GetSoilWaterDemand(NSMCPSCell* p, double* buffer);


EXTERN EXPORT  double NSMCPSCell_GetSoilNutrientMass(NSMCPSCell* pCell, enum NSM_SPECIES_ENUM id, unsigned long layer);


EXTERN EXPORT  unsigned long  NSMCPSCell_GetSoilNutrientConc(NSMCPSCell* pCell, enum NSM_SPECIES_ENUM id, double* buffer);


EXTERN EXPORT  double  NSMCPSCell_ApplyMassFluxToSoilSurface(NSMCPSCell* pCell, enum NSM_SPECIES_ENUM id, double flux, double timeStep);


EXTERN EXPORT  void  NSMCPSCell_LoadNutrientProperties(NSMCPSCell* p, NSMNutrientProperties* np);
#endif
