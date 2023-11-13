#pragma once
#ifndef C_PS_CSCHEMA_SOIL_H
#define C_PS_CSCHEMA_SOIL_H

#include "LINKAGE.h"

#include "C.METEnvironment.h"

#include "C.PS.Environment.Cell.h"

#include "C.PS.CSchema.Flora.h"

#include "C.PS.Environment.Nutrients.h"

#include "C.PS.Environment.Particle.h"

#include "C.PS.Environment.SoilLayer.h"

#include "C.Library.h"







struct NSMPSCSchemaSoil
{
	unsigned long numberSoilLayers;

	NSMListNode* soilLayerList;
};

typedef struct NSMPSCSchemaSoil NSMPSCSchemaSoil;





// S O I L  S C H E M A   P R O T O T Y P E S

EXTERN EXPORT  NSMPSCSchemaSoil*  NSMPSCSchemaSoil_Create();

EXTERN EXPORT  void  NSMPSCSchemaSoil_Init(NSMPSCSchemaSoil* env);

EXTERN EXPORT  unsigned long  NSMPSCSchemaSoil_Add(NSMPSCSchemaSoil* sse, NSMPSSoilLayerEnvironment* sle);

EXTERN EXPORT  unsigned long  NSMPSCSchemaSoil_AddList(NSMPSCSchemaSoil* sc, int numberSoilLayers, NSMPSSoilLayerEnvironment** soilLayerEnvs);

EXTERN EXPORT  unsigned long  NSMPSCSchemaSoil_AddSoilLayerEnvs(NSMPSCSchemaSoil* sse, int numberSoilLayers, ...);

void  NSMPSCSchemaSoil_EmptyList(NSMPSCSchemaSoil* sse);


EXTERN EXPORT  NSMPSCSchemaSoil*  NSMPSCSchemaSoil_Delete(NSMPSCSchemaSoil* sse);

EXTERN EXPORT  unsigned long  NSMPSCSchemaSoil_GetCount(NSMPSCSchemaSoil* sse);

EXTERN EXPORT  NSMPSSoilLayerEnvironment*  NSMPSCSchemaSoil_GetNext(NSMPSCSchemaSoil* sse, NSMPSSoilLayerEnvironment* last);

EXTERN EXPORT  unsigned long  NSMPSCSchemaSoil_GetList(NSMPSCSchemaSoil* sc, unsigned long numberSoilLayersRequested, NSMPSSoilLayerEnvironment** list);
#endif
