#pragma once
#ifndef C_PS_CSCHEMA_FLORA_H
#define C_PS_CSCHEMA_FLORA_H


#include "LINKAGE.h"

#include "C.Library.h"

#include "C.PS.Environment.Plant.h"






struct NSMPSCSchemaFlora
{
	unsigned long numberPlants;

	NSMListNode* plantList;
};

typedef struct NSMPSCSchemaFlora NSMPSCSchemaFlora;





// F L O R A   S C H E M A   P R O T O T Y P E S

EXTERN EXPORT  NSMPSCSchemaFlora*  NSMPSCSchemaFlora_Create();

EXTERN EXPORT  void  NSMPSCSchemaFlora_Init(NSMPSCSchemaFlora* fe);

EXTERN EXPORT  unsigned long  NSMPSCSchemaFlora_Add(NSMPSCSchemaFlora* fe, NSMPSPlantEnvironment* pe);

EXTERN EXPORT  unsigned long  NSMPSCSchemaFlora_AddList(NSMPSCSchemaFlora* fe, int numberPlants, NSMPSPlantEnvironment** plantEnvs);

EXTERN EXPORT  unsigned long  NSMPSCSchemaFlora_AddPlantEnvs(NSMPSCSchemaFlora* fe, int numberPlants, ...);

void  NSMPSCSchemaFlora_EmptyList(NSMPSCSchemaFlora* fe);

EXTERN EXPORT  NSMPSCSchemaFlora*  NSMPSCSchemaFlora_Delete(NSMPSCSchemaFlora* fe);

EXTERN EXPORT  unsigned long  NSMPSCSchemaFlora_GetCount(NSMPSCSchemaFlora* fe);

EXTERN EXPORT  NSMPSPlantEnvironment*  NSMPSCSchemaFlora_GetNext(NSMPSCSchemaFlora* fe, NSMPSPlantEnvironment* last);

EXTERN EXPORT  unsigned long  NSMPSCSchemaFlora_GetList(NSMPSCSchemaFlora* fe, unsigned long numberPlantsRequested, NSMPSPlantEnvironment** list);
#endif
