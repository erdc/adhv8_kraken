#pragma once
#ifndef NUTRIENT_PROPERTIES_H

#define NUTRIENT_PROPERTIES_H


#include <stdlib.h>

#include "LINKAGE.h"



typedef struct 
{
	// Dissolved Transfer Properties   Meters/day
    double NO2_Km;
	double NO3_Km;
	double NH4_Km;
	double ON_Km;
	double OP_Km;
	double DP_Km;

	// Sorption Properties
	double NO2_KOC;
	double NO3_KOC;
	double NH4_KOC;
	double ON_KOC;
	double OP_KOC;
	double DP_KOC;

} NSMNutrientProperties;			
	




// F U N C T I O N   P R O T O T Y P E S

EXTERN EXPORT  NSMNutrientProperties*  NSMNutrientProperties_Create();

EXTERN EXPORT  NSMNutrientProperties*  NSMNutrientProperties_Delete(NSMNutrientProperties*);

EXTERN EXPORT  void  NSMNutrientProperties_Init(NSMNutrientProperties* se);

#endif