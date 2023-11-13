#pragma once
#ifndef C_PS_ENVIRONMENT_PLANT_H
#define C_PS_ENVIRONMENT_PLANT_H

#include "LINKAGE.h"




struct NSMPSPlantEnvironment
{
	char name[200];

	double coverFraction;

    struct NSMPSPlantEnvironment*  next;
};

typedef struct NSMPSPlantEnvironment NSMPSPlantEnvironment;






// P L A N T   E N V I R O N M E N T   P R O T O T Y P E S 


EXTERN EXPORT  NSMPSPlantEnvironment*  NSMPSPlantEnvironment_Create(char* name, double coverFraction);


EXTERN EXPORT  NSMPSPlantEnvironment*  NSMPSPlantEnvironment_Delete(NSMPSPlantEnvironment* pe);


EXTERN EXPORT  void  NSMPSPlantEnvironment_Init(NSMPSPlantEnvironment* pe, char* name, double coverFraction);


EXTERN EXPORT  NSMPSPlantEnvironment*  NSMPSPlantEnvironment_Clone(NSMPSPlantEnvironment* pe);
#endif
