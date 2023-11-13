#pragma once
#ifndef LAND_USE_PROPERTIES_H

#define LAND_USE_PROPERTIES_H


#include <stdlib.h>

#include "LINKAGE.h"



typedef struct 
{

	double soilDetachability;		// Soil Detachability parameter. Units: [kg/m^3]

}  NSMLandUseProperties;	



// F U N C T I O N   P R O T O T Y P E S

EXTERN EXPORT  NSMLandUseProperties* NSMLandUseProperties_Create();

EXTERN EXPORT  NSMLandUseProperties*  NSMLandUseProperties_Delete(NSMLandUseProperties* lup);

EXTERN EXPORT  void NSMLandUSeProperties_Init(NSMLandUseProperties* lup);



#endif