#pragma once
#ifndef Q2PSIL_H
#define Q2PSIL_H


#include "LINKAGE.h"

#include "C.Enumerations.h"

#include "C.Q2.Environment.h"

#include "C.Q2.Overland.Cell.h"

#include "C.PS.Environment.Cell.h"

#include "C.LandUseProperties.h"

#include "C.PS.Cell.h"





// Allow user to disable all overland/soil mass transfer

EXPORT void NSMOverlandCell_EnableSoilExchange(int enable);




// Changed the specie enumeration to exclude word "aquatic"


EXPORT  double NSMGetSedimentWaterMassFlux(NSMOverlandCell* ovCell, 
                    NSMCPSCell* psCell,
                    enum NSM_SPECIES_ENUM specieEnum,
                    NSMNutrientProperties* np,
                    NSMLandUseProperties* lup, 
                    NSMMetEnvironment* me);


// 3-2009: Added to ease computational burden associated with using NSMMetEnvironmet struct!

EXPORT  double NSMGetOverlandSoilMassFlux(NSMOverlandCell* ovCell, 
                    NSMCPSCell* psCell,
                    enum NSM_SPECIES_ENUM specieEnum,
                    NSMNutrientProperties* np,
                    NSMLandUseProperties* lup, 
                    double rainfallRate);

#endif
