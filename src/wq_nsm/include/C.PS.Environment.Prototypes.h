#pragma once


#include "LINKAGE.h"

#include "C.METEnvironment.h"

#include "C.PS.Environment.Cell.h"

#include "C.PS.CSchema.Flora.h"

#include "C.PS.Environment.Nutrients.h"

#include "C.PS.Environment.Particle.h"

#include "C.PS.Environment.SoilLayer.h"

#include "C.PS.CSchema.Soil.h"

#include "C.PS.Cell.h"







//
//// N U T R I E N T    E N V I R O N M E N T    P R O T O T Y P E S
//
//void  initCarbonEnvironment(NSMPSSoilLayerEnvironment* sle);
//
//void  initNitrogenEnvironment(NSMPSSoilLayerEnvironment* sle);
//
//void  initPhosphorusEnvironment(NSMPSSoilLayerEnvironment* sle);
//






// P A R T I C L E   E N V I R O N M E N T   P R O T O T Y P E S 

EXTERN EXPORT  NSMSoilParticleEnvironment* NSMCreateSoilParticleEnvironment();

EXTERN EXPORT  NSMSoilParticleDistributionEnvironment* NSMCreateSoilParticleDistributionEnvironment(long numberOfParticles);










