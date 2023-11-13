#pragma once
#ifndef C_CSCHEMA_PARTICLE_H

#define C_CSCHEMA_PARTICLE_H

#include <stdlib.h>

#include "LINKAGE.h"



typedef struct 
{
	double  aquaticConcentration;			// Aquatic : M/L^3

	double  FOC;							// Particle's Fraction Organic Carbon

	double  particleInteractionParameter;	// DiToro's Vx parameter

} NSMCSchemaParticle;





// P R O T O T Y P E S

EXTERN EXPORT  NSMCSchemaParticle*  NSMCSchemaParticle_Create();

EXTERN EXPORT  void  NSMCSchemaParticle_Init(NSMCSchemaParticle*);

EXTERN EXPORT  NSMCSchemaParticle*  NSMCSchemaParticle_Clone(NSMCSchemaParticle*);

EXTERN EXPORT  NSMCSchemaParticle*  NSMCSchemaParticle_Delete(NSMCSchemaParticle*);


#endif