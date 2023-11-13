#pragma once
#ifndef C_CSCHEMA_PARTICLE_DISTRIBUTION_H

#define C_CSCHEMA_PARTICLE_DISTRIBUTION_H

#include <stdio.h>

#include "LINKAGE.h"

#include "C.Library.h"

#include "C.CSchema.Particle.h"



struct NSMCSchemaParticleDistribution
{
	unsigned long  numberOfParticles;

	double  bindingEffectivenessCoefficient;	// De: DOC binding

	NSMListNode* particleList;
};

typedef struct NSMCSchemaParticleDistribution NSMCSchemaParticleDistribution;






// P R O T O T Y P E S

EXTERN EXPORT  NSMCSchemaParticleDistribution*  NSMCSchemaParticleDistribution_Create();

EXTERN EXPORT  void  NSMCSchemaParticleDistribution_Init(NSMCSchemaParticleDistribution* pcs);

void  NSMCSchemaParticleDistribution_EmptyList(NSMCSchemaParticleDistribution* pcs);

EXTERN EXPORT  unsigned long  NSMCSchemaParticleDistribution_GetCount(NSMCSchemaParticleDistribution* pcs);

EXTERN EXPORT  NSMCSchemaParticle*  NSMCSchemaParticleDistribution_GetNext(NSMCSchemaParticleDistribution* pcs, NSMCSchemaParticle* last);

EXTERN EXPORT  unsigned long  NSMCSchemaParticleDistribution_Add(NSMCSchemaParticleDistribution* pcs, NSMCSchemaParticle* p);

EXTERN EXPORT  unsigned long  NSMCSchemaParticleDistribution_AddList(NSMCSchemaParticleDistribution* pcs, int numberParticles, NSMCSchemaParticle** particles);

EXTERN EXPORT  unsigned long   NSMCSchemaParticleDistribution_AddParticles(NSMCSchemaParticleDistribution* pcs, int numberParticles, ...);

EXTERN EXPORT  NSMCSchemaParticleDistribution*  NSMPSCSchemaParticleDistribution_Delete(NSMCSchemaParticleDistribution* pcs);


#endif
