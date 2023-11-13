#pragma once
#ifndef C_PS_ENVIRONMENT_PARTICLE_H
#define C_PS_ENVIRONMENT_PARTICLE_H

#include "LINKAGE.h"



#define MAX_NUMBER_SOIL_PARTICLES 25


      

struct NSMSoilParticleEnvironment
{
	double  fraction;

	double  bindingEffectivenessCoefficient;

	double  FOC;

	double  particleInteractionParameter;	
};

typedef struct NSMSoilParticleEnvironment NSMSoilParticleEnvironment;







struct  NSMSoilParticleDistributionEnvironment
{
	unsigned long numberOfSoilParticles;

	NSMSoilParticleEnvironment  soilParticles[MAX_NUMBER_SOIL_PARTICLES];
};

typedef struct NSMSoilParticleDistributionEnvironment NSMSoilParticleDistributionEnvironment;
#endif
