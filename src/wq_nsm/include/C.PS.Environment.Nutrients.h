#pragma once
#ifndef C_PS_ENVIRONMENT_NUTRIENTS_H
#define C_PS_ENVIRONMENT_NUTRIENTS_H

#include "LINKAGE.h"





struct NSMPSPool
{
	double  mass;

	double  fractionDissolved;

	double  fractionSolid;

	double  fsum;		// Validity check
};

typedef struct NSMPSPool NSMPSPool;







struct NSMPSCarbonEnvironment
{
	double fractionCarbonInSoil;

	NSMPSPool organicActive;

	NSMPSPool organicStable;
};

typedef struct NSMPSCarbonEnvironment NSMPSCarbonEnvironment;








struct NSMPSNitrogenEnvironment
{
	NSMPSPool NH4;

	NSMPSPool NO3;

	NSMPSPool organicFresh;

	NSMPSPool organicActive;

	NSMPSPool organicStable;
};

typedef struct NSMPSNitrogenEnvironment NSMPSNitrogenEnvironment;








struct NSMPSPhosphorusEnvironment
{
	NSMPSPool organicFresh;

	NSMPSPool organicActive;

	NSMPSPool organicStable;

	NSMPSPool mineralActive;

	NSMPSPool mineralStable;

	NSMPSPool mineralSolution;

	double phosphorusAvailabilityIndex;			// Fraction fertilizer P in solution 
												// after period of fast reaction
};

typedef struct NSMPSPhosphorusEnvironment NSMPSPhosphorusEnvironment;






// N U T R I E N T    E N V I R O N M E N T    P R O T O T Y P E S

struct NSMPSSoilLayerEnvironment;

void  CarbonEnvironment_Init(struct NSMPSSoilLayerEnvironment* sle);

void  NitrogenEnvironment_Init(struct NSMPSSoilLayerEnvironment* sle);

void  PhosphorusEnvironment_Init(struct NSMPSSoilLayerEnvironment* sle);
#endif
