#pragma once
#ifndef C_Q2_CHANNEL_CELL_H
#define C_Q2_CHANNEL_CELL_H

#include <stdio.h>

#include "LINKAGE.h"

#include "C.App.h"

#include "C.METEnvironment.h"

#include "C.NutrientProperties.h"

#include "C.CSchema.Particle.h"

#include "C.CSchema.ParticleDistribution.h"

#include "C.Sediment.h"

#include "C.Q2.Environment.h"





typedef struct 
{
	double waterTemp;		// Temperature of water in cell.  Units: C

	double depth;			// Depth of water in cell.  Units: m

    double minimumDepth;    // Minimum depth to run kinetics

	double vol;				// Volume of water in cell.  Units: m^3

	NSM_AQUATIC_SPECIES_NET aqSpecies;

	double dissolvedOrganicCarbon;

	double dissolvedOrganicCarbon_FOC; 	// Logically = 1.0.

} NSMChannelCell;






// F U N C T I O N   P R O T O T Y P E S

EXPORT  NSMChannelCell*  NSMChannelCell_Create();

EXPORT  void  NSMChannelCell_Init(NSMChannelCell* channelCell);

EXPORT  NSMChannelCell*  NSMChannelCell_Delete(NSMChannelCell* channelCell);



EXPORT  void  NSMChannelCell_EnableSystem(int enable);

EXPORT  void  NSMChannelCell_Kinetics(NSMChannelCell* pCell, NSMMetEnvironment* me, NSM_AQUATIC_ENVIRONMENT* ae);
	
EXPORT  double  NSMChannelCell_CalculateMaximumDissolvedOxygenConc(double waterTempC);

EXPORT  void  NSMChannelCell_SetMinimumDepth(NSMChannelCell* pCell, double depth);


// Two techniques for calculating specie partitioning

EXPORT  void  NSMChannelCell_CalculateSedimentPartitioning(
                    NSMChannelCell* pCell, 
                    NSMNutrientProperties* np, 
                    NSMSedimentEnvironment* se);     // If using NSMSedimentEnvironment


EXPORT  void  NSMChannelCell_CalculateSpecieSedimentPartitioning(
                    NSMChannelCell* pCell, 
                    enum NSM_SPECIES_ENUM specie, 
                    NSMNutrientProperties* np, 
                    NSMSedimentEnvironment* se);    // If using NSMSedimentEnvironment


EXPORT  void  NSMChannelCell_CalculateSpeciePartitioning(
                    NSMChannelCell* pCell, 
                    enum NSM_SPECIES_ENUM specie, 
                    NSMNutrientProperties* np, 
                    NSMCSchemaParticleDistribution* pd);


#endif
