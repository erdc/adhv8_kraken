#pragma once
#ifndef C_Q2_OVERLAND_CELL_H

#define C_Q2_OVERLAND_CELL_H

#include <stdio.h>

#include "LINKAGE.h"

#include "C.App.h"

#include "C.METEnvironment.h"

#include "C.NutrientProperties.h"

#include "C.CSchema.Particle.h"

#include "C.CSchema.ParticleDistribution.h"

#include "C.Sediment.h"

#include "C.Q2.Environment.h"

#include "C.Enumerations.h"




typedef struct 
{
	double waterTemp;		// Temperature of water in cell.  Units: C

	double depth;			// Depth of water in cell.  Units: m

    double minimumDepth;    // Minimum depth to run kinetics

	double vol;				// Volume of water in cell.  Units: m^3

	NSM_AQUATIC_SPECIES_NET aqSpecies;

	double dissolvedOrganicCarbon;

	double dissolvedOrganicCarbon_FOC; 	// Logically = 1.0.

} NSMOverlandCell;





// F U N C T I O N   P R O T O T Y P E S

EXPORT  NSMOverlandCell* NSMOverlandCell_Create();

EXPORT  void  NSM_OVERLAND_Init(NSMOverlandCell* overlandCell);

EXPORT  NSMOverlandCell*  NSMOverlandCell_Delete(NSMOverlandCell* overlandCell);



EXPORT  void  NSMOverlandCell_EnableSystem(int enable);

EXPORT  void  NSMOverlandCell_Kinetics(NSMOverlandCell* pCell, NSMMetEnvironment* me, NSM_AQUATIC_ENVIRONMENT* ae);

EXPORT  double  NSMOverlandCell_CalculateMaximumDissolvedOxygenConc(double waterTempC);

EXPORT  void  NSMOverlandCell_SetMinimumDepth(NSMOverlandCell* pCell, double depth);



// Two techniques for calculating specie partitioning

EXPORT  void  NSMOverlandCell_CalculateSedimentPartitioning(
                    NSMOverlandCell* pCell, 
                    NSMNutrientProperties* np, 
                    NSMSedimentEnvironment* se);             // If using NSMSedimentEnvironment


EXPORT  void  NSMOverlandCell_CalculateSpecieSedimentPartitioning(
                    NSMOverlandCell* pCell,
                    enum NSM_SPECIES_ENUM specie, 
                    NSMNutrientProperties* np, 
                    NSMSedimentEnvironment* se);            // If using NSMSedimentEnvironment


EXPORT  void  NSMOverlandCell_CalculateSpeciePartitioning(
                    NSMOverlandCell* pCell, 
                    enum NSM_SPECIES_ENUM specie, 
                    NSMNutrientProperties* np, 
                    NSMCSchemaParticleDistribution* pd);


#endif