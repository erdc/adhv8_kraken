#pragma once
#ifndef Q2_ENVIRONMENT_H

#define Q2_ENVIRONMENT_H


#include <stdio.h>

#include "LINKAGE.h"

#include "C.Enumerations.h"

#include "C.METEnvironment.h"

#include "C.NutrientProperties.h"

#include "C.CSchema.Particle.h"

#include "C.CSchema.ParticleDistribution.h"









#define NSM_NUMBER_AQUATIC_SPECIES 9




typedef struct
{
	double  conc;
	double  mass;
	double  dmdt;
     
    enum NSM_PARTITIONABLE_ENUM   partitionable; // New in 4.6 & later

	double  fractionDissolved;
	double  fractionBound;

	unsigned long  numberSorbingParticles;
    double*  fractionPerSorbingParticle;

	double  fsum; 	// Sorption validation check. Should be very near 1.0.

}  NSM_AQUATIC_SPECIES;








typedef struct 
{
	long numberSpecies;

	NSM_AQUATIC_SPECIES aquatic_NO2;

	NSM_AQUATIC_SPECIES aquatic_NO3;

	NSM_AQUATIC_SPECIES aquatic_NH4;

	NSM_AQUATIC_SPECIES aquatic_ON;

	NSM_AQUATIC_SPECIES aquatic_OP;

	NSM_AQUATIC_SPECIES aquatic_DP;

	NSM_AQUATIC_SPECIES aquatic_Alg;

	NSM_AQUATIC_SPECIES aquatic_CBOD;

	NSM_AQUATIC_SPECIES aquatic_DO;

	NSM_AQUATIC_SPECIES* list[NSM_NUMBER_AQUATIC_SPECIES];

} NSM_AQUATIC_SPECIES_NET;








typedef struct 
{
	double airTemp;			// Units: C
	
	double solarRad;		// Units: (Watts/m^2) 
	
	double rh;				// Relative humidity. Unitless.
	
	double ap;				// Atmospheric pressure. Units: MilliBars
	
	double windSpeed;		// Mean wind speed above site. Units: m/s

} NSM_MET_OBS;







typedef struct
{
	double rc;      			/* reaction_coefficient */
	
	double etcc;				/* Temperature Correction Coefficient (theta) */

} NSM_COEFFICIENT;






typedef struct
{
	NSM_COEFFICIENT alpha0;		/*TC: alpha0 := Fraction of chloraquatic_OPhil-A to Algal biomass.
								Units: (ug-Chl-A/mg-A). Range: {10-100} */

	NSM_COEFFICIENT alpha1;		/*TC: alpha1 := Fraction of algal biomass that is nitrogen.
								Units: (mg-N/mg-A). Range: {0.07-0.09} */

	NSM_COEFFICIENT alpha2;		/*TC: alpha2 := Fraction of algal biomass that is Phosphorus.
								Units: (mg-N/mg-A). Range: {0.01-0.02} */
	 
	NSM_COEFFICIENT alpha3;		/*TC: alpha3 := O2 production per unit of algal growth.
								Units: (mg-O/mg-A). Range: {1.4-1.8} */

	NSM_COEFFICIENT alpha4;		/*TC: alpha4 := O2 uptake per unit algae respired.
								Units: (mg-O/mg-A). Range: {1.6-2.3} */

	NSM_COEFFICIENT alpha5;		/*TC: alpha5 := O2 uptake per unit NH3 oxidized.
								Units: (mg-O/mg-NH3). Range: {3.0-4.0} */

	NSM_COEFFICIENT alpha6;		/*TC: alpha6 := O2 uptake per unit aquatic_NO2 oxidized.
								Units: (mg-O/mg-aquatic_NO2). Range: {1.0-1.14} */



	NSM_COEFFICIENT beta1;		/*TVSV: beta1 := Ammonia decay (rate constant for oxidation of NH3 (ammonia) to aquatic_NO2): 
								Units: (1/day). Range: {0.1-1.0} */

	NSM_COEFFICIENT beta2;		/*TVSV: beta2 := Nitrite decay (rate constant for the biological oxidation of aquatic_NO2 to aquatic_NO3:
								Units: (1/day). Range: {0.2-2.0} */

	NSM_COEFFICIENT beta3;		/*TVSV: beta3 := Nitrite decay (rate constant for the hydrolysis of organic N to ammonia). 
								Units: (1/day). Range: {0.02-0.4} */

	NSM_COEFFICIENT beta4;		/*TVSV: beta4 := Organic Phosphorus decay (rate constant for the decay of organic-P to dissolved P).
								Units: (1/day). Range: {0.01-0.7} */



	NSM_COEFFICIENT K1;			/*TVSV: K1 := Carbonaceous Deoxygenation rate constant.
								Units: (1/day). Range: {0.02-3.4} */

	NSM_COEFFICIENT K2;			/*TVSV: K2 := Reaeration Rate Coefficient.
								Units: (1/day). Range: {0.0-100.0} */

	NSM_COEFFICIENT K3;			/*TVSV: K3 := aquatic_CBOD settling rate.
								Units: (1/day). Range: {-0.36-0.36} */

	NSM_COEFFICIENT K4;			/*TVSV: K4 := SOD (Benthic Oxygen) uptake. 
								Units: (mg-O/ft^2-day). Range: {variable} */



	NSM_COEFFICIENT mu;			/*TV: mu := Algal growth rate. 
								Units: (1/day). Range: {1.0-3.0} */

	NSM_COEFFICIENT rho;		/*TV: rho := Algal respiration rate. 
								Units (1/day). Range: {0.05,0.5} */



	NSM_COEFFICIENT sigma1;		/*TVSV: sigma3 := Algal settling rate. 
								Units (ft/day).. Range: {variable} */

	NSM_COEFFICIENT sigma2;		/*TVSV: sigma2 := Benthos source rate for dissolved Phosphorus. 
								Units: (mg/ft^2-day). Range: {variable} */

	NSM_COEFFICIENT sigma3;		/*TVSV: sigma3 := Benthos source rate for ammonia Nitrogen. 
								Units: (mg/ft^2-day). Range: {variable} */

	NSM_COEFFICIENT sigma4;		/*TVSV: sigma4 := Organic Nitrogen settling rate: 
								Units: (1/day). Range: {0.001-0.1} */

	NSM_COEFFICIENT sigma5;		/*TVSV: sigma5 := Organic Phosphorus settling rate: 
								Units: (1/day). Range: {0.001-0.1} */


	double Kl;					/* Michaelis-Menton half-saturation constant for light. 
								Range: {0.02-0.10}. Units: (Btu/(ft^2-min)).	 */

	double Kn;					/* Michaelis-Menton half-saturation constant for nitrogen.
								Range: {0.01-0.02}. Units: (mg-N/L).		 */

	double Kp;					/* Michaelis-Menton half-saturation constant for phosphorus. 
								Range: {0.01-0.20}. Units: (mg-P/L).      */


	double Pn;					/* Algal perference factor for ammonia. 
								Range: {0.0-1.0}. Units: nondimensional.	 */


	double lambda0;				/* Non-algal self-shading coefficient. 
								Range: {variable}. Units: (1/ft). */

	double lambda1;				/* Linear algal self-shading coefficient. 
								Range: {0.002-0.02}. Units: ( (1/ft) / (ug Chl-A/L) ). */

	double lambda2;				/* Nonlinear algal self-shading coefficient. 
								Range: {0.0165}. Units: ( (1/ft) / (ug Chl-A/L)^(2/3) ). */

	double min_do_conc;			/* Min [aquatic_DO] for kinetics calculations */
	
} NSM_AQUATIC_ENVIRONMENT;



//  NSM_AQUATIC_ENVIRONMENT Prototypes

EXPORT  NSM_AQUATIC_ENVIRONMENT*  NSM_AQUATIC_ENVIRONMENT_Create();

EXPORT  NSM_AQUATIC_ENVIRONMENT*  NSM_AQUATIC_ENVIRONMENT_Init(NSM_AQUATIC_ENVIRONMENT*);

EXPORT  NSM_AQUATIC_ENVIRONMENT*  NSM_AQUATIC_ENVIRONMENT_Delete(NSM_AQUATIC_ENVIRONMENT*);





//  F U N C T I O N   P R O T O T Y P E S  

EXPORT  void NSM_init();







// Soil-Water exchange environments
//EXPORT  NSMSoilWaterExchangeEnvironment* NSMCreateSoilWaterExchangeEnvironment();
//EXPORT  void NSMInitSoilWaterExchangeEnvironment(NSMSoilWaterExchangeEnvironment* p);

#endif







//                        I N F O   &   N O T E S
//
//
//
//           INTER-COMPARTMENT SPECIE MAPPING
//
//        OVERLAND               SOIL                         INDEX
//        oc.NH4 <-------------> ne.NH4                          0
//        oc.NO3 <-------------> ne.NO3                          1
//        oc.ON  <-------------> ne.OrganicFresh                 2
//        oc.DP  <-------------> pe.MMineralSolution             3
//        oc.OP  <-------------> pe.OrganicFresh                 4
//
//
//        PO4    <-------------> pe.MineralActive
//        Next version drops PO4 though!!!
