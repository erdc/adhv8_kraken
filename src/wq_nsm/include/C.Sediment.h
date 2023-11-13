#pragma once
#ifndef C_SEDIMENT_H

#define C_SEDIMENT_H


#include "LINKAGE.h"

#include "C.Enumerations.h"


typedef struct 
{    
    unsigned long numberSediments;

	double  bindingEffectivenessCoefficient;	// De: DOC binding


    // Following are arrays of size: numberSediments

    int*  tag;                              // Array of ints for host use - NOT USED by NSM

    enum NSM_SEDIMENT_SORBING_ENUM* sorbing;         // Specifies if sediment class is capable of sorbing

	double*  concentration;			        // Aquatic : M/L^3

 	double*  FOC;							// Fraction Organic Carbon

	double*  particleInteractionParameter;	// DiToro's Vx parameter 

    double*  P;                             // Spread fraction - internal

} NSMSedimentEnvironment;




// NSMSedimentEnvironment   P R O T O T Y P E S

EXTERN EXPORT  NSMSedimentEnvironment*  NSMSedimentEnvironment_Create(unsigned long numberSediments);

EXTERN EXPORT  void   NSMSedimentEnvironment_Init(NSMSedimentEnvironment* p, unsigned long numberSediments);

EXTERN EXPORT  NSMSedimentEnvironment*  NSMSedimentEnvironment_Clone(NSMSedimentEnvironment* p);

EXTERN EXPORT  NSMSedimentEnvironment*  NSMSedimentEnvironment_Delete(NSMSedimentEnvironment* p);

#endif



