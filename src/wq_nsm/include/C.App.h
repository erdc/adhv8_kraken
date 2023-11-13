//C.App.h
#pragma once

#ifndef APP_H

#define APP_H

#include "LINKAGE.h"



typedef struct 
{
	int useDefaultVegatativeCover;

	int useDefaultSoilCharacteristics;

	int useDefaultIOLibrary;

} NSMAppEnvironment;



EXPORT  NSMAppEnvironment*  NSMCreateAppEnvironment();

EXPORT  void  NSMInitAppEnvironment(NSMAppEnvironment* p);

EXPORT  void  NSMSetAppEnvironment(NSMAppEnvironment* p);

#endif