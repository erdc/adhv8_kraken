#pragma once


#include <map>

#include "NSMK_CORE.h"

#include "NSMK_PartitionDistribution.h"

#include "NSMK_PoolAttributeGroup.h"

#include "NSMK_Pathways.h"



class Pool;

class StoichiometryRatios;

typedef std::map<enum NSMKPoolEnum, Pool*> PoolMap;



class Pool : public NSMKPool
{
    public:

    PartitionDistribution* partition;

    PartitionDistribution* settling_dmdt;

    PartitionDistribution* resuspension_dmdt;


    public:

    virtual double Partition(const PoolAttributeGroup*) = 0;

    virtual double DmDtKinetic(const NSMKPoolParameters*) = 0;
};

