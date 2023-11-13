#pragma once


#include <map>

#include "NSMK_CORE.h"

#include "NSMK_Pool.h"

#include "NSMK_PartitionDistribution.h"

#include "NSMK_Pathways.h"






class AquaticPool;

class StoichiometryRatios;

typedef std::map<enum NSMKPoolEnum, AquaticPool*> AquaticPoolMap;


class AquaticPool : public Pool
{
    public:

    AquaticReactor* aquaticReactor;   // Ref to container


    public:

    AquaticPool();

    virtual ~AquaticPool();


    public:

    virtual void Configure(const PoolAttributeGroup*);

    virtual double Partition(const PoolAttributeGroup*);

    virtual double DmDtKinetic(const NSMKPoolParameters*) = 0;

    virtual double DmDtSettling(const PoolAttributeGroup*);
};

