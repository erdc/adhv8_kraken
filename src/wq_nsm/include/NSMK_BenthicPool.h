#pragma once


#include <map>

#include "NSMK_CORE.h"

#include "NSMK_Pool.h"

#include "NSMK_PartitionDistribution.h"





class BenthicPool;

class BenthicReactor;

class StoichiometryRatios;

typedef std::map<enum NSMKPoolEnum, BenthicPool*> BenthicPoolMap;




class BenthicPool : public Pool
{
    public:

    BenthicReactor* benthicReactor;   // Ref to container


    public:

    BenthicPool();

    virtual ~BenthicPool();


    public:

    virtual void Configure(const PoolAttributeGroup*);

    virtual double Partition(const PoolAttributeGroup*);

    virtual double DmDtKinetic(const NSMKPoolParameters*) = 0;

    virtual double DmDtResuspension(const PoolAttributeGroup*);
};


