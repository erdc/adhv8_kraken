#pragma once

#include <map>

#include "NSMK_CORE.h"
#include "NSMK_StoichiometryRatios.h"
#include "NSMK_Pool.h"
#include "NSMK_BenthicPool.h"
#include "NSMK_Pathways.h"
#include "NSMK_OxygenConstraints.h"








class BenthicReactor : public NSMKBenthicReactor
{
    // F I E L D S - User-specified 

    public: 

    double area;


    // F I E L D S - NSMK-specified

    public:
    static std::map<NSMKBenthicReactor*, BenthicReactor*> globalObjMap;  

    BenthicPoolMap  poolMap;
   
    PathwayMap  pathwayMap;
  
    StoichiometryRatios*  sr;

    OxygenConstraints  constraints;

        

    // M E T H O D S

    public:

    BenthicReactor(NSMKPoolAttributeGroup* NSMK_pag, double area, int numberParticles);


    NSMKPool* PartitionByPool(const enum NSMKPoolEnum, const PoolAttributeGroup*);

    NSMKPool** Partition(const PoolAttributeGroup*);


    void DmDtResuspension(const PoolAttributeGroup*);

    void DmDtKinetic(const NSMKMetProperties*, const NSMKPoolParameters*);


    public:

    BenthicPool* operator[](enum NSMKPoolEnum);

    const BenthicPool* operator[](enum NSMKPoolEnum) const;


    private:

    //void Configure(const PoolAttributeGroup*);
};
