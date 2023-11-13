#pragma once

#include <map>

#include "NSMK_CORE.h"
#include "NSMK_StoichiometryRatios.h"
#include "NSMK_Pool.h"
#include "NSMK_AquaticPools.h"
#include "NSMK_Pathways.h"
#include "NSMK_OxygenConstraints.h"









class AquaticReactor : public NSMKAquaticReactor
{
    // F I E L D S - User-specified 

    public: 

    double area;

    double slope;

    double hydraulicRadius;



    // F I E L D S - NSMK-specified

    public:

    static std::map<NSMKAquaticReactor*, AquaticReactor*> globalObjMap;  

    AquaticPoolMap  poolMap;    // std::map<enum NSMKPoolEnum, AquaticPool*>
   
    PathwayMap  pathwayMap;     // std::map<enum EnumPathway, Pathway*>
  
    StoichiometryRatios* sr;

    OxygenConstraints  constraints;




    // M E T H O D S

    public:

    AquaticReactor(NSMKPoolAttributeGroup* NSMK_pag, double area, double slope,  double hydraulicRadius, int numberParticles);


    NSMKPool* PartitionByPool(const enum NSMKPoolEnum, const PoolAttributeGroup*);

    NSMKPool** Partition(const PoolAttributeGroup*);


    void DmDtSettling(const PoolAttributeGroup*);

    void DmDtKinetic(const NSMKMetProperties*, const NSMKPoolParameters*);


    public:

    AquaticPool* operator[](enum NSMKPoolEnum);

    const AquaticPool* operator[](enum NSMKPoolEnum) const;


    private:
        
    //void Configure(const PoolAttributeGroup*);
};
