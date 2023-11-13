#pragma once

#include <map>

#include "NSMK_CORE.h"





typedef std::map<enum NSMKPoolEnum, NSMKPoolAttribute*> PoolAttributeMap;

class PoolAttributeGroup : public NSMKPoolAttributeGroup
{
    public:

    PoolAttributeMap  poolAttributeMap;    // std::map<enum NSMKPoolEnum, NSMKPoolAttribute*>
   

    public:

    PoolAttributeGroup();

    ~PoolAttributeGroup();


    public:

    NSMKPoolAttribute* operator[](enum NSMKPoolEnum);

    const  NSMKPoolAttribute* operator[](enum NSMKPoolEnum) const;
};
