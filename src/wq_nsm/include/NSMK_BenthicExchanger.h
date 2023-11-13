#pragma once

#include <map>

#include "NSMK_CORE.h"

#include "NSMK_AquaticReactor.h"

#include "NSMK_BenthicReactor.h"

#include "NSMK_PoolAttributeGroup.h"




class BenthicExchanger : public NSMKBenthicExchanger
{
    // F I E L D S

    public:

    static std::map<NSMKBenthicExchanger*, BenthicExchanger*> globalObjMap;  

    AquaticReactor* ar;

    BenthicReactor* br;


    // C D T O R S

    public:

    BenthicExchanger(NSMKAquaticReactor* pa, NSMKBenthicReactor* pb);

    ~BenthicExchanger();


    // M E T H O D S

    void TransferSettlingToBenthicReactor(PoolAttributeGroup* pag);

    void TransferResuspensionToAqauticReactor(PoolAttributeGroup* pag);

    void DoDiffusiveTransfer(PoolAttributeGroup* pag);
};
