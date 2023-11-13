#pragma once

#include <map>

#include "NSMK_CORE.h"





struct NSMKMetProperties;

struct NSMKPoolParameters;

class AquaticReactor;

class BenthicReactor;




enum EnumPathway
{
    // Aquatic pathways

    EnumPathwayPhytoPhoto,

    EnumPathwayPhytoResp,

    EnumPathwayPhytoDeath,

    EnumPathwayBotAlgPhoto,

    EnumPathwayBotAlgResp,

    EnumPathwayBotAlgDeath,

    EnumPathwayBotAlgExN,

    EnumPathwayBotAlgUpN,

    EnumPathwayBotAlgExP,

    EnumPathwayBotAlgUpP,

    EnumPathwayPOMDiss,

    EnumPathwaySlowCHydr,

    EnumPathwaySlowCOxid,

    EnumPathwayFastCOxid,

    EnumPathwayDONHydr,

    EnumPathwayDONMineralization,

    EnumPathwayNitrif,

    EnumPathwayDenitr,

    EnumPathwayDOPHydr,

    EnumPathwayDOPMineralization,

    EnumPathwayInorganicSuspSolidsSettl,

    EnumPathwayOxReaer,


    // Carbon cycle pools

    EnumPathwayDICAtmosphericFlux,

    EnumPathwayCO2Volatilization,

    EnumPathwayDOCMineralization,

    EnumPathwayPOCHydr,


    // Benthic pathways

    EnumPathwayBenthic_DOC_Mineralization,

    EnumPathwayBenthic_DOC_NO3Denitrification,

    EnumPathwayBenthic_POC_Hydrolysis,

    EnumPathwayBenthic_DON_Mineralization,

    EnumPathwayBenthic_NH4_Nitrification,

    EnumPathwayBenthic_NO3_Denitrification,

    EnumPathwayBenthic_PON_Hydrolysis,

    EnumPathwayBenthic_DOP_Mineralization,

    EnumPathwayBenthic_POP_Hydrolysis 
};










class Pathway
{
    protected: bool  enabled;

    public: 

    double  dmdt;

    Pathway();

    virtual double DmDt(const AquaticReactor*,  const NSMKMetProperties*,  const NSMKPoolParameters*, const NSMKSedProperties*);

    virtual double DmDt(const BenthicReactor*,  const NSMKMetProperties*,  const NSMKPoolParameters*, const NSMKSedProperties*);

    protected:
    static double tCor(TD_COEFFICIENT c, double t);
};

typedef std::map<enum EnumPathway, Pathway*> PathwayMap;





//<AQUATIC PATHWAYS>

class PathwayPhytoPhoto : public Pathway
{
    public:
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};


class PathwayPhytoResp : public Pathway
{
    public:
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};


class PathwayPhytoDeath : public Pathway
{
    public:
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};






class PathwayBotAlgPhoto : public Pathway
{
    public:
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};


class PathwayBotAlgResp : public Pathway
{
    public: 
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};


class PathwayBotAlgDeath : public Pathway
{
    public:
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};

class PathwayBotAlgExN : public Pathway
{
    public:
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};


class PathwayBotAlgUpN : public Pathway
{    
    public:
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};


class PathwayBotAlgExP : public Pathway
{
    public:
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};


class PathwayBotAlgUpP : public Pathway
{
    public:
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};




class PathwayPOMDiss : public Pathway
{
    public:     
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};


class PathwaySlowCHydr : public Pathway
{
    public:
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};


class PathwaySlowCOxid : public Pathway
{
    public:
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};


class PathwayFastCOxid : public Pathway
{
    public:
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};




class PathwayDONHydr : public Pathway
{
    public:
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};



class PathwayDONMineralization : public Pathway
{
     public: 
     virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};



class PathwayNitrif : public Pathway
{
    public: 
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};




class PathwayDenitr : public Pathway
{
    public:
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};






class PathwayDOPHydr : public Pathway
{
    public: 
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};





class PathwayDOPMineralization : public Pathway
{
    public: 
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};





//class PathwayInorganicSuspSolidsSettl : public Pathway
//{
//    public: virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
//
//};


class PathwayOxReaer : public Pathway
{
    public:
    virtual double DmDt(const AquaticReactor* ar, const NSMKMetProperties* mp,  const NSMKPoolParameters* a, const NSMKSedProperties* aspe);
};





class PathwayDICAtmosphericFlux : public Pathway
{
    public:
    static double HenrysConstantCO2;

    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};



class PathwayCO2Volatilization : public Pathway
{
    public:
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};



class PathwayDOCMineralization :  public Pathway
{
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};



class PathwayPOCHydr : public Pathway
{
    virtual double DmDt(const AquaticReactor* ar,  const NSMKMetProperties* mp,  const NSMKPoolParameters* ap, const NSMKSedProperties* asp);
};


//</AQUATIC PATHWAYS>









//<BENTHIC PATHWAYS>

// <Benthic Carbon pathways>

class PathwayBenthic_DOC_Mineralization : public Pathway
{
    virtual double DmDt(const BenthicReactor*,  const NSMKMetProperties*,  const NSMKPoolParameters*, const NSMKSedProperties*);
};


class PathwayBenthic_DOC_NO3Denitrification : public Pathway
{
    virtual double DmDt(const BenthicReactor*,  const NSMKMetProperties*,  const NSMKPoolParameters*, const NSMKSedProperties*);
};


class PathwayBenthic_POC_Hydrolysis : public Pathway
{
    virtual double DmDt(const BenthicReactor*,  const NSMKMetProperties*,  const NSMKPoolParameters*, const NSMKSedProperties*);
};

// </Benthic Carbon pathways>


// <Benthic Nitrogen pathways>

class PathwayBenthic_DON_Mineralization : public Pathway
{
    virtual double DmDt(const BenthicReactor*,  const NSMKMetProperties*,  const NSMKPoolParameters*, const NSMKSedProperties*);
};


class PathwayBenthic_NH4_Nitrification : public Pathway
{
    virtual double DmDt(const BenthicReactor*,  const NSMKMetProperties*,  const NSMKPoolParameters*, const NSMKSedProperties*);
};


class PathwayBenthic_NO3_Denitrification : public Pathway
{
    virtual double DmDt(const BenthicReactor*,  const NSMKMetProperties*,  const NSMKPoolParameters*, const NSMKSedProperties*);
};


class PathwayBenthic_PON_Hydrolysis : public Pathway
{
    virtual double DmDt(const BenthicReactor*,  const NSMKMetProperties*,  const NSMKPoolParameters*, const NSMKSedProperties*);
};

// </Benthic Nitrogen pathways>


// <Benthic Phosphorus pathways>

class PathwayBenthic_DOP_Mineralization : public Pathway
{
    virtual double DmDt(const BenthicReactor*,  const NSMKMetProperties*,  const NSMKPoolParameters*, const NSMKSedProperties*);
};


class PathwayBenthic_POP_Hydrolysis : public Pathway
{
    virtual double DmDt(const BenthicReactor*,  const NSMKMetProperties*,  const NSMKPoolParameters*, const NSMKSedProperties*);
};

// </Benthic Phosphorus pathways>


//</BENTHIC PATHWAYS>