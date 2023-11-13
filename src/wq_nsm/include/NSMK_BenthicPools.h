#pragma once


#include <map>

#include "NSMK_CORE.h"

#include "NSMK_Pool.h"

#include "NSMK_BenthicPool.h"

#include "NSMK_PartitionDistribution.h"

#include "NSMK_Pathways.h"





class BenthicPoolDIC : public BenthicPool
{
    // Dissolved Inorganic Carbon 

    public: 

    BenthicPoolDIC(BenthicReactor* br);

    virtual double Partition(const PoolAttributeGroup* pag, const NSMKSedProperties* asp);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};




class BenthicPoolDIP : public BenthicPool
{    
    // Dissolved Inorganic Phosphorus

    public: 

    BenthicPoolDIP(BenthicReactor* br);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};




class BenthicPoolDOC : public BenthicPool
{
    // Dissolved Organic Carbon 

    public: 

    BenthicPoolDOC(BenthicReactor* br);

    virtual void Configure(const PoolAttributeGroup* pag);

    virtual double Partition(const PoolAttributeGroup* pag, const NSMKSedProperties* asp);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};




class BenthicPoolDON : public BenthicPool	
{
    // Dissolved Organic Nitrogen

    public:

    BenthicPoolDON(BenthicReactor* br);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class BenthicPoolDOP : public BenthicPool
{
    // Dissolved Organic Phosphorus

    public: 

    BenthicPoolDOP(BenthicReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};




class BenthicPoolDOX : public BenthicPool
{
    // Dissolved Oxygen

    public: 

    BenthicPoolDOX(BenthicReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};



class BenthicPoolNH4 : public BenthicPool
{
    // Ammonium Nitrogen

    public:

    BenthicPoolNH4(BenthicReactor*);

    virtual double DmDtKinetic(const NSMKPoolParameters*);
};





class BenthicPoolNO3 : public BenthicPool
{
    // Nitrate Nitrogen

    public:

    BenthicPoolNO3(BenthicReactor*);

    virtual double DmDtKinetic(const NSMKPoolParameters*);
};





class BenthicPoolPOC : public BenthicPool
{
    // Particulate Organic CARBON

    public: 

    BenthicPoolPOC(BenthicReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};




	

class BenthicPoolPON : public BenthicPool
{
    // Particulate Organic Nitrogen

    public: 

    BenthicPoolPON(BenthicReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};

	



class BenthicPoolPOP : public BenthicPool
{
    // Particulate Organic Phosphorus

    public: 

    BenthicPoolPOP(BenthicReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};

	