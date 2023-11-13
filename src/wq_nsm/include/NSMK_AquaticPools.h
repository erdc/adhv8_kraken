#pragma once


#include <map>

#include "NSMK_CORE.h"

#include "NSMK_Pool.h"

#include "NSMK_AquaticPool.h"

#include "NSMK_PartitionDistribution.h"

#include "NSMK_Pathways.h"








class AquaticPoolAB1 : public AquaticPool
{
    // Algae - Bottom substrate group 1

    public: 

    double Pab;			// Ammonium preference - dynamic


    public:

    AquaticPoolAB1(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolAB2 : public AquaticPool
{
    // Algae - Bottom substrate group 2

    public: 

    double Pab;			// Ammonium preference - dynamic


    public:

    AquaticPoolAB2(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolAP1 : public AquaticPool
{
    // Algae - Phytoplankton group 1

    public: 

    double Pap;			// Ammonium preference - dynamic


    public:

    AquaticPoolAP1(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolAP2 : public AquaticPool
{
    // Algae - Phytoplankton group 2

    public: 

    double Pap;			// Ammonium preference - dynamic


    public:

    AquaticPoolAP2(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolCBOD_fast : public AquaticPool
{
    // CBOD fast

    public:

    AquaticPoolCBOD_fast(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolCBOD_slow : public AquaticPool
{
    // CBOD slow

    public: 

    AquaticPoolCBOD_slow(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolPOM : public AquaticPool
{
    // POM 

    public: 

    AquaticPoolPOM(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolDIC : public AquaticPool
{
    // Dissolved Inorganic Carbon

    public: 

    AquaticPoolDIC(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolDIP : public AquaticPool
{    
    // Dissolved InOrganic Phosphorus

    public: 

    AquaticPoolDIP(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolDOC : public AquaticPool
{
    // Dissolved Organic Carbon 

    public: 

    AquaticPoolDOC(AquaticReactor* ar);

    virtual void Configure(const PoolAttributeGroup* pag);

    virtual double Partition(const PoolAttributeGroup* pag, const NSMKSedProperties* asp);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolDON : public AquaticPool	
{
    // Dissolved Organic Nitrogen

    public:

    AquaticPoolDON(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolDOP : public AquaticPool
{
    // Dissolved Organic Phosphorus

    public: 

    AquaticPoolDOP(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolDOX : public AquaticPool
{
    // Dissolved Oxygen

    public: 

    AquaticPoolDOX(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolNH4 : public AquaticPool
{
    // Ammonium Nitrogen

    public:

    AquaticPoolNH4(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolNO3 : public AquaticPool
{
    // Nitrate Nitrogen

    public:

    AquaticPoolNO3(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolPOC : public AquaticPool
{
    // Particulate Organic Carbon 

    public: 

    AquaticPoolPOC(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolPON : public AquaticPool
{
    // Particulate Organic Nitrogen

    public: 

    AquaticPoolPON(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};





class AquaticPoolPOP : public AquaticPool
{
    // Particulate Organic Phosphorus

    public: 

    AquaticPoolPOP(AquaticReactor* ar);

    virtual double DmDtKinetic(const NSMKPoolParameters* ap);
};




