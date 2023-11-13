#pragma once

#include "NSMK_CORE.h"




class PartitionDistribution : public NSMKPartitionDistribution
{
    public:

    double* P;     // Partition multipilers [ used by Partition objs]


    public:

    PartitionDistribution(const int numberParticles = 0);

    virtual ~PartitionDistribution();


    public:

    void Clear();

    double Update();
};


