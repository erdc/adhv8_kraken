#pragma once

#include "NSMK_CORE.h"



struct OxygenConstraints
{
    // F I E L D S

    public:

    static const double k_oxdn;     // Exponential coeff for denitrification [L/mgO2]

    static const double k_oxmc;     // Exponential coeff for DOC mineralization [L/mgO2]

    static const double k_oxmn;     // Exponential coeff for DON mineralization [L/mgO2]

    static const double k_oxmp;     // Exponential coeff for DOP mineralization [L/mgO2]

    static const double k_oxna;     // Exponential coeff for nitrification [L/mgO2]


    static const double K_oxdn;     // DO half-saturation const for denitrification [L/mgO2]

    static const double K_oxmc;     // DO half-saturation const for DOC mineralization [L/mgO2]

    static const double K_oxmn;     // DO half-saturation const for DON mineralization [L/mgO2]

    static const double K_oxmp;     // DO half-saturation const for DOP mineralization [L/mgO2]

    static const double K_oxna;     // DO half-saturation const for nitrification [L/mgO2]




    double Foxmn;      // DON mineralization attenuation due to low oxygen [dimensionless]

    double Foxna;      // Nitrification attenuation due to low oxygen [dimensionless]

    double Foxdn;      // Denitrification attenuation due to low oxygen [dimensionless]

    double Foxmp;      // DOP mineralization attenuation due to low oxygen [dimensionless]

    double Foxmc;      // DOC mineralization attenuation due to low oxygen [dimensionless]

    double Foxb;       // Bottom algae respiration attenuation due to low Oxygen conc [dimensionless]
    
    double Foxc;       //

    double Foxp;       // Attenuation due to low oxygen [dimensionless]


    // C T O R S

    OxygenConstraints();


    // M E T H O D S

    void ComputeAttenuationFactors(enum NSMKEnumOxidationAttenuationOptions, double dissolvedOxygenConc);

};

