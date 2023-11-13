#pragma once



class StoichiometryRatios
{
    public:

    double Rna;     // Nitrogen/Algae ratio
    double Rpa;     // Phosphorus/Algae ratio
    double Rca;     // Carbon/Algae ratio
    double Rda;     // POM/Phytoplankton ratio
    double Roa;     // Oxygen/(Phytoplankton growth via photosynthesis) ratio
    double Ron;     // Oxygen/(Nitrogen consumed via nitrification) ratio
    double Roc;     // Oxygen/(Organic carbon oxidized to carbon dioxide) ratio
    double Rod;     // Oxygen/(POM oxidation) ratio 
    double Rondn;   // Oxygen/(Nitrogen converted via denitrification) ratio
                    //  (AKA CBOD utilized for denitrification process)
    double Rcndn;   // Carbon/(Nitrogen converted via denitrification) ratio
    double Rcco;
    double Rcca;

    StoichiometryRatios();
};






