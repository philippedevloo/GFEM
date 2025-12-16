#include "TPZGFemEnrichment.h"
#include <fstream>
#include <nlohmann/json.hpp>
#include "Elasticity/TPZLinearElasticityConstitutive.h"

typedef std::function<void(TPZGFemEnrichment::MSectorFamily family, const std::vector<REAL> &x, 
                                   std::vector<std::array<REAL, 2>> &phi, 
                                   std::vector<std::array<REAL ,4>> &dphixy)> GFemEnrichmentFunction;
int main() {
    // Instantiate TPZGFemEnrichment object
    TPZGFemEnrichment enrichment;
    
    // Open and read JSON file
    std::ifstream jsonFile("CrackDef.json");
    if (!jsonFile.is_open()) {
        std::cerr << "Error: Could not open CrackDef.json" << std::endl;
        return 1;
    }
    
    // Parse JSON
    nlohmann::json config;
    config = json::parse(jsonFile,nullptr,true,true); // to ignore comments in json file
    jsonFile.close();
    
    // Configure the enrichment object
    enrichment.LoadFromJSON(config);
    GFemEnrichmentFunction enrichFunc = enrichment.GetEnrichmentFunction();

    // enrichment.PlotEnrichmentFunction("Modes", config);

    {
        TPZGeoMesh *gmesh = enrichment.GetGeoMesh();
        std::ofstream out("gmesh_enrichment.txt");
        gmesh->Print(out);
    }

    std::vector<REAL> point = {0.5,0.5};//{M_SQRT1_2, M_SQRT1_2}; // Point (0.5, 0.5)
    {
        std::vector<std::array<REAL, 2>> phi; 
        std::vector<std::array<REAL ,4>> dphixy;
        enrichFunc(TPZGFemEnrichment::EAny, point, phi, dphixy);
        std::cout << "Enrichment function values at point (0.5, 0.5):" << std::endl;
        for (const auto& val : phi) {
            std::cout << "[" << val[0] << ", " << val[1] << "]" << std::endl;
        }
        for (const auto& deriv : dphixy) {
            std::cout << "[" << deriv[0] << ", " << deriv[1] << ", " << deriv[2] << ", " << deriv[3] << "]" << std::endl;
        }
    }
    {
        std::vector<std::array<REAL, 2>> phi; 
        std::vector<std::array<REAL ,4>> dphixy;
        enrichment.ComputeEnrichmentWithDerivatives(TPZGFemEnrichment::EAny, point, phi, dphixy);
        std::cout << "Enrichment function values at point (0.5, 0.5):" << std::endl;
        for (const auto& val : phi) {
            std::cout << "[" << val[0] << ", " << val[1] << "]" << std::endl;
        }
        for (const auto& deriv : dphixy) {
            std::cout << "[" << deriv[0] << ", " << deriv[1] << ", " << deriv[2] << ", " << deriv[3] << "]" << std::endl;
        }
    }
    return 0;
}

