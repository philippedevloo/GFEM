#ifndef TPZGFEMENRICHMENT_H
#define TPZGFEMENRICHMENT_H

#include "pzvec.h"
#include "pzmanvector.h"
#include "TPZSavable.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzfmatrix.h"
#include <string>
#include <functional>

#include "TPZBuildSBFem.h"
#include <nlohmann/json.hpp>

using json = nlohmann::json;


/**
 * @brief Class for computing enrichment functions for GFEM (Generalized Finite Element Method)
 * @ingroup analysis
 */
class TPZGFemEnrichment {

public:
// classifies sectors are upper, lower or neither
    enum MSectorFamily {ENone, EUpper, ELower, EAny};

private:
    /// Radius of influence for enrichment
    REAL fRadius = 0.;
    
    /// Geometric mesh
    TPZAutoPointer<TPZGeoMesh> fGeoMesh = 0;
    
    /// Computational mesh
    TPZCompMesh *fCompMesh = 0;

    /// the SBFem element group
    TPZSBFemElementGroup *fElementGroup = 0;

    /// Which part of the eigenvector to use
    enum EPartOfEigenvector {ENonInitialized,EOnlyReal, ERealPart, EImagPart};
    /// Indexes of the eigenvectors associated with the enrichment
    TPZVec<int> fEigenvectorIndexes;
    /// Part of the eigenvector to use
    TPZVec<EPartOfEigenvector> fPartOfEigenvector;

    /// @brief  Material id associated with each sector
    /// the number of sectors is equal to the number of angles-1
    TPZVec<int> fMaterialIds;

    /// @brief Arc material id for the skeleton elements
    int fArcMatId = -1;

    /// @brief Skeleton material id for the skeleton elements
    int fSkeletonMatId = 4;

    /// @brief Point material id for the crack tip elements
    int fPointMatId = 5;


    /// @brief material id offset for distinguishing the triangular elements from the collapsed elements
    int fMaterialIdOffset = 1000;

    /// @brief vector of angles
    /// the angles need to be provided in ascending order
    /// the last angle has to be equal to the initial angle + 2 Pi
    TPZManVector<REAL,10> fAngles;

    /// maximum sector of upper material
    int fMaxUpper;

    /// minimum sector of lower material
    int fMinLower;

    /// vector of skeleton elements
    TPZVec<TPZInterpolatedElement *> fSkeletonElements;

    int Sector(REAL angle, MSectorFamily sector) const {
        int nangles = fAngles.size();
        auto it = std::upper_bound(&fAngles[0], &fAngles[0]+nangles, angle);
        int pos = it-&fAngles[0]-1;
        if(sector == ELower && angle <= MinLowerAngle()) {
            pos = nangles-2;
        }
        if(sector == EUpper && angle >= MaxUpperAngle()) {
            pos = 0;
        }
        return pos;
    }
    
    REAL Eta(int sector, REAL angle) const {
#ifdef PZDEBUG
        int nangles = fAngles.size();
        if(sector < 0 || sector >= nangles-1) {
            DebugStop();
        }
#endif
        REAL angle1 = fAngles[sector];
        REAL angle2 = fAngles[(sector+1)];
        return 2.*(angle - angle1)/(angle2 - angle1)-1.;
    }

    REAL R(const std::vector<REAL> &x) const {
        REAL r = std::sqrt(x[0]*x[0]+x[1]*x[1]);
        return r;
    }

    REAL Ksi(REAL r) const {
        return r/fRadius;
    }

    void gradX(REAL ksi, REAL eta, int sector, TPZFMatrix<REAL> &gradxieta) const {
        REAL dangle = (fAngles[sector+1]-fAngles[sector]);
        REAL angle = fAngles[sector] + (eta+1.)*dangle/2.;
        REAL cosa = std::cos(angle);
        REAL sina = std::sin(angle);
        REAL r = ksi*fRadius;
        gradxieta(0,0) = fRadius*cosa/2.;
        gradxieta(1,0) = fRadius*sina/2.;
        gradxieta(0,1) = -r*sina*dangle/2.;
        gradxieta(1,1) = r*cosa*dangle/2.;
    }
    
public:

    typedef std::function<void(TPZGFemEnrichment::MSectorFamily family, const std::vector<REAL> &x, 
                                   std::vector<std::array<REAL, 2>> &phi, 
                                   std::vector<std::array<REAL ,4>> &dphixy)> GFemEnrichmentFunction;
    /// Default constructor
    TPZGFemEnrichment() = default;
    
    /// Copy constructor
    TPZGFemEnrichment(const TPZGFemEnrichment &copy) = default;
    
    /// Assignment operator
    TPZGFemEnrichment &operator=(const TPZGFemEnrichment &copy) = default;
    
    /// Destructor
    virtual ~TPZGFemEnrichment();
    
    /// Load configuration from JSON file
    void LoadFromJSON(json &input);
    
    
    /// Compute enrichment function and its derivatives
    void ComputeEnrichmentWithDerivatives(TPZGFemEnrichment::MSectorFamily family, const std::vector<REAL> &x, 
                                          std::vector<std::array<REAL, 2>> &phi, 
                                          std::vector<std::array<REAL ,4>> &dphi);
    
    /// Get enrichment function as std::function
    GFemEnrichmentFunction GetEnrichmentFunction() {
        return [this](TPZGFemEnrichment::MSectorFamily family, const std::vector<REAL> &x, 
                                   std::vector<std::array<REAL, 2>> &phi, 
                                   std::vector<std::array<REAL ,4>> &dphixy) {
            this->ComputeEnrichmentWithDerivatives(family, x, phi, dphixy);
        };
    }

    /// The enrichment function that will be returned by GetEnrichmentFunction
    void EnrichmentFunction(MSectorFamily family, const std::vector<REAL> &x, 
                                   std::vector<std::array<REAL, 2>> &phi, 
                                   std::vector<std::array<REAL ,4>> &dphixy);
    /// Set radius of influence
    void SetRadius(REAL radius) { fRadius = radius; }
    
    /// Get radius of influence
    REAL Radius() const { return fRadius; }
    
    REAL MaxUpperAngle() const {
        return fAngles[fMaxUpper+1];
    }

    REAL MinLowerAngle() const {
        return fAngles[fMinLower];
    }
    /// Get geometric mesh
    TPZGeoMesh *GetGeoMesh() const { return fGeoMesh.operator->(); }
    
    /// Get computational mesh
    TPZCompMesh *GetCompMesh() const { return fCompMesh; }

    /// Plot the enrichment function
    void PlotEnrichmentFunction(const std::string &filename, json &input);
    
    /**
     * Rotation utilities (made public as requested).
     * These are stateless static helpers and safe to expose.
     */
    static void ComputeRotationMatrix(REAL angle, TPZFNMatrix<4,REAL> &rotmat);
    static void RotateVector(REAL angle, std::array<REAL,2> &vec);
    static void RotateTensor(REAL angle, std::array<REAL,4> &tensor);
    
private:

    /// @brief  read the geometric mesh for the file name defined in the json input
    /// @param input json object with the input data
    /// @return pointer to the geometric mesh
    TPZGeoMesh *CreateGeometry(json &input);


    /// @brief Create a computational mesh based on the SBFem hybrid H1 space
    /// @param gmesh autopointer to the geometric mesh
    /// @param input json object with the input data
    /// @return pointer to the computational mesh
    TPZCompMesh *CreateSBFemSpace(TPZBuildSBFem &builder, json &input);

    /// @brief configure the SBFem builder according to the json input
    /// @param builder reference to the SBFem builder
    /// @param input json object with the input data
    void ConfigureSBFemBuilder (TPZBuildSBFem &builder, json &input);
    /// @brief insert material objects into the H1 computational mesh according to the json input
    /// @param cmesh reference to the computational mesh
    /// @param input json object with the input data
    void InsertMaterialObjectsH1(TPZCompMesh &cmesh, json &input);

    /// @brief Compute the stiffness matrix 
    void ComputeSBFemApproximationSpace();

    /// @brief return the eigenvalues corresponding to singular modes
    /// store the result in fEigenvectorIndexes and fPartOfEigenvector
    void IdentifySingularModes();

    /// @brief Load an eigenmode from a file into the computational mesh
    void LoadEigenMode(int64_t modeindex);

    /// @brief Plot the solution of the computational mesh
    /// @param cmesh reference to the computational mesh
    /// @param output json object with the output data
    void PlotSolution(const std::string &rootname, TPZCompMesh &cmesh, int refstep, json &input);

    /// @brief extend the enrichment function for angles outside the provided range
    void ExtendEnrichmentFunction(MSectorFamily family,
                                 const std::vector<REAL> &x,
                                 std::vector<std::array<REAL, 2>> &phi,
                                 std::vector<std::array<REAL ,4>> &dphixy);

    /* ...existing code... */
};

#endif // TPZGFEMENRICHMENT_H