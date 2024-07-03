#ifndef TPZGFEMCOMPELH1_H
#define TPZGFEMCOMPELH1_H

#include "TPZCompElH1.h"
#include "TPZGFemCompMesh.h"

template<class TSHAPE>
class TPZGFemCompElH1 : public TPZCompElH1<TSHAPE> {
public:
    // Constructor: use the same parameters as TPZCompElH1 and pass them to the TPZCompElH1 constructor
    TPZGFemCompElH1(TPZCompMesh &mesh, TPZGeoEl *gel);

    /// computes the shape functions in the master element AND its derivatives
    void ComputeShape(TPZVec<REAL> &intpoint, TPZMaterialData &data) override;

    void ComputeRequiredData(TPZMaterialDataT<STATE> &data,
									 TPZVec<REAL> &qsi) override;
    // Add any new methods specific to TPZGFemCompElH1
    void SetColor(GFemcolors color) {
        fColor = color;
    }

    GFemcolors Color() const {
        return fColor;
    }

    // Initialize the active connects data structure
    void IdentifyActiveConnects();
protected:

    // indexes of the active connects associated with the elements
    TPZManVector<int> fActiveConnectIndexes;

    // number of H1 functions associated with the element
    int64_t fNumH1Functions;
    // color associated with the element
    GFemcolors fColor;

private:
    // Add any new member variables specific to TPZGFemCompElH1
};

#endif // TPZGFEMCOMPELH1_H