#ifndef TPZGFEMCOMPMESH_H
#define TPZGFEMCOMPMESH_H

#include "pzcmesh.h"
#include "TPZFrac2d.h"

#include <map>
#include <functional>

typedef std::function<void(const TPZVec<REAL> &x, GFemcolors color, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &gradphi)> GFemShapeFunctionType;

inline void BlackWhite(const TPZVec<REAL> &x, GFemcolors color, TPZFMatrix<STATE> &phi, TPZFMatrix<STATE> &gradphi)
{
    phi.Zero();
    int nr = phi.Rows();
    if (color == Black) {
        for(int i=0; i<nr; i++) phi(i,i) = 1.0;
        gradphi.Zero();
    } else {
        phi.Zero();
        gradphi.Zero();
    }

}
class TPZGFemCompMesh : public TPZCompMesh {

public:
    TPZGFemCompMesh(TPZGeoMesh *gmesh) : TPZCompMesh(gmesh) {
        
    }
    // associate a general function with a connect index
    std::map<int, GFemShapeFunctionType> fShapeFunctionMap;

    // Objects that define the two dimensional fracture
    TPZFrac2d fFrac;

    // Initialize the active connects data structure
    void IdentifyActiveConnectsByElement();

    // initialize the fShapeFunctionMap data structure
    void InitializeShapeFunctionMap();

    // Identify the color of the elements
    void GetColorMap(std::map<int64_t,GFemcolors> &colormap);

    /// @brief Draw the element colors
    void DrawElementColors(const std::string &filename);
};

#endif // TPZGFEMCOMPMESH_H