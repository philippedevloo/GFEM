#ifndef TPZGFEMCOMPMESH_H
#define TPZGFEMCOMPMESH_H

#include "pzcmesh.h"

#include <map>
#include <functional>

enum GFemcolors {Black, White};
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

    // Initialize the active connects data structure
    void IdentifyActiveConnects();
};

#endif // TPZGFEMCOMPMESH_H