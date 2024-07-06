
#include "TPZGFemCompElH1.h"
#include "TPZMaterialDataT.h"
#include "TPZMaterialT.h"
#include "pzaxestools.h"
// Implement the TPZGFemCompElH1 class here
// this is where we implement the methods of the class

// Constructor: use the same parameters as TPZCompElH1 and pass them to the TPZCompElH1 constructor
template<class TSHAPE>
TPZGFemCompElH1<TSHAPE>::TPZGFemCompElH1(TPZCompMesh &mesh, TPZGeoEl *gel, GFemcolors color) : TPZCompElH1<TSHAPE>(mesh, gel) {
    // identify the number of H1 functions associated with the element
    fNumH1Functions = this->NShapeF();
    // initialize the color of the element
    fColor = color;
    TPZGFemCompMesh *gfemmesh = dynamic_cast<TPZGFemCompMesh *>(&mesh);
    if (!gfemmesh) {
        DebugStop();
    }
}

// Initialize the active connects data structure
template<class TSHAPE>
void TPZGFemCompElH1<TSHAPE>::IdentifyActiveConnects() {
    TPZGFemCompMesh *gfemmesh = dynamic_cast<TPZGFemCompMesh *>(this->Mesh());
    if (!gfemmesh) {
        DebugStop();
    }
    int nc = this->NConnects();
    fActiveConnectIndexes.resize(nc);
    int count = 0;
    for (int ic = 0; ic < nc; ic++) {
        int64_t cindex = this->ConnectIndex(ic);
        if(gfemmesh->fShapeFunctionMap.find(cindex) != gfemmesh->fShapeFunctionMap.end()) {
            fActiveConnectIndexes[count] = cindex;
            count++;
        }
    }
    fActiveConnectIndexes.resize(count);
}

#include "TPZGFemElasticity2D.h"

// computes the shape functions in the master element AND its derivatives
template<class TSHAPE>
void TPZGFemCompElH1<TSHAPE>::ComputeShape(TPZVec<REAL> &intpoint, TPZMaterialData &data) {
    // call the base class method to compute the shape functions
    int eldim = TSHAPE::Dimension;
    if(data.fShapeType != TPZMaterialData::EScalarShape) DebugStop();
    TPZCompElH1<TSHAPE>::ComputeShape(intpoint, data);
    if(fActiveConnectIndexes.size() == 0) {
        return;
    }
    TPZMaterialT<STATE> *mat = dynamic_cast<TPZMaterialT<STATE> *>(this->Material());
    int nstate = mat->NStateVariables();
    TPZGFemElasticity2D::AtomicToVec(nstate, data);
    TPZFMatrix<REAL> dphixyz(3*nstate,data.dphix.Cols());
    TPZGFemElasticity2D::VecAxesToXYZ(nstate, data.axes, data.dphix, dphixyz);

    if(0 && eldim ==2) {
        data.phi.Print("phi =", std::cout, EMathematicaInput);
        dphixyz.Print("dphixyz =", std::cout, EMathematicaInput);
    }
    // compute the shape functions associated with the element
    TPZGFemCompMesh *gfemmesh = dynamic_cast<TPZGFemCompMesh *>(this->Mesh());
    if (!gfemmesh) {
        DebugStop();
    }
    int dim = gfemmesh->Dimension();
    TPZFNMatrix<9,REAL> GVal(nstate,nstate);
    TPZFNMatrix<9,REAL> DValxyz(3*nstate,nstate);

    TPZFNMatrix<20,REAL> gphi = data.phi;
    TPZFNMatrix<40,REAL> gdphi = data.dphix;
    int nc = this->NConnects();
    int count = 0;
    int phicount = 0;
    for (int ic = 0; ic < nc; ic++) {
        int64_t cindex = this->ConnectIndex(ic);
        if(cindex != fActiveConnectIndexes[count]) {
            phicount += this->Connect(ic).NShape();
            continue;
        }
        TPZConnect &c = this->Connect(ic);
        int nshape = c.NShape();
        auto iter = gfemmesh->fShapeFunctionMap.find(cindex);
        if(iter == gfemmesh->fShapeFunctionMap.end()) {
            DebugStop();
        }
        iter->second(data.x, fColor, GVal, DValxyz);
        for(int i=0; i<nshape; i++) {
            REAL phikeep = data.phi(0,phicount*nstate);
            TPZFNMatrix<3,REAL> dphikeep(3,1);
            for(int j=0; j<3; j++) {
                dphikeep(j,0) = dphixyz(j,phicount*nstate);
            }
            // dphixyz.Print("dphixyz =", std::cout, EMathematicaInput);
            // dphikeep.Print("dphikeep =", std::cout, EMathematicaInput);
            for(int ist = 0; ist<nstate; ist++) {
                int funcpos = phicount*nstate+ist;
                for(int jst = 0; jst < nstate; jst++) {
                    data.phi(jst,funcpos) = GVal(jst,ist)*phikeep;
                    for(int j=0; j<3; j++) {
                        dphixyz(jst*3+j,funcpos) = DValxyz(jst*3+j,ist)*phikeep+GVal(jst,ist)*dphikeep(j,0);
                    }
                }
            }
            phicount++;
        }
        count++;
        if(count == fActiveConnectIndexes.size()) break;
    }
    if(count != fActiveConnectIndexes.size()) {
        DebugStop();
    }
    if(0 && eldim ==2) {
        data.phi.Print("phiA =", std::cout, EMathematicaInput);
        dphixyz.Print("dphixyzA =", std::cout, EMathematicaInput);
    }

    TPZGFemElasticity2D::VecXYZToAxes(nstate, data.axes, dphixyz, data.dphix);
    if(0 && eldim ==2) {
        dphixyz.Print("dphixA =", std::cout, EMathematicaInput);
    }
    if(nstate == 1) 
    {   
        data.fShapeType = TPZMaterialData::EScalarShape;
        data.phi.Transpose();
    }
    else data.fShapeType = TPZMaterialData::EVecShape;
}

template<class TSHAPE>
void TPZGFemCompElH1<TSHAPE>::ComputeRequiredData(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi) {
    TPZGFemCompMesh *gfemmesh = dynamic_cast<TPZGFemCompMesh *>(this->Mesh());
    if (!gfemmesh) {
        DebugStop();
    }
    TPZMaterialT<STATE> *mat = dynamic_cast<TPZMaterialT<STATE> *>(this->Material());
    int nstate = mat->NStateVariables();
    int dim = gfemmesh->Dimension();
    data.phi.Redim(fNumH1Functions, 1);
    data.dphix.Redim(dim, fNumH1Functions);
    data.fShapeType = TPZMaterialData::EScalarShape;
    // bring the derivatives back to the axes of the element
    TPZCompElH1<TSHAPE>::ComputeRequiredData(data, qsi);
}

#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"

using namespace pzshape;

template class TPZGFemCompElH1<TPZShapeLinear>;
template class TPZGFemCompElH1<TPZShapeQuad>;
template class TPZGFemCompElH1<TPZShapeTriang>;
template class TPZGFemCompElH1<TPZShapeTetra>;