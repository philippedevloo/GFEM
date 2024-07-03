
#include "TPZGFemCompMesh.h"
#include "TPZGFemCompElH1.h"

#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"

using namespace pzshape;

// Implement the TPZGFemCompMesh class here
// this is where we implement the methods of the class
void TPZGFemCompMesh::IdentifyActiveConnects() {
    int64_t nel = fElementVec.NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fElementVec[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        auto elt = gel->Type();
        switch(elt) {
            case EPoint:
                DebugStop();
                break;
            case EOned:
            {
                auto oned = dynamic_cast<TPZGFemCompElH1<TPZShapeLinear> *>(cel);
                if (!oned) {
                    DebugStop();
                }
                oned->IdentifyActiveConnects();
            }
                break;
            case ETriangle:
            {
                auto tri = dynamic_cast<TPZGFemCompElH1<TPZShapeTriang> *>(cel);
                if (!tri) {
                    DebugStop();
                }
                tri->IdentifyActiveConnects();
            }
                break;
            case EQuadrilateral:
            {
                auto quad = dynamic_cast<TPZGFemCompElH1<TPZShapeQuad> *>(cel);
                if (!quad) {
                    DebugStop();
                }
                quad->IdentifyActiveConnects();
            }
            case ETetraedro:
            {
                auto tet = dynamic_cast<TPZGFemCompElH1<TPZShapeTetra> *>(cel);
                if (!tet) {
                    DebugStop();
                }
                tet->IdentifyActiveConnects();
            }
                break;
            case EPrisma:
            case ECube:
            case EPiramide:
                std::cout << "Element type not implemented" << std::endl;
                DebugStop();
                break;
            default:
                DebugStop();
        }
    }
}