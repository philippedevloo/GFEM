
#include "TPZGFemCompMesh.h"
#include "TPZGFemCompElH1.h"

#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"

using namespace pzshape;


// Implement the TPZGFemCompMesh class here
// this is where we implement the methods of the class
void TPZGFemCompMesh::IdentifyActiveConnectsByElement() {
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
                break;
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

// Identify the color of the elements
void TPZGFemCompMesh::GetColorMap(std::map<int64_t,GFemcolors> &colormap) {
    int64_t nel = fElementVec.NElements();
    int numblack = 0;
    int numwhite = 0;
    int numgrey = 0;
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fElementVec[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        auto elt = gel->Type();
        GFemcolors color = Nocolor;
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
                color = oned->Color();
            }
                break;
            case ETriangle:
            {
                auto tri = dynamic_cast<TPZGFemCompElH1<TPZShapeTriang> *>(cel);
                if (!tri) {
                    DebugStop();
                }
                color = tri->Color();
            }
                break;
            case EQuadrilateral:
            {
                auto quad = dynamic_cast<TPZGFemCompElH1<TPZShapeQuad> *>(cel);
                if (!quad) {
                    DebugStop();
                }
                color = quad->Color();
            }
                break;
            case ETetraedro:
            {
                auto tet = dynamic_cast<TPZGFemCompElH1<TPZShapeTetra> *>(cel);
                if (!tet) {
                    DebugStop();
                }
                color = tet->Color();
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
        colormap[cel->Index()] = color;
        if(color == Black) {
            numblack++;
        } else if(color == White) {
            numwhite++;
        } else if(color == Grey) {
            numgrey++;
        } else if(color == Nocolor) {
            DebugStop();
        }
    }
    std::cout << "Number of black elements " << numblack << std::endl;
    std::cout << "Number of white elements " << numwhite << std::endl;
    std::cout << "Number of grey elements " << numgrey << std::endl;
}

#include "TPZMaterial.h"
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
// Initialize the fShapeFunctionMap data structure
void TPZGFemCompMesh::InitializeShapeFunctionMap() {
    TPZGeoMesh *gmesh = Reference();
    gmesh->ResetReference();
    LoadReferences();
    std::map<int64_t,GFemcolors> colormap;
    // get the color of each element
    GetColorMap(colormap);
    std::set<int> matids;
    for(auto it : fMaterialVec) {
        matids.insert(it.first);
    }

    // associate a value with the function
    std::map<int64_t,int> connectvalue;
    int64_t nel = fElementVec.NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fElementVec[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        int nsides = gel->NCornerNodes();
        for (int is = 0; is < nsides; is++) {
            // colors of all connected elements
            std::set<GFemcolors> colors = {colormap[cel->Index()]};
            int64_t cindex = cel->ConnectIndex(is);
            if(fShapeFunctionMap.find(cindex) != fShapeFunctionMap.end()) {
                continue;
            }
            TPZGeoElSide gelside(gel, is);
            TPZStack<TPZGeoElSide> equal;
            gelside.AllNeighbours(equal);
            int neq = equal.size();
            bool activeside = true;
            for(auto neigh : equal) {
                if(neigh.Element()->HasSubElement()) {
                    continue;
                }
                int neighmatid = neigh.Element()->MaterialId();
                if(matids.find(neighmatid) == matids.end()) {
                    continue;
                }
                TPZCompEl *neighcel = neigh.Element()->Reference();
                if(!neigh.Element()->Reference()) {
                    activeside = false;
                    break;
                }
                int64_t neighcelindex = neighcel->Index();
                if(colormap.find(neighcelindex) == colormap.end()) {
                    DebugStop();
                }
                colors.insert(colormap[neighcel->Index()]);
            }
            if(!activeside) {
                continue;
            }
            if(colors.size() == 2 && colors.find(Black) != colors.end() && colors.find(White) != colors.end()) {
                fShapeFunctionMap[cindex] = BlackWhite;
                connectvalue[cindex] = -2;
            } else {
                fShapeFunctionMap[cindex] = fFrac.DeslocFunction();
                connectvalue[cindex] = 2;
            }
        }
    }
    // add the internal connects for elements neighbouring the fracture element
    {
        extern int cutmat;
        extern int fracedge;
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (!gel) {
                continue;
            }
            if(gel->HasSubElement()) {
                continue;
            }
            if(gel->MaterialId() != cutmat) {
                continue;
            }
            int nsides = gel->NSides();
            int ncorners = gel->NCornerNodes();
            for (int is = ncorners; is < nsides; is++) {
                TPZGeoElSide gelside(gel, is);
                if(gelside.HasNeighbour(fracedge)) {
                    continue;
                }
                // colors of all connected elements
                std::set<GFemcolors> colors = {colormap[gel->Index()]};
                TPZStack<TPZGeoElSide> equal;
                gelside.AllNeighbours(equal);
                int neq = equal.size();
                for(auto neigh : equal) {
                    if(neigh.Element()->HasSubElement()) {
                        continue;
                    }
                    int neighmatid = neigh.Element()->MaterialId();
                    if(matids.find(neighmatid) == matids.end()) {
                        continue;
                    }
                    TPZCompEl *neighcel = neigh.Element()->Reference();
                    TPZConnect &c = neighcel->Connect(neigh.Side());
                    if(c.NShape()*c.NState() == 0) {
                        break;
                    }
                    int64_t neighconnectindex = neighcel->ConnectIndex(neigh.Side());
                    if(fShapeFunctionMap.find(neighconnectindex) != fShapeFunctionMap.end()) {
                        break;
                    }
                    fShapeFunctionMap[neighconnectindex] = BlackWhite;
                }
            }
        }
    }
    int nconnects = NConnects();
    for (int ic = 0; ic < nconnects; ic++) {
        TPZConnect &c = ConnectVec()[ic];
        if(fShapeFunctionMap.find(ic) == fShapeFunctionMap.end()) {
            c.SetCondensed(true);
        }
    }
    // plot the connect colors
    {
        TPZFMatrix<STATE> &sol = Solution();
        for( int ic = 0; ic < nconnects; ic++) {
            TPZConnect &c = ConnectVec()[ic];
            if(c.NShape() == 0) {
                continue;
            }
            int64_t seqnum = c.SequenceNumber();
            int64_t pos = Block().Position(seqnum);
            if(connectvalue.find(ic) == connectvalue.end()) {
                sol(pos,0) = 0.;
            } else  {
                sol(pos,0) = connectvalue[ic];
            }
        }
        TPZStack<std::string> fields;
        TPZMaterial *mat = FindMaterial(1);
        int nstate = mat->NStateVariables();
        int dim = Dimension();
        if(nstate == 1) {
            fields.Push("Pressure");
        } else if (nstate == 2) {
            fields.Push("Displacement");
        } else if (nstate == 3) {
            DebugStop();
        }
        std::string plotfile = "connectcolor";
        int vtkRes = 0;
        auto vtk = TPZVTKGenerator(this, fields, plotfile, vtkRes);
        vtk.SetNThreads(0);
        vtk.Do();
    }
    IdentifyActiveConnectsByElement();
}

/// @brief Draw the element colors
void TPZGFemCompMesh::DrawElementColors(const std::string &filename)
{
    std::map<int64_t,GFemcolors> colormap;
    // get the color of each element
    GetColorMap(colormap);
    TPZVec<REAL> elcolor(NElements(),-1);
    int64_t nel = fElementVec.NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fElementVec[el];
        if (!cel) {
            continue;
        }
        elcolor[el] = colormap[cel->Index()];
    }
    std::ofstream out(filename);
    TPZVTKGeoMesh::PrintCMeshVTK(this, out, elcolor, "elementcolor");
}
