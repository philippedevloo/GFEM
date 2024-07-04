/**
 * @file This file implements an error estimator in space H1.
 */

#include "pzlog.h"

#include "TPZAcademicGeoMesh.h"
#include "pzinterpolationspace.h"
#include "pzcheckgeom.h"
#include "pzstepsolver.h"
#include <pzbuildmultiphysicsmesh.h>
#include "pzgeoelbc.h"
#include "pzskylstrmatrix.h"
#include "pzskylmat.h"
#include "TPZLinearAnalysis.h"
#include "TPZBndCondT.h"
#include <TPZNullMaterial.h>
#include <TPZMultiphysicsCompMesh.h>
#include "TPZSSpStructMatrix.h"
#include "TPZSYSMPMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "TPZGmshReader.h"
#include "TPZAnalyticSolution.h"
#include <TPZSimpleTimer.h>
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "Projection/TPZL2ProjectionCS.h"
#include "tpzpoint.h"
#include "pzgeopoint.h"
#include "TPZElementMatrixT.h"              

#include <fstream>
#include <ctime>
#include <cstdio>
#include <cmath>

#include "TPZGenGrid2D.h"

#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"

#include "TPZGFemCompMesh.h"

#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include <iostream>

using namespace pzgeom;
/// @brief Read a gmsh file and return a geometric mesh
/// @param filename 
/// @return a geometric mesh
TPZGeoMesh *ReadGmsh(const std::string &filename);

/// @brief Include the missing boundary elements
/// @param gmesh 
void AdjustGeoMesh(TPZGeoMesh *gmesh);

// @brief return the coordinates of the nodes at the bottom of the domain
TPZFMatrix<REAL> BottomRadiusNodes(TPZGeoMesh *gmesh, TPZVec<REAL> &angles);

/// @brief Create a computational mesh with H1 elements
TPZCompMesh *CreateH1CompMesh(TPZGeoMesh *gmesh);

/// @brief Create a computational mesh with GFem elements
TPZGFemCompMesh *CreateGFemCompMesh(TPZGeoMesh *gmesh);

/// @brief Create the computational "multiphysics" mesh with only HDiv elements
TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(TPZCompMesh *cmeshH1, TPZGFemCompMesh *cmeshGFem);

/// @brief Simulate the NACA profile using H1 approximation
void Simulate(TPZCompMesh *cmesh, std::string name);

/// @brief print the results of the analysis
void PrintResults(TPZCompMesh *cmesh, std::string plotfile);

/// @brief Clean up the multiphysics mesh
void CleanUp(TPZMultiphysicsCompMesh *cmesh_m);

/// @brief Insert material objects into the computational mesh
void InsertMaterialObjectsDarcy(TPZCompMesh *cmesh);

/// @brief Insert material objects into the computational mesh
void InsertMaterialObjectsElasticity3D(TPZCompMesh *cmesh);

/// @brief Insert material objects into the computational mesh
void InsertMaterialObjectsElasticity3DMF(TPZMultiphysicsCompMesh *cmesh);

/// @brief Insert material objects int the multiphysics mesh
void InsertMaterialObjectsDarcyMF(TPZMultiphysicsCompMesh *cmesh_m);

void exportData(TPZCompMesh *cmesh, int64_t id);

int volmat = 1;//Domain
int BCt = 2;//PressIn
int BCb = 3;//PressOut
int middle = 10;//Middle
int fracedge = 12;//Fracture edge
int Hole = 15;//Hole
int noslip = 4;
int outer = 5;
int pointmat1 = -1;
int pointmat2 = -2;
int pointmat3 = -3;
enum simultype {H1,GFem,GFemNofrac};
simultype simtype = GFemNofrac;

int main() {

    gRefDBase.InitializeRefPatterns(3);
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    TPZGeoMesh *gmesh = 0;
    
    if(1) {
        gmesh = ReadGmsh("Reference.msh");
        std::ofstream out2("gmeshbefore.vtk"); 
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
        TPZManVector<REAL,3> angles;
        auto coords = BottomRadiusNodes(gmesh,angles);
        std::ofstream out3("angles.nb");
        coords.Print("coords = ",out3,EMathematicaInput);
        out3 << "angles = {" << angles << "}" << std::endl;
        //TPZGeoMesh *gmesh = ReadGmsh("quadmesh.msh");
        AdjustGeoMesh(gmesh);
    } else 
    {
        int nel = 5;
        TPZAcademicGeoMesh acadgmesh(nel,TPZAcademicGeoMesh::ETetrahedra);
        TPZVec<int> bcids(6,outer);
        bcids[0] = BCt;
        bcids[5] = BCb;
        acadgmesh.SetBCIDVector(bcids);
        gmesh = acadgmesh.CreateGeoMesh();
    }
    {
        std::ofstream out("gmesh.txt");
        gmesh->Print(out);

        int64_t nel = gmesh->NElements();
        TPZVec<REAL> area(nel,0.);
        for(int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->Dimension() == 3)
                area[el] = 1./gel->Volume();
        }
        std::ofstream out2("gmesh.vtk"); 
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2,area,true);
    }

    if(0)
    {
        TPZCheckGeom check(gmesh);
        int uniform = 0;
        check.UniformRefine(uniform);

        {
            std::ofstream out4("gmeshfine.txt");
            gmesh->Print(out4);

            std::ofstream out5("gmeshfine.vtk"); 
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out5);
        }
    }
    // indicating the flux order
    int defaultporder = 1;
    int64_t nel = gmesh->NElements();
    auto cmeshH1 = CreateH1CompMesh(gmesh);
    Simulate(cmeshH1,"H1");
    return 0;
    auto cmeshGFem = CreateGFemCompMesh(gmesh);
    auto cmesh_m = CreateMultiphysicsMesh(cmeshH1,cmeshGFem);
    Simulate(cmesh_m,"GFem");
    CleanUp(cmesh_m);
    return 0;   
}

/// @brief Read a gmsh file and return a geometric mesh
/// @param filename 
/// @return a geometric mesh

// int volmat = 1;//Domain
// int BCt = 2;//Pressin
// int BCb = 3;//Pressout
// int Hole = 200;//Hole


TPZGeoMesh *ReadGmsh(const std::string &meshfilename)
{
    TPZGmshReader gmsh;
    gmsh.SetVerbose(1);
    gmsh.GetDimNamePhysical()[2]["PressIn"] = BCt;  
    gmsh.GetDimNamePhysical()[2]["PressOut"] = middle;
    gmsh.GetDimNamePhysical()[2]["Hole"] = Hole;
    gmsh.GetDimNamePhysical()[2]["NoSlip"] = noslip;  
    gmsh.GetDimNamePhysical()[2]["NoPenetration"] = outer;
    gmsh.GetDimNamePhysical()[3]["Domain"] = volmat;
    auto gmesh = gmsh.GeometricGmshMesh(meshfilename);
    return gmesh;
}

void AdjustNodeHeight(TPZGeoMesh *gmesh){
    int64_t nelem = gmesh->NElements();
    for(int64_t el = 0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->Dimension() != 2) continue;
        int matid = gel->MaterialId();
        if(matid != outer) continue;
        if(gel->HasSubElement()) continue;
        int nnodes = gel->NNodes();
        for(int64_t in = 0; in<nnodes; in++) {
            TPZManVector<REAL,3> xco(3);
            gel->NodePtr(in)->GetCoordinates(xco);
            REAL zlarge = xco[2]/0.05;
            int level = round(zlarge);
            REAL betterh = level*0.05;
            if(fabs(level-zlarge) > 0.2) {
                std::cout << "Node " << in << " height " << xco[2] << " better height " << betterh << std::endl;
            }
            gel->NodePtr(in)->SetCoord(2,betterh);
        }
    }
}

TPZFMatrix<REAL> BottomRadiusNodes(TPZGeoMesh *gmesh, TPZVec<REAL> &anglevec) {
    std::map<REAL,int64_t> angles;
    int64_t nelem = gmesh->NElements();
    for(int64_t el = 0; el <nelem; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->Dimension() != 2) continue;
        int matid = gel->MaterialId();
        if(matid != outer) continue;
        if(gel->HasSubElement()) continue;
        for(int in = 0; in < gel->NCornerNodes(); in++) {
            TPZManVector<REAL,3> xco(3);
            gel->NodePtr(in)->GetCoordinates(xco);
            // if(fabs(xco[2]) > 1.e-6) continue;
            REAL norm = 0.;
            for(int i = 0; i<2; i++) norm += xco[i]*xco[i];
            norm = sqrt(norm);
            for(int i = 0; i<2; i++) xco[i] /= norm;
            REAL anglel = atan2(xco[1],xco[0]);
            if(angles.find(anglel) != angles.end()) continue;
            angles[anglel] = gel->NodeIndex(in);
        }
    }
    TPZFMatrix<REAL> result(angles.size(),3,0.);
    anglevec.Resize(angles.size());
    int64_t count = 0;
    for(auto it : angles) {
        TPZManVector<REAL,3> xco(3);
        gmesh->NodeVec()[it.second].GetCoordinates(xco);
        for(int i = 0; i<3; i++) result(count,i) = xco[i];
        anglevec[count] = it.first;
        count++;
    }
    return result;
}
void AdjustAngles(TPZGeoMesh *gmesh) {
    int64_t nelem = gmesh->NElements();
    for(int64_t el = 0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->Dimension() != 2) continue;
        int matid = gel->MaterialId();
        if(matid != outer) continue;
        if(gel->HasSubElement()) continue;
        for(int in = 0; in < gel->NCornerNodes(); in++) {
            TPZManVector<REAL,3> xco(3);
            gel->NodePtr(in)->GetCoordinates(xco);
            REAL norm = 0.;
            for(int i = 0; i<2; i++) norm += xco[i]*xco[i];
            norm = sqrt(norm);
            for(int i = 0; i<2; i++) xco[i] /= norm;
            REAL anglel = atan2(xco[1],xco[0]);
            xco[0] *= norm;
            xco[1] *= norm;
            REAL anglebig = 32.*anglel/M_PI;
            REAL angledegree = anglebig*180./32.;
            int anglecount = round(anglebig);
            REAL betterangle = anglecount*M_PI/32;
            if(fabs(anglebig-anglecount) > 0.25) {
                std::cout << "Node " << gel->NodeIndex(in) << " angle " << anglel << " better angle " << betterangle << " anglebig " << anglebig<< " xco " << xco << std::endl;
            }
            xco[0] = cos(betterangle)*norm;
            xco[1] = sin(betterangle)*norm;
            gel->NodePtr(in)->SetCoord(0,xco[0]);
            gel->NodePtr(in)->SetCoord(1,xco[1]);
        }
    }
}

REAL Radius(TPZGeoMesh *gmesh) {
    int64_t nelem = gmesh->NElements();
    REAL radius = 0.;
    for(int64_t el = 0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->Dimension() != 2) continue;
        int matid = gel->MaterialId();
        if(matid != outer) continue;
        TPZManVector<REAL,3> xco(3);
        gel->NodePtr(0)->GetCoordinates(xco);
        REAL norm = 0.;
        for(int i = 0; i<2; i++) norm += xco[i]*xco[i];
        norm = sqrt(norm);
        radius = norm;
    }
    return radius;
}

void AdjustRadius(TPZGeoMesh *gmesh) {
    REAL radius = Radius(gmesh);
    int64_t nelem = gmesh->NElements();
    for(int64_t el = 0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->Dimension() != 2) continue;
        int matid = gel->MaterialId();
        if(matid != outer) continue;
        for(int in = 0; in < gel->NCornerNodes(); in++) {
            TPZManVector<REAL,3> xco(3);
            gel->NodePtr(in)->GetCoordinates(xco);
            REAL norm = 0.;
            for(int i = 0; i<2; i++) norm += xco[i]*xco[i];
            norm = sqrt(norm);
            for(int i = 0; i<2; i++) xco[i] *= radius/norm;
            gel->NodePtr(in)->SetCoord(0,xco[0]);
            gel->NodePtr(in)->SetCoord(1,xco[1]);
        }
    }
}
bool HasVertical(TPZGeoEl *gel) {
    for(int in = 0; in < gel->NCornerNodes(); in++) {
        TPZManVector<REAL,3> xco(3),xco2(3);
        gel->NodePtr(in)->GetCoordinates(xco);
        gel->NodePtr((in+1)%gel->NCornerNodes())->GetCoordinates(xco2);
        REAL norm = 0.;
        for(int i = 0; i<2; i++) norm += (xco[i]-xco2[i])*(xco[i]-xco2[i]);
        norm = sqrt(norm);
        if(norm < 1.e-8) return true;
    }
    return false;
}

int HorizontalSide(TPZGeoEl *gel) {
    for(int is = gel->NCornerNodes(); is<gel->NSides()-1; is++) {
        TPZManVector<REAL,3> xco1(3),xco2(3);
        int s1 = gel->SideNodeLocIndex(is,0);
        int s2 = gel->SideNodeLocIndex(is,1);
        gel->NodePtr(s1)->GetCoordinates(xco1);
        gel->NodePtr(s2)->GetCoordinates(xco2);
        if(fabs(xco1[2]-xco2[2]) < 1.e-6) return is;
    }
    return -1;

}

void DivideNecessaryElements(TPZGeoMesh *gmesh) {
    int64_t nelem = gmesh->NElements();
    for(int64_t el = 0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->Dimension() != 2) continue;
        int matid = gel->MaterialId();
        if(matid != outer) continue;
        if(gel->HasSubElement()) continue;
        if(HasVertical(gel)) continue;
        int horizside = HorizontalSide(gel);
        if(horizside == -1) DebugStop();
        {
            TPZStack<TPZGeoElSide> neighbours;
            TPZGeoElSide gelside(gel,horizside);
            neighbours.push_back(gelside);
            gelside.AllNeighbours(neighbours);
            for(auto neighbour : neighbours) {
                TPZGeoEl *neighgel = neighbour.Element();
                if(neighgel->HasSubElement()) continue;
                TPZManVector<int,7> sidestorefine(neighgel->NSides(),0);
                sidestorefine[neighbour.Side()] = 1;
                auto refpat = TPZRefPatternTools::PerfectMatchRefPattern(neighgel,sidestorefine);
                neighbour.Element()->SetRefPattern(refpat);
                TPZStack<TPZGeoEl *> sons;
                neighbour.Element()->Divide(sons);
            }
        }
    }
    AdjustAngles(gmesh);
    AdjustRadius(gmesh);

}
void ProjectNodes(TPZGeoMesh *gmesh) {
    AdjustNodeHeight(gmesh);
    AdjustAngles(gmesh);
    AdjustRadius(gmesh);
    DivideNecessaryElements(gmesh);
    AdjustAngles(gmesh);
    AdjustRadius(gmesh);
    int64_t nelem = gmesh->NElements();
    for(int64_t el = 0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->Dimension() != 2) continue;
        int matid = gel->MaterialId();
        if(matid != outer) continue;
        if(gel->HasSubElement()) continue;
        TPZManVector<REAL,3> normal(3,0.);
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZManVector<REAL,3> xi(2,0.);
        gelside.Normal(xi, normal);
        std::cout << "Normal " << normal << std::endl;
    }

}

/// @brief Include the missing boundary elements
/// @param gmesh
void AdjustGeoMesh(TPZGeoMesh *gmesh)
{
    using namespace pzgeom;
    TPZVec<int64_t> nodeindex = {0};
    TPZGeoElRefPattern<TPZGeoPoint> *p1 = new TPZGeoElRefPattern<TPZGeoPoint>(0, nodeindex, pointmat1, *gmesh);
    nodeindex[0] = 1;
    TPZGeoElRefPattern<TPZGeoPoint> *p2 = new TPZGeoElRefPattern<TPZGeoPoint>(1, nodeindex, pointmat2, *gmesh);
    nodeindex[0] = 2;
    TPZGeoElRefPattern<TPZGeoPoint> *p3 = new TPZGeoElRefPattern<TPZGeoPoint>(2, nodeindex, pointmat3, *gmesh);
    gmesh->BuildConnectivity();

    int numouter = 0;
    int numedge = 0;
    int dim = gmesh->Dimension();
    std::set<int> matids = {BCt,outer};
    int64_t nel = gmesh->NElements();
    for(int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->MaterialId() == volmat) {
            int firstside = gel->FirstSide(2);
            for(int is = firstside; is<gel->NSides()-1; is++) {
                TPZGeoElSide gelside(gel,is);
                auto res = gelside.Neighbour().HasNeighbour(matids);
                if(res) continue;
                res = gelside.Neighbour().HasNeighbour(volmat);
                // found a boundary element
                if(res == gelside){
                    TPZGeoElBC(gelside,BCb);
                    numouter++;
                }
            }
        }
        if(gel->MaterialId() == Hole) {
            int firstside = gel->FirstSide(1);
            for(int is = firstside; is<gel->NSides()-1; is++) {
                TPZGeoElSide gelside(gel,is);
                if(gelside.Neighbour().HasNeighbour(Hole) == gelside) {
                    TPZGeoElBC(gelside,fracedge);
                    numedge++;
                }
            }
        }
    }
    ProjectNodes(gmesh);

    std::cout << "Created " << numouter << " bottom boundary elements\n";
    std::cout << "Created " << numedge << " edge boundary elements\n";
}

void BuildBlueRedElements(TPZGeoMesh *gmesh, std::set<int64_t> &blue, std::set<int64_t> &red) {
    int64_t nel = gmesh->NElements();
    int dim = gmesh->Dimension();
    std::set<TPZGeoEl *> cutelements;
    for(int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->MaterialId() == Hole) {
            cutelements.insert(gel);
            TPZGeoElSide gelside(gel);
            TPZGeoElSide neighbour = gelside.Neighbour();
            int nfound = 0;
            // look at the neighbours of the cut elements
            // compute the parameter transformation. If the det is positive
            while(neighbour != gelside)
            {
                int64_t neighindex = neighbour.Element()->Index();
                if(neighbour.Element()->MaterialId() == volmat)
                {
                    if(neighbour.IsNeighbourCounterClockWise(gelside)) {
                        blue.insert(neighindex);
                    }
                    else {
                        red.insert(neighindex);
                    }
                    // TPZTransform<REAL> tr(dim-1);
                    // gelside.SideTransform3(neighbour, tr);
                    // REAL det = tr.Mult()(0,0);
                    // if(dim == 3) {
                    //     auto &a = tr.Mult();
                    //     det = a(0,0)*a(1,1)-a(0,1)*a(1,0);
                    // }
                    // if(det > 0.) {
                    //     blue.insert(neighindex);
                    // } else {
                    //     red.insert(neighindex);
                    // }
                    nfound++;
                }
                neighbour = neighbour.Neighbour();
            }
            if(nfound != 2) DebugStop();
        }
    }
    // loop over the dim-1 sides of the cut elements
    for(auto gelcut : cutelements) {
        for(int side = gelcut->FirstSide(dim-2); side < gelcut->FirstSide(dim-1); side++) {
            TPZGeoElSide gelside(gelcut,side);
            if(gelside.HasNeighbour(fracedge)) continue;
            bool change = true;
            while(change) {
                change = false;
                // loop over the neighbours of the dim-1 sides of the gelcut element
                for(auto neighbour = gelside.Neighbour(); neighbour != gelside; neighbour++) {
                    int64_t neighindex = neighbour.Element()->Index();
                    TPZGeoEl *neiggel = neighbour.Element();
                    int neighdim = neiggel->Dimension();
                    // skip point elements
                    if(neighdim == 0) continue;
                    if(blue.find(neighindex) != blue.end() || red.find(neighindex) != red.end()) {
                        continue;
                    }
                    // loop over the face sides of the neighbour element. 
                    // if there is a neighbour of a colour, adopt that colour
                    bool found = false;
                    int lastside;
                    if(neighdim == 1) lastside = neiggel->NSides();
                    else lastside = neiggel->FirstSide(dim);
                    for(int faceside = neighbour.Element()->FirstSide(1); faceside < lastside; faceside++) {
                        TPZGeoElSide face(neighbour.Element(),faceside);
                        for(auto neighbourface = face.Neighbour(); neighbourface != face; neighbourface++) {
                            int64_t neighfaceindex = neighbourface.Element()->Index();
                            if(blue.find(neighfaceindex) != blue.end()) {
                                blue.insert(neighindex);
                                found = true;
                                change = true;
                            }
                            if(red.find(neighfaceindex) != red.end()) {
                                red.insert(neighindex);
                                found = true;
                                change=true;
                            }
                            if(found) break;
                        }
                        if(found) break;
                    }
                }
                // if no change was made, the all neighours should be either blue or red
                if(change == false) {
                    for(auto neighbour = gelside.Neighbour(); neighbour != gelside; neighbour++) {
                        int64_t neighindex = neighbour.Element()->Index();
                        if(blue.find(neighindex) == blue.end() && red.find(neighindex) == red.end()) {
                            DebugStop();
                        }
                    }
                }
            }
        }
    }
        // loop over the dim-1 sides of the cut elements
    for(auto gelcut : cutelements) {
        for(int side = 0; side < gelcut->FirstSide(dim-2); side++) {
            TPZGeoElSide gelside(gelcut,side);
            if(gelside.HasNeighbour(fracedge)) continue;
            bool change = true;
            while(change) {
                change = false;
                // loop over the neighbours of the dim-1 sides of the gelcut element
                for(auto neighbour = gelside.Neighbour(); neighbour != gelside; neighbour++) {
                    int64_t neighindex = neighbour.Element()->Index();
                    TPZGeoEl *neiggel = neighbour.Element();
                    int neighdim = neiggel->Dimension();
                    // skip point elements
                    if(neighdim == 0) continue;
                    if(blue.find(neighindex) != blue.end() || red.find(neighindex) != red.end()) {
                        continue;
                    }
                    // loop over the face sides of the neighbour element. 
                    // if there is a neighbour of a colour, adopt that colour
                    bool found = false;
                    int lastside;
                    if(neighdim == 1) lastside = neiggel->NSides();
                    else lastside = neiggel->FirstSide(dim);
                    for(int faceside = neighbour.Element()->FirstSide(1); faceside < lastside; faceside++) {
                        TPZGeoElSide face(neighbour.Element(),faceside);
                        for(auto neighbourface = face.Neighbour(); neighbourface != face; neighbourface++) {
                            int64_t neighfaceindex = neighbourface.Element()->Index();
                            if(blue.find(neighfaceindex) != blue.end()) {
                                blue.insert(neighindex);
                                found = true;
                                change = true;
                            }
                            if(red.find(neighfaceindex) != red.end()) {
                                red.insert(neighindex);
                                found = true;
                                change=true;
                            }
                            if(found) break;
                        }
                        if(found) break;
                    }
                }
                // if no change was made, the all neighours should be either blue or red
                if(change == false) {
                    for(auto neighbour = gelside.Neighbour(); neighbour != gelside; neighbour++) {
                        int64_t neighindex = neighbour.Element()->Index();
                        if(blue.find(neighindex) == blue.end() && red.find(neighindex) == red.end()) {
                            DebugStop();
                        }
                    }
                }
            }
        }
    }

    // print the mesh with the colouring
    {
        TPZVec<int> color(gmesh->NElements(),0);
        for(auto el : blue) {
            if(!gmesh->Element(el)->HasSubElement()) color[el] = 1;
        }
        for(auto el : red) {
            if(!gmesh->Element(el)->HasSubElement()) color[el] = -1;
        }
        std::ofstream out("color.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, color);
    }
}
/// @brief Create a computational mesh with H1 elements
TPZCompMesh *CreateH1CompMesh(TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    switch(simtype)
    {
        case H1:
            InsertMaterialObjectsDarcy(cmesh);
            break;
        case GFem:
        case GFemNofrac:
            InsertMaterialObjectsElasticity3D(cmesh);
            break;
    }

    // cmesh->AutoBuild(matidsh1);
    cmesh->SetDefaultOrder(1);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    std::set<int> matids = {volmat,BCb,BCt,pointmat1,pointmat2,pointmat3};
    {
        int64_t nel = gmesh->NElements();
        TPZStack<int64_t > gelstack;
        for(int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->HasSubElement()) continue;
            int matid = gel->MaterialId();
            if(matids.find(matid) != matids.end()) {
                gelstack.Push(gel->Index());
            }
        }
        cmesh->ApproxSpace().BuildMesh(*cmesh, gelstack);
    }
    cmesh->SetDefaultOrder(1);
    cmesh->ExpandSolution();
    cmesh->ComputeNodElCon();
    {
        std::fstream out("cmesh.txt");
        cmesh->Print(out);
    }
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

#include "TPZGFemCompElH1.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"

/// @brief Create a computational mesh with L2 elements
 TPZGFemCompMesh *CreateGFemCompMesh(TPZGeoMesh *gmesh)
 {
    gmesh->ResetReference();
    TPZGFemCompMesh *cmesh = new TPZGFemCompMesh(gmesh);
    cmesh->SetDefaultOrder(1);

    cmesh->SetAllCreateFunctionsContinuous();
    switch(simtype)
    {
        case H1:
            InsertMaterialObjectsDarcy(cmesh);
            break;
        case GFem:
            InsertMaterialObjectsElasticity3D(cmesh);
            break;
        case GFemNofrac:
            break;
    }
     //o método ApproxSpace().CreateDisconnectedElements(true) é chamado para transformar os elementos em elementos discontinuos, criando assim um espaço L2 de aproximação.

    std::set<int64_t> blue, red;
    BuildBlueRedElements(gmesh, blue, red);
    for(auto el : blue) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->HasSubElement()) continue;
        int matid = gel->MaterialId();
        if(!cmesh->FindMaterial(matid)) continue;
        switch(gel->Type())
        {
            case EOned:
            {
                auto *cel = new TPZGFemCompElH1< pzshape::TPZShapeLinear >(*cmesh, gel);
                cel->SetColor(Black);
            }
            break;
            case ETriangle:
            {
                auto *cel = new TPZGFemCompElH1< pzshape::TPZShapeTriang >(*cmesh, gel);
                cel->SetColor(Black);
            }
            break;
            case EQuadrilateral:
            {
                auto *cel = new TPZGFemCompElH1< pzshape::TPZShapeQuad >(*cmesh, gel);
                cel->SetColor(Black);
            }
            break;
            case ETetraedro:
            {
                auto *cel = new TPZGFemCompElH1< pzshape::TPZShapeTetra >(*cmesh, gel);
                cel->SetColor(Black);
            }
            break;
            default:
                DebugStop();
        }
    }
    // identify the connect that need to be active
    {
        int64_t nel = gmesh->NElements();
        for(int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(gel->HasSubElement()) continue;
            if(gel->MaterialId() == Hole) {
                int nsides = gel->NSides();
                for(int is = 0; is<nsides; is++) {
                    TPZGeoElSide gelside(gel,is);
                    if(gelside.HasNeighbour(fracedge)) continue;
                    TPZGeoElSide neighbour = gelside.Neighbour();
                    for(;neighbour  != gelside; neighbour++)
                    {
                        TPZGeoEl *neighgel = neighbour.Element();
                        if(!neighgel->Reference()) continue;
                        TPZCompEl *cel = neighbour.Element()->Reference();
                        if(!cel) DebugStop();
                        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                        int64_t cindex = intel->ConnectIndex(neighbour.Side());
                        if(cmesh->fShapeFunctionMap.find(cindex) == cmesh->fShapeFunctionMap.end())
                        {
                            cmesh->fShapeFunctionMap[cindex] = BlackWhite;
                        }
                    }
                }
            }
        }
    }
    // set all other connects to inactive
    {
        int64_t nconnects = cmesh->NConnects();
        for(int64_t ic=0; ic<nconnects; ic++)
        {
            if(cmesh->fShapeFunctionMap.find(ic) == cmesh->fShapeFunctionMap.end())
            {
                cmesh->ConnectVec()[ic].SetCondensed(true);
            }
        }
    }
    cmesh->ExpandSolution();
    cmesh->IdentifyActiveConnects();
    {
        std::ofstream out("cmeshGFem.txt");
        cmesh->Print(out);
        std::ofstream out2("cmeshGFem.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out2);
    }
    return cmesh;
 }

#include "TPZGFemDarcyFlow.h"

/// @brief Create a computational "multiphysics" mesh with only HDiv elements
 TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(TPZCompMesh* cmeshH1,TPZGFemCompMesh* cmeshGFem){
//
    auto gmesh = cmeshH1->Reference();
    TPZMultiphysicsCompMesh *cmesh_m = new TPZMultiphysicsCompMesh(gmesh);

    const int dim = cmeshH1->Dimension();
    const int pord = cmeshH1->GetDefaultOrder();
    cmesh_m->SetDimModel(dim);
    cmesh_m->SetDefaultOrder(pord);
         
    cmesh_m->SetAllCreateFunctionsMultiphysicElem();
    switch (simtype)
    {
        case H1:
            InsertMaterialObjectsDarcyMF(cmesh_m);
            break;
        case GFem:
        case GFemNofrac:
            InsertMaterialObjectsElasticity3DMF(cmesh_m);
            break;
    }

    TPZManVector<int, 2> active_approx_spaces(2, 1);
    TPZManVector<TPZCompMesh *, 2> mesh_vec(2);
    mesh_vec[0] = cmeshH1;
    mesh_vec[1] = cmeshGFem;

    cmesh_m->BuildMultiphysicsSpace(active_approx_spaces, mesh_vec);

    cmesh_m->ExpandSolution();
    return cmesh_m;
 }

/// @brief Simulate the NACA profile using H1 approximation
void Simulate(TPZCompMesh *cmesh, std::string name)
{

    TPZLinearAnalysis an(cmesh,RenumType::EMetis);
    TPZSkylineStructMatrix<STATE> strmat(cmesh);
    // TPZSSpStructMatrix<STATE> strmat(cmesh);
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Assemble();
    if (1)
    {
        int nel = cmesh->NElements();
        TPZFMatrix<STATE> &sol = cmesh->Solution();
        for (int iel = 0; iel < nel; iel++)
        {
            TPZCompEl *cel = cmesh->ElementVec()[iel];
            if (!cel)
                continue;
            TPZGeoEl *gel = cel->Reference();
            if (!gel)
                continue;
            for(int in = 0; in < gel->NCornerNodes(); in++) {
                TPZManVector<REAL,3> xco(3);
                gel->NodePtr(in)->GetCoordinates(xco);
                TPZConnect &c = cel->Connect(in);
                int64_t seqnum = c.SequenceNumber();
                int64_t pos = cmesh->Block().Position(seqnum);
                sol(pos+2,0) = xco[2];
            }
        }
        std::ofstream out("cmeshFixed.txt");
        cmesh->Print(out);
        PrintResults(cmesh,"Imposed");
        TPZMatrixSolver<STATE> *solve = dynamic_cast<TPZMatrixSolver<STATE> *>(an.Solver());
        TPZAutoPointer<TPZMatrix<STATE>> matptr = solve->Matrix();
        TPZFMatrix<STATE> &rhs = an.Rhs();
        TPZFMatrix<STATE> &meshsol = cmesh->Solution();
        TPZFMatrix<STATE> loc(rhs.Rows(),1,0.);
        matptr->Multiply(meshsol,loc);
        loc-=rhs;
        an.LoadSolution(loc);
        PrintResults(cmesh,"Residual");
        std::ofstream out2("residual.txt");
        std::ofstream out3("rhs.txt");
        std::ofstream out4("mat.txt");
        std::ofstream out5("sol.txt");
        rhs.Print("rhs",out3,EMathematicaInput);
        meshsol.Print("sol",out5,EMathematicaInput);
        matptr->Print("mat",out4,EMathematicaInput);
        an.PrintVectorByElement(out2,loc, 1.0e-10);
        int nvar = an.Solution().Rows();

        exportData(cmesh,7);
    }

//  Solve the system of equations and save the solutions: phi_0 e phi_1 
    an.Solve();
    if (1)
    {
        std::ofstream out("cmeshH1.txt");
        cmesh->Print(out);
    }
    std::cout << "--------- PostProcess " << name << " ---------" << std::endl;
    //printa na tela "--------- PostProcess ---------", indicando que a simulação está em processamento.
    TPZMultiphysicsCompMesh *cmesh_m = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
    PrintResults(cmesh,name);
    //chama a função PrintResults para realizar o pós-processamento dos resultados. Essa função provavelmente gera saídas com os resultados da simulação.

}

void PrintResults(TPZCompMesh *cmesh, std::string plotfile)
//Define a função do tipo void chamada PrintResults, que recebe como parâmetros TPZLinearAnalysis &an e  TPZCompMesh *cmesh
{
 
    std::cout << "--------- Post Process ---------" << std::endl;
    //printa para o usuário "--------- Post Process ---------" indicando o início da fase de pós-processamento. 
    TPZSimpleTimer postProc("Post processing time");
    //declara uma variável chamada postProc, do tipo TPZSimpleTimer, chamando um construtor com uma string como argumento, igual a "Post processing time".
    //define o nome base do arquivo de saída para o pós-processamento. O nome base é "postprocess".
    constexpr int vtkRes{0};
    //define a variável do tipo inteiro denominada vtkRes, do tipo constexpr, que significa que é uma expressão constante, ou seja,  vtkRes é um valor constante e não pode ser alterado. Ainda, {0} indica o valor associado a essa constante, e portanto não será alterado, com valor determinado na hora de compilação.
    //define a resolução para o formato de arquivo VTK. Neste caso, a resolução é definida como 0, o que geralmente significa que a resolução será automática.
    TPZStack<std::string> fields;
    TPZMaterial *mat = cmesh->FindMaterial(volmat);
    int nstate = mat->NStateVariables();
    int dim = cmesh->Dimension();
    switch(simtype)
    {
        case H1:
            plotfile += "_darcy";
            fields.Push("Pressure");
            fields.Push("Flux");
            break;
        case GFem:
        case GFemNofrac:
            plotfile += "_elast";
            {
                if(dim == 2) {
                    fields.Push("SigmaX");
                    fields.Push("SigmaY");
                    fields.Push("SigmaZ");
                    fields.Push("TauXY");
                    fields.Push("Displacement");
                } else if (dim == 3) {
                    fields.Push("StressX");
                    fields.Push("StressY");
                    fields.Push("StressZ");
                    fields.Push("Displacement");
                    fields.Push("VonMises");
                    fields.Push("Strain");
                }
            }
            break;
        default:
            DebugStop();
    }
    //nesse conjunto de linhas de código, temos que TPZVec é uma estrutura do tipo vetor que contém como argumento uma variável chamda "fields" que é uma lista de strings, que, pelo que se chamam, são relacionadas à pressão e ao fluxo.
    //cria um vetor de strings chamado fields que contém os nomes dos campos que serão pós-processados. Neste caso, os campos incluem "Pressure" (pressão) e "Flux" (fluxo). Esses campos representam propriedades do problema que desejamos visualizar após a simulação.
    std::cout << plotfile << std::endl;
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    //essa linha de código declara uma variável chamada vtk do tipo auto, o que significa que o compilador irá deduzir o tipo que ela terá a depender do que ela é igual. No caso, ela é igual a função TPZVTKGenerator, de parâmetros cmesh, fields, plotfile, vtkRes.
    //cria um objeto vtk da classe TPZVTKGenerator, que é usado para gerar arquivos VTK a partir dos dados da malha computacional cmesh. Os argumentos passados para o construtor incluem a malha computacional, os campos a serem pós-processados, o nome base do arquivo de saída (plotfile) e a resolução VTK (vtkRes).
    vtk.SetNThreads(0);
    //define o número de threads a serem usadas durante o pós-processamento. A variável global_nthread provavelmente contém o número desejado de threads.
    vtk.Do();
    //inicia o processo de geração dos arquivos VTK. Esta função gera arquivos de saída contendo informações sobre os campos especificados na malha computacional.
    std::cout << "Total time = " << postProc.ReturnTimeDouble()/1000. << " s" << std::endl;
    //imprime o tempo gasto no pós-processamento, convertido para segundos.
    
    return;
    //a função é concluída e retorna.
}

void CleanUp(TPZMultiphysicsCompMesh *cmesh_m)
{
    auto gmesh = cmesh_m->Reference();
    {
        int64_t nnod = cmesh_m->NConnects();
        for(int64_t i=0; i<nnod; i++)
        {
            TPZConnect &c = cmesh_m->ConnectVec()[i];
            if(c.HasDependency())
            {
                c.RemoveDepend();
            }
        }    
    }

    TPZVec<TPZCompMesh*> meshvec = cmesh_m->MeshVector();
    {
        int64_t nnod = meshvec[0]->NConnects();
        for(int64_t i=0; i<nnod; i++)
        {
            TPZConnect &c = meshvec[0]->ConnectVec()[i];
            if(c.HasDependency())
            {
                c.RemoveDepend();
            }
        }    
    }
    delete cmesh_m;
    delete meshvec[0];
    delete meshvec[1];
    delete gmesh;
}

/// @brief Insert material objects into the computational mesh
void InsertMaterialObjectsDarcy(TPZCompMesh *cmesh) {
    int dim = cmesh->Dimension();
    TPZDarcyFlow *material = new TPZDarcyFlow(volmat,dim);
    cmesh->InsertMaterialObject(material);

    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);

    val2[0] = -1.;
    auto bnd3 = material->CreateBC(material, BCb, 0, val1, val2);
    cmesh->InsertMaterialObject(bnd3);
    val2[0] = 1.;
    auto bnd4 = material->CreateBC(material, BCt, 0, val1, val2);
    cmesh->InsertMaterialObject(bnd4);


}

/// @brief Insert material objects int the multiphysics mesh
void InsertMaterialObjectsDarcyMF(TPZMultiphysicsCompMesh *cmesh_m){
    int dim = cmesh_m->Dimension();
    std::set<int> materialIDs;
    TPZGFemDarcyFlow *material = new TPZGFemDarcyFlow(volmat,dim);
    //  material->SetBigNumber(10e8);
     cmesh_m->InsertMaterialObject(material);
     materialIDs.insert(volmat);

     // 2. Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);

    val2[0] = 1.;
    auto bnd3 = material->CreateBC(material, BCt, 0, val1, val2);
    cmesh_m->InsertMaterialObject(bnd3);
    val2[0] = -1.;
    auto bnd4 = material->CreateBC(material, BCb, 0, val1, val2);
    cmesh_m->InsertMaterialObject(bnd4);

}

#include "Elasticity/TPZElasticity3D.h"

void disp(const TPZVec<REAL> &x, TPZVec<REAL> &disp, TPZFMatrix<REAL> &grad) {
    disp[2] = 1.;
    disp[1] = 0;
    disp[0] = 0;
    grad(0,0) = 0.;
    grad(1,0) = 0.;
    grad(2,0) = 0.;
    grad(0,1) = 0.;
    grad(1,1) = 0.;
    grad(2,1) = 1.;
    grad(0,2) = 0.;
    grad(1,2) = -1.;
    grad(2,2) = 0.;
}
/// @brief Insert material objects int the multiphysics mesh
void InsertMaterialObjectsElasticity3D(TPZCompMesh *cmesh_m){
    int dim = cmesh_m->Dimension();
    std::set<int> materialIDs;
    STATE E = 100., nu = 0.;
    TPZManVector<STATE> force(3,0.);
    TPZElasticity3D *material = new TPZElasticity3D(volmat,E,nu,force);
    force[2] = 1.;
    material->SetPostProcessingDirection(force);
    //  material->SetBigNumber(10e8);
     cmesh_m->InsertMaterialObject(material);
     materialIDs.insert(volmat);

     // 2. Boundary Conditions
    TPZFMatrix<STATE> val1(dim,dim,0.);
    TPZManVector<STATE> val2(dim,0.);

    // 2.3 Neumann traction condition on the top
    val2[2] = 100.;
    auto bnd3 = material->CreateBC(material, BCb, 1, val1, val2);
    cmesh_m->InsertMaterialObject(bnd3);
    val2.Fill(0.);
    // 2.4 Mixed condition on the bottom
    val2[2] = -100.;
    auto bnd4 = material->CreateBC(material, BCt, 1, val1, val2);
    // bnd4->SetForcingFunctionBC(disp,1);
    cmesh_m->InsertMaterialObject(bnd4);

    val1(0,0) = 1.;
    val1(1,1) = 1.;
    val1(2,2) = 1.;
    //point mats
    val2[2] = 0.;
    auto bnd5 = material->CreateBC(material, pointmat1, 0, val1, val2);
    cmesh_m->InsertMaterialObject(bnd5);

    val1(1,1) = 0.;
    auto bnd6 = material->CreateBC(material, pointmat2, 0, val1, val2);
    cmesh_m->InsertMaterialObject(bnd6);

    val1(0,0) = 0.;
    auto bnd7 = material->CreateBC(material, pointmat3, 0, val1, val2);
    cmesh_m->InsertMaterialObject(bnd7);
}

#include "TPZGFemElasticity3D.h"

/// @brief Insert material objects int the multiphysics mesh
void InsertMaterialObjectsElasticity3DMF(TPZMultiphysicsCompMesh *cmesh_m){
    int dim = cmesh_m->Dimension();
    std::set<int> materialIDs;
    STATE E = 100., nu = 0.;
    TPZManVector<STATE> force(3,0.);
    TPZGFemElasticity3D *material = new TPZGFemElasticity3D(volmat,E,nu,force);
    //  material->SetBigNumber(10e8);
     cmesh_m->InsertMaterialObject(material);
     materialIDs.insert(volmat);

     // 2. Boundary Conditions
    TPZFMatrix<STATE> val1(dim,dim,0.);
    TPZManVector<STATE> val2(dim,0.);

    // 2.3 Neumann traction condition on the top
    val2[2] = 100.;
    auto bnd3 = material->CreateBC(material, BCt, 1, val1, val2);
    cmesh_m->InsertMaterialObject(bnd3);
    val2.Fill(0.);
    // 2.4 Mixed condition on the bottom
    //val1(1,1) = 1000.;
    val2.Fill(0.);
    auto bnd4 = material->CreateBC(material, BCb, 0, val1, val2);
    cmesh_m->InsertMaterialObject(bnd4);
}

void exportData(TPZCompMesh *cmesh, int64_t id)
{
    std::ofstream out("element.nb");

    for (int64_t i = 0; i < cmesh->NElements(); i++)
    {
        TPZCompEl *cel = cmesh->ElementVec()[i];
        if (!cel)
            continue;
        
        TPZGeoEl *gel = cel->Reference();
        if (!gel)
            continue;
        int ncon = cel->NConnects();
        for (int ic = 0; ic < ncon; ic++)
        {
            TPZConnect &c = cel->Connect(ic);
            int64_t c_index = cel->ConnectIndex(ic);
            if (c_index == id)
            {
                int ncorner = gel->NCornerNodes();
                TPZFMatrix<REAL> coords(ncorner,3,0.), sol(ncorner*3,1,0.);
                for (int in = 0; in < gel->NCornerNodes(); in++)
                {
                    TPZManVector<REAL,3> xco(3);
                    gel->NodePtr(in)->GetCoordinates(xco);
                    sol(3*in+2,0) = xco[2];
                    for (int d = 0; d < 3; d++)
                    {
                        coords(in,d) = xco[d];
                    }
                }
                TPZElementMatrixT<REAL> ek, ef;
                cel->CalcStiff(ek, ef);
                out << "ic = " << ic << ";\n";
                coords.Print("coords=",out,EMathematicaInput);
                sol.Print("sol=",out,EMathematicaInput);
                ek.fMat.Print("stiff=",out,EMathematicaInput);
                ef.fMat.Print("force=",out,EMathematicaInput);
                TPZFMatrix<REAL> rhs;
                ek.fMat.Multiply(sol,rhs);
                rhs-=ef.fMat;
                rhs.Print("residual=",out,EMathematicaInput);
            }
        }
    }
}