/**
 * @file This file implements an error estimator in space H1.
 */

#include "pzlog.h"
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
//#include "TPZSYSMPMatrix.h"
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

/// @brief Refine towards the fracture edge
void RefineTowardsFrac(TPZGeoMesh *gmesh,int nref);

/// @brief Identify the fracture edge coordinates
void FractureEnds(TPZGeoMesh *gmesh, TPZVec<REAL> &fracend);

/// @brief Choose the elements to be enriched based on their distance to the fracture
/// @param gmesh geometric mesh that contains the elements
/// @param radius distance from the fracture edge
/// @param enrichedelements elements to be enriched
void ElementsToEnrich(TPZGeoMesh *gmesh, TPZVec<REAL> &fracend, REAL radius, std::set<int64_t> &enrichedelements);

/// @brief Create a computational mesh with H1 elements
TPZCompMesh *CreateH1CompMesh(TPZGeoMesh *gmesh);

/// @brief Create a computational mesh with GFem elements
TPZGFemCompMesh *CreateGFemCompMesh(TPZGeoMesh *gmesh);

/// @brief Create an H1 mirror mesh from the GFem mesh
TPZCompMesh *CreateH1MirrorMesh(TPZGFemCompMesh *cmeshGFem);

/// @brief Assemble the equations and compute the connect restraints
void ComputeConnectRestraints(TPZMultiphysicsCompMesh *cmesh);

/// @brief Create the computational "multiphysics" mesh with only HDiv elements
TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(TPZCompMesh *cmeshH1, TPZGFemCompMesh *cmeshGFem,
    TPZCompMesh *cmeshH1Mirror = nullptr);

/// @brief Simulate the NACA profile using H1 approximation
void Simulate(TPZCompMesh *cmesh, std::string name);

/// @brief print the results of the analysis
void PrintResults(TPZCompMesh *cmesh, std::string plotfile);

/// @brief Clean up the multiphysics mesh
void CleanUp(TPZMultiphysicsCompMesh *cmesh_m);

/// @brief Insert material objects into the computational mesh
void InsertMaterialObjectsDarcy(TPZCompMesh *cmesh);

/// @brief Insert material objects into the computational mesh
void InsertMaterialObjectsElasticity2D(TPZCompMesh *cmesh);

/// @brief Insert material objects into the computational mesh
void InsertMaterialObjectsElasticity2DMF(TPZMultiphysicsCompMesh *cmesh);

/// @brief Insert material objects int the multiphysics mesh
void InsertMaterialObjectsDarcy2DMF(TPZMultiphysicsCompMesh *cmesh_m);

int volmat = 1;
int BCb = 2;
int BCr = 3;
int BCt = 4;
int BCl = 5;
int cutmat = 6;
int fracedge = 7;
REAL edgelength = 0.1;
// Darcy = darcy simulation with GFem
// DarcyNoFrac = darcy simulation using multiphysics but no gfem elements
// DarcyDiscontinuous = darcy simulation with H1 discontinuous elements
enum simultype {Darcy,Elast,DarcyNofrac,ElastNoFrac,DarcyOrthogonal};
simultype simtype = DarcyOrthogonal;

int main() {
    gRefDBase.InitializeRefPatterns(3);

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    TPZGeoMesh *gmesh = ReadGmsh("SquareFrac.msh");
    //TPZGeoMesh *gmesh = ReadGmsh("quadmesh.msh");

    {
        std::ofstream out("gmesh.txt");
        gmesh->Print(out);

        std::ofstream out2("gmesh.vtk"); 
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
    }

    TPZCheckGeom check(gmesh);
    int uniform = 1;
    check.UniformRefine(uniform);

    RefineTowardsFrac(gmesh,2);

    {
        std::ofstream out4("gmeshfine.txt");
        gmesh->Print(out4);

        std::ofstream out5("gmeshfine.vtk"); 
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out5);
    }

    // indicating the flux order
    int defaultporder = 1;
    int64_t nel = gmesh->NElements();
    auto cmeshH1 = CreateH1CompMesh(gmesh);
    std::string plotfile = "H1";
    switch (simtype) {
        case Darcy:
        case DarcyNofrac:
        case DarcyOrthogonal:
            plotfile += "_Darcy";
            break;
        case Elast:
        case ElastNoFrac:
            plotfile += "_Elast";
            break;
    }
    Simulate(cmeshH1,plotfile);
    auto cmeshGFem = CreateGFemCompMesh(gmesh);
    if(simtype != DarcyOrthogonal) {
        auto cmesh_m = CreateMultiphysicsMesh(cmeshH1,cmeshGFem);
        plotfile = "GFem";
        switch (simtype) {
            case Darcy:
            case DarcyNofrac:
                plotfile += "_Darcy";
                break;
            case Elast:
            case ElastNoFrac:
                plotfile += "_Elast";
                break;
            case DarcyOrthogonal:
                DebugStop();
                break;
        }
        Simulate(cmesh_m,plotfile);
        CleanUp(cmesh_m);

    } else {
        auto cmeshH1Mirror = CreateH1MirrorMesh(cmeshGFem);
        auto cmesh_m = CreateMultiphysicsMesh(cmeshH1,cmeshGFem,cmeshH1Mirror);
        plotfile = "GFem";
        switch (simtype) {
            case Darcy:
            case DarcyNofrac:
            case Elast:
            case ElastNoFrac:
                DebugStop();
                break;
            case DarcyOrthogonal:
                plotfile += "_DarcyOrthogonal";
                break;
        }
        ComputeConnectRestraints(cmesh_m);
        Simulate(cmesh_m,plotfile);
        CleanUp(cmesh_m);
    }
    return 0;   
}

/// @brief Read a gmsh file and return a geometric mesh
/// @param filename 
/// @return a geometric mesh
// Physical Surface("dom") = {1};
// //+
// Physical Line("bcB") = {1};
// Physical Line("bcR") = {2};
// Physical Line("bcT") = {3};
// Physical Line("bcL") = {4};

TPZGeoMesh *ReadGmsh(const std::string &meshfilename)
{
    TPZGmshReader gmsh;
    gmsh.SetVerbose(1);
    gmsh.GetDimNamePhysical()[1]["BcB"] = BCb;
    gmsh.GetDimNamePhysical()[1]["BcR"] = BCr;  
    gmsh.GetDimNamePhysical()[1]["BcT"] = BCt;
    gmsh.GetDimNamePhysical()[1]["BcL"] = BCl;
    gmsh.GetDimNamePhysical()[1]["Fracture"] = cutmat;
    gmsh.GetDimNamePhysical()[1]["FractureEdge"] = fracedge;
    gmsh.GetDimNamePhysical()[2]["dom"] = volmat;
    auto gmesh = gmsh.GeometricGmshMesh(meshfilename);
    return gmesh;
}

/// @brief Identify the fracture edge coordinates
void FractureEnds(TPZGeoMesh *gmesh, TPZVec<REAL> &fracend)
{
    int64_t nel = gmesh->NElements();
    for(int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->MaterialId() == fracedge) {
            TPZManVector<REAL,3> center(0,0.);
            gel->X(center,fracend);
            return;
        }
    }
}

/// @brief Refine towards the fracture edge
void RefineTowardsFrac(TPZGeoMesh *gmesh,int nref)
{
    std::set<int> matids = {fracedge};
    for(int iref = 0; iref<nref; iref++)
    {
        TPZRefPatternTools::RefineDirectional(gmesh, matids);
    }
}


void BuildBlueRedElements(TPZGeoMesh *gmesh, std::set<int64_t> &blue, std::set<int64_t> &red) {
    int64_t nel = gmesh->NElements();
    int dim = gmesh->Dimension();
    std::set<TPZGeoEl *> cutelements;
    for(int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->HasSubElement()) continue;
        if(gel->MaterialId() == cutmat) {
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
                    int neighorient = 1;
                    {
                        TPZGeoEl *gel = neighbour.Element();
                        TPZManVector<REAL,3> center(dim,0.);
                        gel->CenterPoint(gel->NSides()-1,center);
                        TPZFNMatrix<9,REAL> gradx(3,dim);
                        gel->GradX(center, gradx);
                        REAL z = gradx(0,0)*gradx(1,1)-gradx(1,0)*gradx(0,1);
                        if(z < 0.) {
                            neighorient = -1;
                        }
                    }
                    bool counterclockwise = neighbour.IsNeighbourCounterClockWise(gelside);
                    if(neighorient == -1) counterclockwise = !counterclockwise;
                    if(counterclockwise) {
                        blue.insert(neighindex);
                    }
                    else {
                        red.insert(neighindex);
                    }

                    // TPZTransform<REAL> tr(1);
                    // gelside.SideTransform3(neighbour, tr);
                    // if(tr.Mult()(0,0) > 0.) {
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
                        if(neighbour.Element()->HasSubElement()) continue;
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
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, color,true);
    }
}

#include "pzvec_extras.h"
/// @brief Choose the elements to be enriched based on their distance to the fracture
/// @param gmesh geometric mesh that contains the elements
/// @param radius distance from the fracture edge
/// @param enrichedelements elements to be enriched
void ElementsToEnrich(TPZGeoMesh *gmesh, TPZVec<REAL> &fracend, REAL radius, std::set<int64_t> &enrichedelements) {
    int64_t nel = gmesh->NElements();
    for(int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->HasSubElement()) continue;
        if(gel->Dimension() != gmesh->Dimension()) continue;
        TPZManVector<REAL,3> center(gel->Dimension(),0.);
        gel->CenterPoint(gel->NSides()-1,center);
        TPZManVector<REAL,3> xcenter(3,0.);
        gel->X(center,xcenter);
        REAL dist = Norm(xcenter-fracend);
        if(dist < radius) {
            enrichedelements.insert(el);
        }
    }
}



/// @brief Create a computational mesh with H1 elements
TPZCompMesh *CreateH1CompMesh(TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    switch(simtype) {
        case Darcy:
        case DarcyNofrac:
            InsertMaterialObjectsDarcy(cmesh);
            break;
        case Elast:
        case ElastNoFrac:
            InsertMaterialObjectsElasticity2D(cmesh);
            break;
        case DarcyOrthogonal:
            InsertMaterialObjectsDarcy(cmesh);
            break;
    }
    cmesh->ExpandSolution();
    cmesh->ComputeNodElCon();
    {
        std::fstream out("cmesh.txt");
        cmesh->Print(out);
    }

    // cmesh->AutoBuild(matidsh1);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    std::set<int> matids = {volmat,BCb,BCt,BCr,BCl};
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

void CreateGFemCompElements(TPZGFemCompMesh *cmesh, std::set<int64_t> &elements, GFemcolors color){
    TPZGeoMesh *gmesh = cmesh->Reference();
    for(auto el : elements) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->HasSubElement()) continue;
        int matid = gel->MaterialId();
        if(!cmesh->FindMaterial(matid)) continue;
        if(gel->Reference()) continue;
        switch(gel->Type())
        {
            case EOned:
            {
                auto *cel = new TPZGFemCompElH1< pzshape::TPZShapeLinear >(*cmesh, gel,color);
            }
            break;
            case ETriangle:
            {
                auto *cel = new TPZGFemCompElH1< pzshape::TPZShapeTriang >(*cmesh, gel, color);
            }
            break;
            case EQuadrilateral:
            {
                auto *cel = new TPZGFemCompElH1< pzshape::TPZShapeQuad >(*cmesh, gel, color);
            }
            break;
            default:
                DebugStop();
        }
    }
}


#include "TPZFrac2d.h"

/// @brief Create a computational mesh with L2 elements
 TPZGFemCompMesh *CreateGFemCompMesh(TPZGeoMesh *gmesh)
 {
    gmesh->ResetReference();
    TPZGFemCompMesh *cmesh = new TPZGFemCompMesh(gmesh);
    int dim = 2;

    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(1);

    cmesh->SetAllCreateFunctionsContinuous();
    switch(simtype) {
        case Darcy:
            InsertMaterialObjectsDarcy(cmesh);
            break;
        case DarcyNofrac:
            break;
        case Elast:
            InsertMaterialObjectsElasticity2D(cmesh);
            break;
        case ElastNoFrac:
            break;
        case DarcyOrthogonal:
            InsertMaterialObjectsDarcy(cmesh);
            break;
    }
    if(simtype == DarcyNofrac || simtype == ElastNoFrac) {
        return cmesh;
    }
     //o método ApproxSpace().CreateDisconnectedElements(true) é chamado para transformar os elementos em elementos discontinuos, criando assim um espaço L2 de aproximação.

    std::set<int64_t> blue, red;
    BuildBlueRedElements(gmesh, blue, red);
    CreateGFemCompElements(cmesh, blue, Black);
    CreateGFemCompElements(cmesh, red, White);
    TPZManVector<REAL> fracend(3,0.);
    FractureEnds(gmesh,fracend);
    std::set<int64_t> enrichedelements;
    ElementsToEnrich(gmesh, fracend, edgelength*3, enrichedelements);
    CreateGFemCompElements(cmesh, enrichedelements, Grey);
    cmesh->DrawElementColors("GFemColor.vtk");
    cmesh->ExpandSolution();
    cmesh->InitializeShapeFunctionMap();
    {
        std::ofstream out("cmeshGFem.txt");
        cmesh->Print(out);
    }
    return cmesh;
 }

/// @brief Create an H1 mirror mesh from the GFem mesh
TPZCompMesh *CreateH1MirrorMesh(TPZGFemCompMesh *cmeshGFem) {
    TPZGeoMesh *gmesh = cmeshGFem->Reference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    switch(simtype) {
        case Darcy:
        case DarcyNofrac:
            InsertMaterialObjectsDarcy(cmesh);
            break;
        case Elast:
        case ElastNoFrac:
            InsertMaterialObjectsElasticity2D(cmesh);
            break;
        case DarcyOrthogonal:
            InsertMaterialObjectsDarcy(cmesh);
            break;
    }
    cmesh->ExpandSolution();
    cmesh->ComputeNodElCon();
    {
        std::fstream out("cmesh.txt");
        cmesh->Print(out);
    }

    // cmesh->AutoBuild(matidsh1);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    std::set<int> matids = {volmat,BCb,BCt,BCr,BCl};
    {
        int64_t nel = cmeshGFem->NElements();
        TPZStack<int64_t > gelstack;
        for(int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmeshGFem->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if(gel->HasSubElement()) DebugStop();
            int matid = gel->MaterialId();
            if(matids.find(matid) != matids.end()) {
                gelstack.Push(gel->Index());
            } else {
                DebugStop();
            }
        }
        cmesh->ApproxSpace().BuildMesh(*cmesh, gelstack);
    }
    cmesh->ExpandSolution();
    cmesh->ComputeNodElCon();
    {
        std::fstream out("cmesh.txt");
        cmesh->Print(out);
    }
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

#include "TPZGFemDarcyFlow.h"

/// @brief Create a computational "multiphysics" mesh with only HDiv elements
 TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(TPZCompMesh* cmeshH1,TPZGFemCompMesh* cmeshGFem, TPZCompMesh *cmeshMirror){
//
    auto gmesh = cmeshH1->Reference();
    TPZMultiphysicsCompMesh *cmesh_m = new TPZMultiphysicsCompMesh(gmesh);

    const int dim = cmeshH1->Dimension();
    const int pord = cmeshH1->GetDefaultOrder();
    cmesh_m->SetDimModel(dim);
    cmesh_m->SetDefaultOrder(pord);
         
    cmesh_m->SetAllCreateFunctionsMultiphysicElem();
    switch(simtype) {
        case Darcy:
        case DarcyNofrac:
        case DarcyOrthogonal:
            InsertMaterialObjectsDarcy2DMF(cmesh_m);
            break;
        case Elast:
        case ElastNoFrac:
            InsertMaterialObjectsElasticity2DMF(cmesh_m);
            break;
    }

    int nmeshes = 2;
    if(cmeshMirror) nmeshes = 3;
    TPZManVector<int, 3> active_approx_spaces(nmeshes, 1);
    TPZManVector<TPZCompMesh *, 3> mesh_vec(nmeshes, 0);
    mesh_vec[0] = cmeshH1;
    mesh_vec[1] = cmeshGFem;
    if(cmeshMirror) {
        mesh_vec[1] = cmeshMirror;
        mesh_vec[2] = cmeshGFem;
    }

    cmesh_m->BuildMultiphysicsSpace(active_approx_spaces, mesh_vec);
    cmesh_m->ExpandSolution();

    {
        std::cout << "Number of equations " << cmesh_m->NEquations() << std::endl;
        int64_t neqsol = cmesh_m->Solution().Rows();
        std::cout << "Solution size " << neqsol << std::endl;
        int64_t neq = cmesh_m->NEquations();
        TPZBlock &block = cmesh_m->Block();
        int64_t ncon = block.NBlocks();
        std::set<int64_t> allseq;
        for(int64_t ic=0; ic<ncon; ic++){
            TPZConnect &c = cmesh_m->ConnectVec()[ic];
            int csz = c.NShape()*c.NState();
            int64_t ib = c.SequenceNumber();
            allseq.insert(ib);
            int sz = block.Size(ib);
            bool cactive = true;
            if(c.HasDependency() || c.IsCondensed()) cactive = false;
            if(csz != sz) {
                std::cout << "Connect " << ib << " size " << sz << " nshape " << c.NShape() << " nstate " << c.NState() << std::endl;
                c.Print(*cmesh_m);
                DebugStop();
            }
            int64_t firsteq = block.Position(ib);
            if(cactive && firsteq+sz > neq) DebugStop();
            if(firsteq+sz > neqsol) DebugStop();
        }
        if(allseq.size() != ncon) DebugStop();
    }
    if(nmeshes == 3) {
        int64_t ncon0 = cmeshH1->NConnects();
        int64_t ncon1 = cmeshMirror->NConnects();
        for (int64_t ic = 0; ic < ncon1; ic++) {
            TPZConnect &c = cmesh_m->ConnectVec()[ic+ncon0];
            c.SetCondensed(true);
        }
    }

    cmesh_m->ExpandSolution();
    cmesh_m->CleanUpUnconnectedNodes();
    {
        std::cout << "Number of equations " << cmesh_m->NEquations() << std::endl;
        int64_t neqsol = cmesh_m->Solution().Rows();
        std::cout << "Solution size " << neqsol << std::endl;
        int64_t neq = cmesh_m->NEquations();
        TPZBlock &block = cmesh_m->Block();
        int64_t ncon = block.NBlocks();
        std::set<int64_t> allseq;
        for(int64_t ic=0; ic<ncon; ic++){
            TPZConnect &c = cmesh_m->ConnectVec()[ic];
            int csz = c.NShape()*c.NState();
            int64_t ib = c.SequenceNumber();
            allseq.insert(ib);
            int sz = block.Size(ib);
            bool cactive = true;
            if(c.HasDependency() || c.IsCondensed()) cactive = false;
            if(csz != sz) {
                std::cout << "Connect " << ib << " size " << sz << " nshape " << c.NShape() << " nstate " << c.NState() << std::endl;
                c.Print(*cmesh_m);
                DebugStop();
            }
            int64_t firsteq = block.Position(ib);
            if(cactive && firsteq+sz > neq) DebugStop();
            if(firsteq+sz > neqsol) DebugStop();
        }
        if(allseq.size() != ncon) DebugStop();
    }

    return cmesh_m;
 }

/// @brief Simulate the NACA profile using H1 approximation
void Simulate(TPZCompMesh *cmesh, std::string name)
{

    TPZLinearAnalysis an(cmesh,RenumType::EMetis);
    //TPZSkylineStructMatrix<STATE> strmat(cmesh);
    TPZSSpStructMatrix<STATE> strmat(cmesh);
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Assemble();

    int nvar = an.Solution().Rows();
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
    if(nstate == 1) {
        fields.Push("Pressure");
        fields.Push("Flux");
    } else if (nstate == 2) {
        fields.Push("SigmaX");
        fields.Push("SigmaY");
        fields.Push("SigmaZ");
        fields.Push("TauXY");
        fields.Push("Displacement");
    } else if (nstate == 3) {
        DebugStop();
    };
    //nesse conjunto de linhas de código, temos que TPZVec é uma estrutura do tipo vetor que contém como argumento uma variável chamda "fields" que é uma lista de strings, que, pelo que se chamam, são relacionadas à pressão e ao fluxo.
    //cria um vetor de strings chamado fields que contém os nomes dos campos que serão pós-processados. Neste caso, os campos incluem "Pressure" (pressão) e "Flux" (fluxo). Esses campos representam propriedades do problema que desejamos visualizar após a simulação.
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
    int dim = 2;
    TPZDarcyFlow *material = new TPZDarcyFlow(volmat,dim);
    cmesh->InsertMaterialObject(material);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(2);

    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);

    auto bnd1 = material->CreateBC(material, BCl, 1, val1, val2);
    TPZBndCondT<STATE> *bndcond = dynamic_cast<TPZBndCondT<STATE> *>(bnd1);
    cmesh->InsertMaterialObject(bnd1);
    auto bnd2 = material->CreateBC(material, BCr, 1, val1, val2);
    cmesh->InsertMaterialObject(bnd2);
    val2[0] = -1.;
    auto bnd3 = material->CreateBC(material, BCb, 0, val1, val2);
    cmesh->InsertMaterialObject(bnd3);
    val2[0] = 1.;
    auto bnd4 = material->CreateBC(material, BCt, 0, val1, val2);
    cmesh->InsertMaterialObject(bnd4);

    TPZGFemCompMesh *gfem = dynamic_cast<TPZGFemCompMesh *>(cmesh);
    if(gfem) {
        TPZManVector<REAL,3> fracend(3,0.),fracdir(3,0.);
        FractureEnds(cmesh->Reference(),fracend);
        fracdir[0] = 1.;
        REAL nu = 0.3;
        int planestress = 1;
        gfem->fFrac.SetFracData(fracend,fracdir,nu,Scalar2d, planestress);
    }

}

/// @brief Insert material objects int the multiphysics mesh
void InsertMaterialObjectsDarcy2DMF(TPZMultiphysicsCompMesh *cmesh_m){
    int dim = cmesh_m->Dimension();
    std::set<int> materialIDs;
    TPZGFemDarcyFlow *material = new TPZGFemDarcyFlow(volmat,dim);
    //  material->SetBigNumber(10e8);
     cmesh_m->InsertMaterialObject(material);
     materialIDs.insert(volmat);

     // 2. Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);

    // 2.1 Neumann condition on the profile
    auto bnd1 = material->CreateBC(material, BCl, 1, val1, val2);
    cmesh_m->InsertMaterialObject(bnd1);

    // // 2.2 Neumann condition on the infity boundary
    auto bnd2 = material->CreateBC(material, BCr, 1, val1, val2);
    cmesh_m->InsertMaterialObject(bnd2);
    val2[0] = 1.;
    auto bnd3 = material->CreateBC(material, BCt, 0, val1, val2);
    cmesh_m->InsertMaterialObject(bnd3);
    val2[0] = -1.;
    auto bnd4 = material->CreateBC(material, BCb, 0, val1, val2);
    cmesh_m->InsertMaterialObject(bnd4);

}

#include "Elasticity/TPZElasticity2D.h"

/// @brief Insert material objects int the multiphysics mesh
void InsertMaterialObjectsElasticity2D(TPZCompMesh *cmesh_m){
    int dim = cmesh_m->Dimension();
    std::set<int> materialIDs;
    STATE E = 100., nu = 0.3;
    STATE fx = 0., fy = 0.;
    int planestress = 1;
    TPZElasticity2D *material = new TPZElasticity2D(volmat,E,nu,fx,fy,planestress);
    //  material->SetBigNumber(10e8);
     cmesh_m->InsertMaterialObject(material);
     materialIDs.insert(volmat);

     // 2. Boundary Conditions
    TPZFMatrix<STATE> val1(dim,dim,0.);
    TPZManVector<STATE> val2(dim,0.);

    // 2.1 mixed on the left
    val1(0,0) = 1.;
    auto bnd1 = material->CreateBC(material, BCl, 2, val1, val2);
    cmesh_m->InsertMaterialObject(bnd1);
    val1.Zero();

    // 2.2 Neumann condition on the right
    auto bnd2 = material->CreateBC(material, BCr, 1, val1, val2);
    cmesh_m->InsertMaterialObject(bnd2);
    // 2.3 Neumann traction condition on the top
    val2[1] = 1.;
    auto bnd3 = material->CreateBC(material, BCt, 1, val1, val2);
    cmesh_m->InsertMaterialObject(bnd3);
    val2.Fill(0.);
    // 2.4 Mixed condition on the bottom
    val1(1,1) = 1000.;
    auto bnd4 = material->CreateBC(material, BCb, 2, val1, val2);
    cmesh_m->InsertMaterialObject(bnd4);

    TPZGFemCompMesh *gfem = dynamic_cast<TPZGFemCompMesh *>(cmesh_m);
    if(gfem) {
        TPZManVector<REAL,3> fracend(3,0.),fracdir(3,0.);
        FractureEnds(cmesh_m->Reference(),fracend);
        fracdir[0] = 1.;
        gfem->fFrac.SetFracData(fracend,fracdir,nu,Vector2d, planestress);
    }
}

#include "TPZGFemElasticity2D.h"

/// @brief Insert material objects int the multiphysics mesh
void InsertMaterialObjectsElasticity2DMF(TPZMultiphysicsCompMesh *cmesh_m){
    int dim = cmesh_m->Dimension();
    std::set<int> materialIDs;
    STATE E = 100., nu = 0.3;
    STATE fx = 0., fy = 0.;
    TPZGFemElasticity2D *material = new TPZGFemElasticity2D(volmat,E,nu,fx,fy);
    //  material->SetBigNumber(10e8);
     cmesh_m->InsertMaterialObject(material);
     materialIDs.insert(volmat);

     // 2. Boundary Conditions
    TPZFMatrix<STATE> val1(dim,dim,0.);
    TPZManVector<STATE> val2(dim,0.);

    // 2.1 mixed on the left
    val1(0,0) = 1.;
    auto bnd1 = material->CreateBC(material, BCl, 2, val1, val2);
    cmesh_m->InsertMaterialObject(bnd1);
    val1.Zero();

    // 2.2 Neumann condition on the right
    auto bnd2 = material->CreateBC(material, BCr, 1, val1, val2);
    cmesh_m->InsertMaterialObject(bnd2);
    // 2.3 Neumann traction condition on the top
    val2[1] = 1.;
    auto bnd3 = material->CreateBC(material, BCt, 1, val1, val2);
    cmesh_m->InsertMaterialObject(bnd3);
    val2.Fill(0.);
    // 2.4 Mixed condition on the bottom
    val1(1,1) = 1000.;
    auto bnd4 = material->CreateBC(material, BCb, 2, val1, val2);
    cmesh_m->InsertMaterialObject(bnd4);
}

#include "TPZGFemOrthogonal.h"

void ComputeConnectRestraints(TPZMultiphysicsCompMesh *cmesh_m)
{
    cmesh_m->CleanUpUnconnectedNodes();
    {
        std::cout << "Number of equations " << cmesh_m->NEquations() << std::endl;
        int64_t neqsol = cmesh_m->Solution().Rows();
        std::cout << "Solution size " << neqsol << std::endl;
        int64_t neq = cmesh_m->NEquations();
        TPZBlock &block = cmesh_m->Block();
        int64_t ncon = block.NBlocks();
        std::set<int64_t> allseq;
        for(int64_t ic=0; ic<ncon; ic++){
            TPZConnect &c = cmesh_m->ConnectVec()[ic];
            int csz = c.NShape()*c.NState();
            int64_t ib = c.SequenceNumber();
            allseq.insert(ib);
            int sz = block.Size(ib);
            bool cactive = true;
            if(c.HasDependency() || c.IsCondensed()) cactive = false;
            if(csz != sz) {
                std::cout << "Connect " << ib << " size " << sz << " nshape " << c.NShape() << " nstate " << c.NState() << std::endl;
                c.Print(*cmesh_m);
                DebugStop();
            }
            int64_t firsteq = block.Position(ib);
            if(cactive && firsteq+sz > neq) DebugStop();
            if(firsteq+sz > neqsol) DebugStop();
        }
        if(allseq.size() != ncon) DebugStop();
    }
    TPZLinearAnalysis an(cmesh_m,RenumType::ENone);
    //TPZSkylineStructMatrix<STATE> strmat(cmesh_m);
    TPZSSpStructMatrix<STATE> strmat(cmesh_m);
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Assemble();
    TPZMatrixSolver<STATE> &matsolver = an.MatrixSolver<STATE>();
    TPZMatrix<STATE> *global_mat = matsolver.Matrix().operator->();
    TPZGFemOrthogonal ortho(cmesh_m,global_mat);
    ortho.OrthogonalizeConnects();
    global_mat->Zero();
    an.Assemble();
    ortho.VerifyOrthogonality();
    int64_t gfemconnectindex = ortho.LargestEigenvalueRatio();
    ortho.DrawOrthogonalization(gfemconnectindex,"Ortho");
}