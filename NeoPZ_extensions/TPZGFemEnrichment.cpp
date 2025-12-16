#include "TPZGFemEnrichment.h"
#include "pzgeoelbc.h"
#include "pzgeoel.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "pzgmesh.h"

#include "pzinterpolationspace.h"
#include "pzintel.h"


TPZGFemEnrichment::~TPZGFemEnrichment() {
}

#include <list>
#include <tuple>
#include "TPZGeoElement.h"
#include "TPZRefLinear.h"
#include "pzgeotriangle.h"
#include "tpzgeoblend.h"
#include "tpzarc3d.h"

    /// Load configuration from JSON file
void TPZGFemEnrichment::LoadFromJSON(json &input) {
    // read radius of influence
    fRadius = input["radius"];
    // read angles
    fAngles.clear();
    for(auto &item : input["angles"]) {
        REAL angle = item;
        fAngles.push_back(angle);
    }
    for(int i=0; i<fAngles.size()-1; i++) {
        if(fAngles[i+1] <= fAngles[i]) DebugStop();
    }
    if(fabs(fAngles[fAngles.size()-1]-fAngles[0]-2.*M_PI) > 1.e-6) {
        std::cout << "Last angle " << fAngles[fAngles.size()-1] << " first angle " << fAngles[0] << std::endl;
        std::cout << "Difference error " << fAngles[fAngles.size()-1]-fAngles[0]-M_2_PI << std::endl;
        DebugStop();
    }
    fAngles[fAngles.size()-1] = fAngles[0] + 2.*M_PI;
    fMaterialIds.clear();
    if(input.find("sector material ids") == input.end()) {
        std::cout << "Error! sector material ids not found in input JSON." << std::endl;
        DebugStop();
    }
    for(auto &item : input["sector material ids"]) {
        int matid = item;
        fMaterialIds.push_back(matid);
    }
    if(fMaterialIds.size() != fAngles.size()-1) {
        std::cout << "Number of material ids " << fMaterialIds.size() << " number of sectors " << fAngles.size()-1 << std::endl;
        DebugStop();
    }
    // ensure angles are sorted
    // set max upper and min lower sectors
    fMaxUpper = fAngles.size()-2;
    fMinLower = 1;
    int sector = 0;
    CreateGeometry(input);
    // std::cout << "nel " << fGeoMesh->NElements() << std::endl;
    TPZBuildSBFem builder(fGeoMesh);
    std::cout << "before configure sbfem builder" << std::endl;
    ConfigureSBFemBuilder(builder, input);
    fCompMesh = CreateSBFemSpace(builder, input);
    ComputeSBFemApproximationSpace();
    IdentifySingularModes();
}


TPZGeoMesh *TPZGFemEnrichment::CreateGeometry(json &input) {

    std::list<std::tuple<REAL,REAL,int>> sectors;
    for(int i=0; i< fMaterialIds.size(); i++) {
        sectors.push_back({fAngles[i],fAngles[i+1],fMaterialIds[i]});
    }
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    // we assume the centernode is at 0,0
    int64_t centerindex = gmesh->NodeVec().AllocateNewElement();
    TPZManVector x(3,0.);
    gmesh->NodeVec()[centerindex].Initialize(x,*gmesh);
    int pointmatid = fPointMatId;
    if(input.find("point_singularity") != input.end()) {
        pointmatid = input["point_singularity"];
        fPointMatId = pointmatid;
    }
    {
        TPZManVector<int64_t> nodeindex = {centerindex};
        int64_t index;
        gmesh->CreateGeoElement(EPoint,nodeindex,pointmatid,index);
    }


    // generate the arc nodes and arc elements
    int arcmatid = fArcMatId;
    int64_t ncirclenodes = sectors.size()+1;
    x[0] = fRadius;
    TPZManVector<int64_t> nodeindexes(ncirclenodes,-1);
    nodeindexes[0] = gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[nodeindexes[0]].Initialize(x,*gmesh);
    int count = 1;
    for(auto &it : sectors) {
        nodeindexes[count] = gmesh->NodeVec().AllocateNewElement();
        REAL prevangle = std::get<0>(it);
        REAL angle = std::get<1>(it);
        x[0] = fRadius*std::cos(angle);
        x[1] = fRadius*std::sin(angle);
        gmesh->NodeVec()[nodeindexes[count]].Initialize(x,*gmesh);
        int matid = std::get<2>(it);
        TPZManVector<int64_t> lineindexes(3, -1);
        lineindexes[0] = nodeindexes[count-1];
        lineindexes[1] = nodeindexes[count];
        lineindexes[2] = lineindexes[0];
        {
            // create the arc elements
            auto *gel = new TPZGeoElement<pzgeom::TPZArc3D,pzrefine::TPZRefLinear> (lineindexes,arcmatid,*gmesh);
            TPZManVector<REAL> ax1(3,0.), ax2(3,0.);
            gmesh->NodeVec()[lineindexes[0]].GetCoordinates(ax1);
            ax1[0] /= fRadius;
            ax1[1] /= fRadius;
            ax2[0] = -ax1[1];
            ax2[1] = ax1[0];
            TPZManVector<REAL> center(3,0.);
            TPZFMatrix<REAL> rotationmatrix(3,3,0.);
            rotationmatrix.Identity();

            for(int i=0; i<3; i++) {
                rotationmatrix(i,0) = ax1[i];
                rotationmatrix(i,1) = ax2[i];
            }
            REAL openingangle = angle - prevangle;
            gel->Geom().SetAttributes(gmesh, center, rotationmatrix, fRadius, openingangle);
        }
        // create the triangular element
        TPZManVector<int64_t> triindexes(lineindexes);
        triindexes[2] = centerindex;
        int64_t index;
        gmesh->CreateGeoBlendElement(ETriangle,triindexes,matid+fMaterialIdOffset,index);
        prevangle = angle;
        count++;
    }

    gmesh->BuildConnectivity();
    
    fGeoMesh = gmesh;
    return gmesh;
}

/// @brief configure the SBFem builder according to the json input
/// @param builder reference to the SBFem builder
/// @param input json object with the input data
void TPZGFemEnrichment::ConfigureSBFemBuilder (TPZBuildSBFem &builder, json &input) {
    int skeletonmatid = fSkeletonMatId;
    if(input.find("skeletonMatid") != input.end()) {
        skeletonmatid = input["skeletonMatid"];
        fSkeletonMatId = skeletonmatid;
    }
    
    builder.SetSkeletonMatid(skeletonmatid);


    std::map<int,int> matidtranslation;
    std::set<int> sbfemids;
    for (const auto &item : fMaterialIds) {
        matidtranslation[fMaterialIdOffset + item] = item;
        sbfemids.insert(item);
    }
    builder.SetMatIdTranslation(matidtranslation);
    std::set<int> boundmatids;
    boundmatids.insert(fArcMatId);
    builder.SetBoundaryMatIds(boundmatids);
    int skeletonporder = 6;
    if(input.find("skeleton porder") != input.end()) {
        skeletonporder = input["skeleton porder"];
    }
    builder.SetSkeletonPOrder(skeletonporder);
    int nref = 1;
    if(input.find("nref_skeleton") != input.end()) {
        nref = input["nref_skeleton"];
    }
    builder.SetNRefSkeleton(nref);
}

#include "TPZBndCondT.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "Elasticity/TPZElasticity2D.h"
#include "TPZNullMaterial.h"

void TPZGFemEnrichment::InsertMaterialObjectsH1(TPZCompMesh &cmesh, json &input)
{
    int dim = cmesh.Dimension();
    TPZMaterialT<STATE> *mat = nullptr;
    int nstate = 0;
    std::string problemtype = "elastic";
    if(input.find("problem_type") != input.end()) {
        problemtype = input["problem_type"];
        std::cout << "Problem type " << problemtype << std::endl;
    }

    for (const auto &item : input["materials"]) {
        int matid = item["matid"];
        if(problemtype == "Scalar") {
            nstate = 1;
            REAL perm = item["permeability"];
            auto darcy = new TPZDarcyFlow(matid, dim);
            //        auto darcy = new TPZDarcyFlow(matid, dim);
            darcy->SetConstantPermeability(perm);
            cmesh.InsertMaterialObject(darcy);
            if(!mat) mat = darcy;
        } else if (problemtype == "Elastic") {
            nstate = 2;
            bool planeStress = true;
            REAL fx(0.), fy(0.);
            {
                if(item["tension state"] == "plane stress") {
                    planeStress = true;
                } else if (item["tension state"] == "plane strain") {
                    planeStress = false;
                } else {
                    DebugStop();
                }
            }
            if(item.find("isotropic") != item.end()) {
                REAL young = item["isotropic"]["young"];
                REAL poisson = item["isotropic"]["poisson"];
                TPZElasticity2D *matelast = new TPZElasticity2D(matid,young,poisson,fx,fy,planeStress);
                if(!mat) mat = matelast;
                cmesh.InsertMaterialObject(matelast);
            } else if(item.find("orthotropic") != item.end()) {
                REAL Ex = item["orthotropic"]["Ex"];
                REAL Ey = item["orthotropic"]["Ey"];
                REAL nuxy = item["orthotropic"]["nuxy"];
                REAL Gxy = item["orthotropic"]["Gxy"];
                TPZLinearElasticityConstitutive orthotropic;
                orthotropic.SetOrthotropicProperties(Ex,Ey,nuxy,Gxy);
                orthotropic.SetPlaneStress(planeStress);
                TPZElasticity2D *matelast = new TPZElasticity2D(matid);
                matelast->SetConstitutiveLaw(orthotropic);
                if(!mat) mat = matelast;
                cmesh.InsertMaterialObject(matelast);
            } else {
                DebugStop();
            }
            

            nstate = 2;
            
        } else {
            std::cout << "Problem type " << problemtype << " is unknown\n";
            DebugStop();
        }
    }
    for (const auto &item : input["boundary_conditions"]) {
        int matid = item["matid"];
        std::string typeName = item["name"];
        std::vector<double> valvec = item["value"];
        
        if(valvec.size() != nstate) DebugStop();

        int type = item["type"];
        TPZFNMatrix<4,STATE> val1(nstate,nstate,0.);
        TPZManVector<STATE,3> val2(nstate,0.);
        for(int i = 0; i<nstate; i++) val2[i] = valvec[i];
        if (type == 0) {
        } else if (type == 1) {
        } else if (type == 2) {
            DebugStop();
        }
        TPZBndCondT<STATE> *bc = mat->CreateBC(mat, matid, type, val1, val2);
        cmesh.InsertMaterialObject(bc);
    }

    TPZNullMaterial<STATE> *arcmat = new TPZNullMaterial<STATE>(fArcMatId, dim,nstate);
    cmesh.InsertMaterialObject(arcmat);
    int skeletonMatId = fSkeletonMatId;
    TPZNullMaterial<STATE> *skelmat = new TPZNullMaterial<STATE>(skeletonMatId, dim,nstate);
    cmesh.InsertMaterialObject(skelmat);
}

#include "TPZSBFemElementGroup.h"

/// @brief Create a computational mesh based on the SBFem hybrid H1 space
/// @param gmesh autopointer to the geometric mesh
/// @param input json object with the input data
/// @return pointer to the computational mesh
TPZCompMesh *TPZGFemEnrichment::CreateSBFemSpace(TPZBuildSBFem &builder, json &input) {
    TPZAutoPointer<TPZGeoMesh> gmesh = fGeoMesh;
    if(!gmesh) {
        DebugStop();
    }
    int dim = gmesh->Dimension();
    ConfigureSBFemBuilder(builder, input);
    int64_t singular_partition = fPointMatId;
    if(input.find("point_singularity") != input.end()) {
        singular_partition = input["point_singularity"];
        fPointMatId = singular_partition;
    }
    // two steps :
    // 1 create standard elements for all elements that do not touch the singularity
    int singular_matid = fPointMatId;
    TPZStack<int64_t> volels;
    std::map<int,int> matidtranslation = builder.GetMatIdTranslation();
    
    int64_t singular_node = -1;
    TPZGeoEl *singular_element = 0;
    int64_t nel = gmesh->NElements();
    std::set<int64_t> singular_group;
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        // std::cout << " el " << el << " matid " << gel->MaterialId() << std::endl;
        if(!gel || gel->MaterialId() != singular_matid) continue;
        singular_node = gel->NodeIndex(0);
        singular_element = gel;
        break;
    }
    if(singular_node == -1) DebugStop();
    TPZGeoElSide singularside (singular_element);
    std::set<int64_t> to_delete;
    for (auto gelside = singularside.Neighbour(); gelside != singularside; gelside++) {
        TPZGeoEl *gel = gelside.Element();
        if(!gel || gel->HasSubElement()) continue;
        if(gel->Dimension() == dim) {
            singular_group.insert(gel->Index());
        } else if(gel->Dimension() == dim-1) {
            // boundary element touching the singularity
            to_delete.insert(gel->Index());
        }
    }
    for(auto it : to_delete) {
        TPZGeoEl *gel = gmesh->Element(it);
        gel->RemoveConnectivities();
        delete gel;
    }

    nel = gmesh->NElements();
    for (auto it : singular_group) {
        volels.Push(it);
    }
    singular_partition = builder.AddPartition(singular_group, singular_node);
    builder.AddSkeletonElements();
    // 2 create collapsed elements towards the singular point
    // builder.SetNRefSkeleton(0);
    int nref = builder.GetNRefSkeleton();
    builder.DivideSkeleton(nref);
    int nsectors = fAngles.size()-1;
    int ndivisions = (1 << nref);
    int nsectorscreated = nsectors * ndivisions;
    TPZManVector<REAL> anglesextended(nsectorscreated+1,-1);
    for(int isector = 0; isector < nsectors; isector++) {
        anglesextended[isector*ndivisions] = fAngles[isector];
        anglesextended[(isector+1)*ndivisions] = fAngles[isector+1];
        REAL delangle = (fAngles[isector+1]-fAngles[isector])/REAL(ndivisions);
        for(int idiv = 1; idiv < ndivisions; idiv++) {
            anglesextended[isector*ndivisions + idiv] = fAngles[isector] + delangle*REAL(idiv);
        }
    }
    REAL middleangle = 0.5*(anglesextended[0]+anglesextended[nsectorscreated]);
    for(int i=0; i< nsectorscreated; i++) {
        REAL angle0 = anglesextended[i];
        REAL angle1 = anglesextended[i+1];
        REAL center = 0.5*(angle0+angle1);
        if(center < middleangle) {
            fMinLower = i+1;
        } else {
            break;
        }
    }
    for(int i=nsectorscreated-1; i >=0; i--) {
        REAL angle0 = anglesextended[i];
        REAL angle1 = anglesextended[i+1];
        REAL center = 0.5*(angle0+angle1);
        if(center > middleangle){
            fMaxUpper = i-1;
        } else {
            break;
        }
    }
    fAngles = anglesextended;
    // std::cout << "Number of sectors after refinement " << nsectorscreated << " angles : ";
    // for(auto angle : fAngles) {
    //     std::cout << angle << " ";
    // }
    // std::cout << std::endl;
    // auto matidtranslation = builder.GetMatIdTranslation();
    // std::set<int> matids;
    // for(auto it : matidtranslation) {
    //     matids.insert(it.first);
    // }
    // builder.GenerateCollapsedGeometricElements(matids);
    TPZCompMesh *cmeshptr = new TPZCompMesh(gmesh);
    TPZCompMesh &cmesh = *cmeshptr;
    InsertMaterialObjectsH1(cmesh, input);
    builder.BuildComputationMesh(cmesh);
    
    {
        int64_t nel = cmesh.NElements();
        for(int64_t el = 0; el < nel ; el++) {
            TPZCompEl *cel = cmesh.Element(el);
            if(!cel) continue;
            TPZSBFemElementGroup *sbfemgroup = dynamic_cast<TPZSBFemElementGroup *>(cel);
            if(!sbfemgroup) continue;
            fElementGroup = sbfemgroup;
        }
    }
    auto &elvec = fElementGroup->GetElGroup();
    if(elvec.size() != fAngles.size()-1) {
        std::cout << "Number of SBFem elements " << elvec.size() << " number of sectors " << fAngles.size()-1 << std::endl;
        DebugStop();
    }
    for(auto &el : elvec) {
        // std::cout << " Element index " << el->Index() << " matid " << el->Material()->Id() << std::endl;
        TPZSBFemVolume *sbfemvol = dynamic_cast<TPZSBFemVolume *>(el);
        if(!sbfemvol) {
            DebugStop();
        }
        TPZGeoEl *gel = sbfemvol->Reference();
        int face = gel->FirstSide(1);
        TPZGeoElSide gelside(gel,face);
        TPZGeoElSide neighbour = gelside.Neighbour();
        // std::cout << " Neighbour element index " << neighbour.Element()->Index() << " matid " << neighbour.Element()->MaterialId() << std::endl;

        TPZCompEl *neighcel = neighbour.Element()->Reference();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(neighcel);
        if(!intel) {
            DebugStop();
        }
        {
            int sector = fSkeletonElements.size();
            TPZGeoEl *gel = intel->Reference();
            TPZManVector<REAL,3> x(3,0.);
            gel->NodePtr(0)->GetCoordinates(x);
            REAL angle = fAngles[sector];
            TPZManVector<REAL,3> expectedx(3,0.);
            expectedx[0] = fRadius*std::cos(angle);;
            expectedx[1] = fRadius*std::sin(angle);
            REAL dist = 0.;
            for(int i=0; i<3; i++) {
                dist += (x[i]-expectedx[i])*(x[i]-expectedx[i]);
            }
            if(dist > 1.e-6) {
                std::cout << "Sector " << sector << " angle " << angle << " expectedx " << expectedx[0] << " " << expectedx[1] << " got ";
                for(int i=0; i<3; i++) {
                    std::cout << x[i] << " ";
                }
                std::cout << " dist error " << dist << std::endl;
                DebugStop();
            }
            gel->NodePtr(1)->GetCoordinates(x);
            angle = fAngles[sector+1];
            expectedx[0] = fRadius*std::cos(angle);
            expectedx[1] = fRadius*std::sin(angle);
            dist = 0.;
            for(int i=0; i<3; i++) {
                dist += (x[i]-expectedx[i])*(x[i]-expectedx[i]);
            }
            if(dist > 1.e-6) {
                std::cout << "Sector " << sector << " angle " << angle << " expectedx " << expectedx[0] << " " << expectedx[1] << " got ";
                for(int i=0; i<3; i++) {
                    std::cout << x[i] << " ";
                }
                std::cout << " dist error " << dist << std::endl;
                DebugStop();
            }
        }
        fSkeletonElements.push_back(intel);
    }
    // std::set<int> matidsskel = builder.GetBoundaryMatIds();
    // matidsskel.insert(builder.GetSkeletonMatid());
    // cmesh.SetDefaultOrder(builder.GetSkeletonPOrder());
    // cmesh.SetAllCreateFunctionsContinuous();
    // cmesh.AutoBuild(matidsskel);
    // builder.CreateVolumetricElements(cmesh);
    // builder.CreateElementGroups(cmesh);
    fCompMesh = cmeshptr;
    return cmeshptr;
}


#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZLinearAnalysis.h"
#include "TPZElementMatrixT.h"

/// @brief Compute the stiffness matrix and load the singular mode as the mesh solution
void TPZGFemEnrichment::ComputeSBFemApproximationSpace() {
    if(!fElementGroup) {
        DebugStop();
    }
    TPZElementMatrixT<STATE> ek,ef;
    fElementGroup->CalcStiff(ek,ef);
}

/// @brief load the singular mode as the mesh solution
void TPZGFemEnrichment::IdentifySingularModes() {

    REAL range[] = {0.01,0.9};
    /// find the SBFEmElementGroup

    if(!fElementGroup) {
        DebugStop();
    }
    // find an eigenvalue larger than zero and smaller than 1
    TPZManVector<std::complex<REAL>,10> eigenvalues;
    eigenvalues = fElementGroup->EigenValues();
    int64_t neigen = eigenvalues.size();
    int64_t selected = -1;
    for (int64_t ieigen = 0; ieigen<neigen; ieigen++) {
        REAL realpart = eigenvalues[ieigen].real();
        if(abs(realpart) > range[0] && abs(realpart) < range[1]) {
            REAL imagpart = eigenvalues[ieigen].imag();
            if (abs(imagpart) < 1e-6) {
                fEigenvectorIndexes.push_back(ieigen);
                fPartOfEigenvector.push_back(EPartOfEigenvector::EOnlyReal);
            } else {
#ifdef PZDEBUG
                if(ieigen == neigen-1 || fabs(eigenvalues[ieigen].imag()-eigenvalues[ieigen+1].imag()) > 1.e-6) {
                    DebugStop();
                }
#endif
                fEigenvectorIndexes.push_back(ieigen);
                fPartOfEigenvector.push_back(EPartOfEigenvector::ERealPart);
                fPartOfEigenvector.push_back(EPartOfEigenvector::EImagPart);
                ieigen++;
            }
        }
    }
    std::cout << "Found " << fEigenvectorIndexes.size() << " singular modes in the range (" << range[0] << "," << range[1] << ")\n";
    std::cout << "Selected Eigenvalues:\n";
    for (int i=0; i<fEigenvectorIndexes.size(); i++) {
        auto modeindex = fEigenvectorIndexes[i];
        std::cout << "  " << modeindex << " : " << eigenvalues[modeindex] << " part of eigenvector " << fPartOfEigenvector[i] << "\n";
    }
}

#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"

/// @brief Plot the solution of the computational mesh
/// @param cmesh reference to the computational mesh
/// @param output json object with the output data
void TPZGFemEnrichment::PlotSolution(const std::string &rootname, TPZCompMesh &cmesh, int refstep, json &output) {
    std::set<int> matiderror;
    if (output.find("matid_translation") != output.end()) {
        for (const auto &item : output["matid_translation"]) {
            int matid = item["to"];
            matiderror.insert(matid);
        }
    } else {
        DebugStop();
    }    // make a directory based on the exact solution
    
    std::string dirname = "PostProcess";
    if(output.find("directory") != output.end()) {
        std::string exactname = output["directory"];
        dirname = exactname;
    }
    int matid = *matiderror.rbegin();
    TPZMaterial *mat = cmesh.FindMaterial(matid);
    if (!mat) {
        std::cout << "Material id " << matid << " not found in the computational mesh\n";
        DebugStop();
    }
    TPZDarcyFlow *darcy = dynamic_cast<TPZDarcyFlow *>(mat);
    TPZElasticity2D *elast = dynamic_cast<TPZElasticity2D *>(mat);
    if (!darcy && !elast) {
        std::cout << "Material id " << matid << " is neither a Darcy flow or elasticity material\n";
        DebugStop();
    }
    std::system((std::string("mkdir -p ") + dirname).c_str());
    // std::string plotfile = dirname + std::string("/Solution.vtk");
    std::string plotfile = dirname + "/" + rootname;

    int skeletonporder = output["skeleton porder"];
    int nref_skeleton = output["nref_skeleton"];

    plotfile +=  "_rS" + std::to_string(nref_skeleton) + "_pS" + std::to_string(skeletonporder);

    TPZStack<std::string> postprocess;
    if (output.find("PostProcess") != output.end()) {
        for (const auto &item : output["PostProcess"]) {
            std::string field = item;
            
            if(darcy && darcy->VariableIndex(field) != -1) {
                postprocess.Push(field);
            }
            if(elast && elast->VariableIndex(field) != -1) {
                postprocess.Push(field);
            }
        }
    }
    int dim = 2;
    int res = 5;
    std::cout << "Post processing in file " << plotfile << std::endl;
    TPZVTKGenerator vtk(&cmesh, postprocess, plotfile, res, dim);
    vtk.SetStep(refstep);
    vtk.Do();
    {
        plotfile += ".geom.vtk";
        std::ofstream out(plotfile);
        TPZVTKGeoMesh::PrintGMeshVTK(cmesh.Reference(), out);
    }
}

/// @brief Load an eigenmode from a file into the computational mesh
void TPZGFemEnrichment::LoadEigenMode(int64_t modeindex) {

    if(!fElementGroup) {
        DebugStop();
    }
    fElementGroup->LoadEigenVector(modeindex);
}

/// The enrichment function that will be returned by GetEnrichmentFunction
void TPZGFemEnrichment::EnrichmentFunction(TPZGFemEnrichment::MSectorFamily family, 
        const std::vector<REAL> &x, 
        std::vector<std::array<REAL, 2>> &phi, 
                std::vector<std::array<REAL ,4>> &dphixy) {
        REAL angle = std::atan2(x[1],x[0]);
        REAL minangle = fAngles[0];
        REAL maxangle = fAngles[fAngles.size()-1];
        while(angle < minangle) {
            angle += 2.*M_PI;
        }
        while(angle > maxangle) {
            angle -= 2.*M_PI;
        }
        REAL refangle = angle;
        if(angle > minangle+M_PI && family == EUpper) {
            angle = minangle;
        }

        if(angle < maxangle -M_PI && family == ELower) {
            angle = maxangle;
        }

        int sector = Sector(angle, family);
        REAL r = R(x);
        REAL ksi = Ksi(r);
        REAL eta = Eta(sector, angle);
        TPZSBFemVolume *vol = dynamic_cast<TPZSBFemVolume *>(fElementGroup->GetElGroup()[sector]);
        TPZManVector<REAL> qsi = {eta, 1.-2.*ksi};
#ifdef PZDEBUG
        {
            if(!vol) {
                DebugStop();
            }
            TPZGeoEl *gel = vol->Reference();
            TPZManVector<REAL,3> locx(3,0.),etav = {eta,(-2.*ksi+1)};
            gel->X(etav,locx);
            REAL diff = (locx[0]-x[0])*(locx[0]-x[0]) +
                        (locx[1]-x[1])*(locx[1]-x[1]);
            if(diff > 1.e-6) {
                int nnodes = gel->NNodes();
                for(int inode = 0; inode < nnodes; inode++) {
                    TPZManVector<REAL,3> nodeco(3,0.);
                    gel->NodePtr(inode)->GetCoordinates(nodeco);
                    std::cout << " Node " << inode << " coordinates " << nodeco[0] << " " << nodeco[1] << " " << nodeco[2] << std::endl;
                }
                std::cout << "Computed point " << locx[0] << " " << locx[1] << " expected point " << x[0] << " " << x[1] << " diff " << diff << std::endl;
                DebugStop();
            }
        }
#endif
        TPZMaterialDataT<REAL> data;
        int neigen = fEigenvectorIndexes.size();
        for(int ieigen = 0; ieigen<neigen; ieigen++) {
            int eigenindex = fEigenvectorIndexes[ieigen];
            LoadEigenMode(eigenindex);
            vol->ComputeSolution(qsi, data, false);
            phi.push_back({data.sol[0][0], data.sol[0][1]});
            TPZFNMatrix<4,REAL> dsolxy(2,2,0.);
            TPZAxesTools<REAL>::Axes2XYZ(data.dsol[0], dsolxy, data.axes);
            dphixy.push_back({dsolxy(0,0), dsolxy(1,0), dsolxy(0,1), dsolxy(1,1)});
        }
    }

    /// @brief Plot the enrichment function
    void TPZGFemEnrichment::PlotEnrichmentFunction(const std::string &filename, json &input) {
        GFemEnrichmentFunction enrichFunc = GetEnrichmentFunction();
        std::string plotfile = filename;
        TPZManVector<std::string,10> postprocnames;
        if(input.find("PostProcess") != input.end()) {
            for (const auto &item : input["PostProcess"]) {
                std::string field = item;
                postprocnames.push_back(field);
            }
        }
        std::string directoryname = input["directory"];
        std::system((std::string("mkdir -p ") + directoryname).c_str());
        plotfile = directoryname + "/" + plotfile;
        int resolution =3;
        if(input.find("PostProcessResolution") != input.end()) {
            resolution = input["PostProcessResolution"];
        }
        int dim = 2;
        int neigen = fEigenvectorIndexes.size();
        for(int ieigen = 0; ieigen<neigen; ieigen++) {
            LoadEigenMode(fEigenvectorIndexes[ieigen]);
            TPZVTKGenerator vtk(fCompMesh, postprocnames, plotfile, resolution, dim);
            vtk.SetStep(ieigen);
            vtk.Do();
        }
    }


        /// Compute enrichment function and its derivatives
void TPZGFemEnrichment::ComputeEnrichmentWithDerivatives(TPZGFemEnrichment::MSectorFamily family, const std::vector<REAL> &x, 
                                          std::vector<std::array<REAL, 2>> &phi, 
                                          std::vector<std::array<REAL ,4>> &dphi) {

    REAL angle = std::atan2(x[1],x[0]);
    REAL minangle = fAngles[0];
    REAL maxangle = fAngles[fAngles.size()-1];
    while(angle < minangle) {
        angle += 2.*M_PI;
    }
    while(angle > maxangle) {
        angle -= 2.*M_PI;
    }
    REAL refangle = angle;
        // if(sector == ELower && pos <= fMaxUpper) return nangles-1;
        // if(sector == EUpper && pos >= fMinLower) return 0;

    int sector = Sector(angle, family);
    REAL MaxUpper = MaxUpperAngle();;
    REAL MinLower = MinLowerAngle();
    if(angle >= MaxUpper && family == EUpper) {
        angle -= 2.*M_PI;
        refangle = minangle;
        std::cout << "Upper Adjusted angle to " << angle << " refangle " << refangle << std::endl;
    }

    if(angle <= MinLower && family == ELower) {
        angle += 2.*M_PI;
        refangle = maxangle;
        std::cout << "Lower Adjusted angle to " << angle << " refangle " << refangle << std::endl;
    }

    REAL deltheta = fAngles[sector+1]-fAngles[sector];
    REAL r = R(x);
    REAL ksi = Ksi(r);
    TPZFNMatrix<4,REAL> GradXInvT= {
        {cos(angle)/fRadius,-2.*sin(angle)/(ksi*deltheta*fRadius)},
        {sin(angle)/fRadius,2.*cos(angle)/(ksi*deltheta*fRadius)}
    };
    // GradXInvT.Print("GradXInvT");
    TPZSBFemVolume *vol = dynamic_cast<TPZSBFemVolume *>(fElementGroup->GetElGroup()[sector]);
    auto &PhiSBFem = vol->GetPhi();
    TPZManVector<REAL,1> eta = {2.*(angle - fAngles[sector])/(fAngles[sector+1]-fAngles[sector]) -1.};
    auto intel = fSkeletonElements[sector];
#ifdef PZDEBUG
    {
        if(!intel) {
            DebugStop();
        }
        TPZGeoEl *gel = intel->Reference();
        TPZManVector<REAL,3> locx(3,0.);
        gel->X(eta,locx);
        REAL diff = (locx[0]-fRadius*cos(angle))*(locx[0]-fRadius*cos(angle)) +
                    (locx[1]-fRadius*sin(angle))*(locx[1]-fRadius*sin(angle));
        if(diff > 1.e-6) {
            std::cout << "Computed point " << locx[0] << " " << locx[1] << " expected point " << fRadius*cos(angle) << " " << fRadius*sin(angle) << " diff " << diff << std::endl;
            std::cout << "Angle " << angle << " eta " << eta[0] << std::endl;
            std::cout << "Sector " << sector << " angle range " << fAngles[sector] << " " << fAngles[sector+1] << std::endl;
            std::cout << "middle angle " << 0.5*(fAngles[sector]+fAngles[sector+1]) << " deltheta " << deltheta <<std::endl;
            std::cout << "fRadius " << fRadius << std::endl;
            std::cout << "Norm(x) " << sqrt(x[0]*x[0]+x[1]*x[1]) << std::endl;
            DebugStop();
        }
    }
#endif
    int nshape = intel->NShapeF();
    TPZFMatrix<REAL> phiskel(nshape,1,0.), dphiskel(1,nshape,0.);
    if(eta[0] > 1) {
        TPZManVector<REAL,1> eta1 = {1.0};
        intel->Shape(eta1, phiskel, dphiskel);
        for(int ishape = 0; ishape<nshape; ishape++) {
            phiskel(ishape,0) += (eta[0] - 1.)*dphiskel(0,ishape);
        }   
    } else if (eta[0] < -1) {
        TPZManVector<REAL,1> eta1 = {-1.0};
        intel->Shape(eta1, phiskel, dphiskel);
        for(int ishape = 0; ishape<nshape; ishape++) {
            phiskel(ishape,0) += (eta[0] + 1.)*dphiskel(0,ishape);
        }
    } else {
        intel->Shape(eta, phiskel, dphiskel);
    }
    intel->Shape(eta, phiskel, dphiskel);
    int64_t neigen = fEigenvectorIndexes.size();
    phi.resize(neigen);
    dphi.resize(neigen);
    TPZVec<CSTATE> eigenvalues(neigen);
    
    auto computeEnrichmentMode = [this, &PhiSBFem, ksi, &phiskel, &dphiskel, nshape, &GradXInvT]
                                  (int ieigen, TPZManVector<std::complex<REAL>,2> &sol, TPZManVector<std::complex<REAL>,4> &dsolxy) {
        int eigenindex = fEigenvectorIndexes[ieigen];
        std::complex<REAL> lambda = -fElementGroup->EigenValues()[eigenindex];
        // std::cout << " lambda[" << eigenindex << "] = " << lambda << std::endl;
        TPZManVector<std::complex<REAL>,4> dsolksieta(4,0.);
        std::complex<REAL> explambda = std::pow(ksi,lambda);
        std::complex<REAL> dexplambda_dksi = lambda*std::pow(ksi,lambda-1.);
        
        for(int ishape = 0; ishape<nshape; ishape++) {
            std::complex<REAL> phivalx = PhiSBFem(2*ishape, eigenindex)*explambda;
            std::complex<REAL> phivaly = PhiSBFem(2*ishape+1, eigenindex)*explambda;
            sol[0] += phiskel(ishape,0)*phivalx;
            // derivative of sol[0] respect to ksi
            dsolksieta[0] += PhiSBFem(2*ishape, eigenindex)*dexplambda_dksi*phiskel(ishape,0);
            // derivative of sol[0] respect to eta
            dsolksieta[1] += PhiSBFem(2*ishape, eigenindex)*explambda*dphiskel(0,ishape);
            sol[1] += phiskel(ishape,0)*phivaly;
            // derivative of sol[1] respect to ksi
            dsolksieta[2] += PhiSBFem(2*ishape+1, eigenindex)*dexplambda_dksi*phiskel(ishape,0);
            // derivative of sol[1] respect to eta
            dsolksieta[3] += PhiSBFem(2*ishape+1, eigenindex)*explambda*dphiskel(0,ishape);
        }
        // std::cout << " sol = (" << sol[0] << "," << sol[1] << ")\n";
        // std::cout << " dsolksieta = (" << dsolksieta[0] << "," << dsolksieta[1] << "," << dsolksieta[2] << "," << dsolksieta[3] << ")\n";
        // compute derivatives respect to x and y
        dsolxy[0] = GradXInvT(0,0)*dsolksieta[0] + GradXInvT(0,1)*dsolksieta[1];
        dsolxy[1] = GradXInvT(1,0)*dsolksieta[0] + GradXInvT(1,1)*dsolksieta[1];
        dsolxy[2] = GradXInvT(0,0)*dsolksieta[2] + GradXInvT(0,1)*dsolksieta[3];
        dsolxy[3] = GradXInvT(1,0)*dsolksieta[2] + GradXInvT(1,1)*dsolksieta[3];
    };
    
    for(int ieigen = 0; ieigen<neigen; ieigen++) {
        eigenvalues[ieigen] = fElementGroup->EigenValues()[fEigenvectorIndexes[ieigen]];
        TPZManVector<std::complex<REAL>,2> sol(2,0.);
        TPZManVector<std::complex<REAL>,4> dsolxy(4,0.);
        computeEnrichmentMode(ieigen, sol, dsolxy);
        
        if(fPartOfEigenvector[ieigen] == EPartOfEigenvector::EOnlyReal) {
            phi[ieigen] = {sol[0].real(), sol[1].real()};
            dphi[ieigen] = {dsolxy[0].real(), dsolxy[1].real(), dsolxy[2].real(), dsolxy[3].real()};
        } else if (fPartOfEigenvector[ieigen] == EPartOfEigenvector::ERealPart) {
            // the eigenvector and eigenvalue are complex conjugates, the multiplying coefficients are also complex conjugates
            // the sum of both gives the real part (the sum the multiplying coefficients are real)
            // compute also the imaginary part (the difference of the multiplying coefficients are imaginary)
            TPZManVector<std::complex<REAL>,2> sol2(2,0.);
            TPZManVector<std::complex<REAL>,4> dsolxy2(4,0.);
            computeEnrichmentMode(ieigen+1, sol2, dsolxy2);
            phi[ieigen] = {(sol[0]+sol2[0]).real(), (sol[1]+sol2[1]).real()};
            dphi[ieigen] = {(dsolxy[0]+dsolxy2[0]).real(), (dsolxy[1]+dsolxy2[1]).real(), 
                (dsolxy[2]+dsolxy2[2]).real(), (dsolxy[3]+dsolxy2[3]).real()};
            phi[ieigen+1] = {(sol[0]-sol2[0]).imag(), (sol[1]-sol2[1]).imag()};
            dphi[ieigen+1] = {(dsolxy[0]-dsolxy2[0]).imag(), (dsolxy[1]-dsolxy2[1]).imag(), 
                (dsolxy[2]-dsolxy2[2]).imag(), (dsolxy[3]-dsolxy2[3]).imag()};
            ieigen++;
        } else if (fPartOfEigenvector[ieigen] == EPartOfEigenvector::EImagPart) {
        } else {
            DebugStop();
        }
    }
    // if(fabs(angle - refangle) > 1.e-6) {
    //     ExtendEnrichmentFunction(family, x, phi, dphi);
    // }
    // Compute the enrichment function values
    // EnrichmentFunction(family, x, phi, dphi);
}


/// @brief extend the enrichment function for angles outside the provided range
void TPZGFemEnrichment::ExtendEnrichmentFunction(MSectorFamily family,
                                 const std::vector<REAL> &x,
                                 std::vector<std::array<REAL, 2>> &phi,
                                 std::vector<std::array<REAL ,4>> &dphixy) {
    // Use the static helper methods declared in the header.
    REAL angle = std::atan2(x[1],x[0]);

    REAL minangle = fAngles[0];
    REAL maxangle = fAngles[fAngles.size()-1];

    while(angle < minangle) {
        angle += 2.*M_PI;
    }
    while(angle > maxangle) {
        angle -= 2.*M_PI;
    }
    REAL refangle = 0.;;
    REAL delta_angle = 0.;
    if(angle > minangle+M_PI && family == EUpper) {
        angle -= 2.*M_PI;
        refangle = minangle;
        delta_angle = angle-minangle; // negative angle
    }
    TPZFNMatrix<4,REAL> rotate(2,2);

    if(angle < maxangle -M_PI && family == ELower) {
        angle += 2.*M_PI;
        refangle = maxangle;
        delta_angle = angle - maxangle; // positive angle
    }

    REAL r = R(x);
    std::vector<REAL> locx = {r*std::cos(refangle), r*std::sin(refangle)};
    ComputeEnrichmentWithDerivatives(family, locx, phi, dphixy);
    int nphi = phi.size();
    for(int i=0; i<nphi; i++) {
        // bring phi and dphixy to the polar coordinate system aligned with refangle
        TPZGFemEnrichment::RotateVector(refangle, phi[i]);
        TPZGFemEnrichment::RotateTensor(refangle, dphixy[i]);
        // extend phi using the derivatives the theta direction
        phi[i][0] += delta_angle * dphixy[i][1];
        phi[i][1] += delta_angle * dphixy[i][3];
        // bring back to the global coordinate system
        TPZGFemEnrichment::RotateVector(-refangle-delta_angle, phi[i]);
        TPZGFemEnrichment::RotateTensor(-refangle-delta_angle, dphixy[i]);
    }
}

// Static helper implementations
void TPZGFemEnrichment::ComputeRotationMatrix(REAL angle, TPZFNMatrix<4,REAL> &rotmat) {
    rotmat.Resize(2,2);
    REAL c = cos(angle);
    REAL s = sin(angle);
    rotmat(0,0) = c;
    rotmat(0,1) = -s;
    rotmat(1,0) = s;
    rotmat(1,1) = c;
}

void TPZGFemEnrichment::RotateVector(REAL angle, std::array<REAL,2> &vec) {
    TPZFNMatrix<4,REAL> rotmat(2,2);
    TPZGFemEnrichment::ComputeRotationMatrix(angle, rotmat);
    std::array<REAL,2> vecorig(vec);
    vec[0] = rotmat(0,0)*vecorig[0]+rotmat(0,1)*vecorig[1];
    vec[1] = rotmat(1,0)*vecorig[0]+rotmat(1,1)*vecorig[1];
}

void TPZGFemEnrichment::RotateTensor(REAL angle, std::array<REAL,4> &tensor) {
    TPZFNMatrix<4,REAL> rotmat(2,2), tensorm(2,2), tenorrot(2,2);
    tensorm(0,0) = tensor[0];
    tensorm(1,0) = tensor[2];
    tensorm(0,1) = tensor[1];
    tensorm(1,1) = tensor[3];
    TPZGFemEnrichment::ComputeRotationMatrix(angle, rotmat);
    for(int i=0; i<2; i++) {
        for(int j=0; j<2; j++) {
            tenorrot(i,j) = 0.;
            for(int k=0; k<2; k++) {
                tenorrot(i,j) += rotmat(i,k)*tensorm(k,j);
            }
        }
    }
    for(int i=0; i<2; i++) {
        for(int j=0; j<2; j++) {
            tensorm(i,j) = 0.;
            for(int k=0; k<2; k++) {
                tensorm(i,j) += tenorrot(i,k)*rotmat(j,k);
            }
        }
    }
    tensor[0] = tensorm(0,0);
    tensor[1] = tensorm(0,1);
    tensor[2] = tensorm(1,0);
    tensor[3] = tensorm(1,1);
}
