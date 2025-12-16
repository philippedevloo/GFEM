#include "pzlog.h"
#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "TPZGmshReader.h"
#include "TPZBuildSBFem.h"

#include "TPZBndCondT.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "Elasticity/TPZElasticity2D.h"
#include "TPZAnalyticSolution.h"

#include "tpzarc3d.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"


#include <nlohmann/json.hpp>

using json = nlohmann::json;


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
void ComputeSingularMode(TPZCompMesh &cmesh);

/// @brief return the index of the eigenvectors corresponding to singular modes
std::list<int64_t> GetSingularModes(TPZCompMesh &cmesh);

/// @brief Load an eigenmode from a file into the computational mesh
void LoadEigenMode(TPZCompMesh &cmesh, int64_t modeindex);

/// @brief Plot the solution of the computational mesh
/// @param cmesh reference to the computational mesh
/// @param output json object with the output data
void PlotSolution(const std::string &rootname, TPZCompMesh &cmesh, int refstep, json &input);

std::string gProblem;
std::map<std::string,TLaplaceExample1> gExact1;
std::map<std::string,TElasticity2DAnalytic> gExact2;


/// @brief Add an exact solution to the map
void AddExact(const std::string &problemtype, const std::string &exactname) {
    if(problemtype == "Scalar") {
        if(gExact1.find(exactname) != gExact1.end()) return;
        TLaplaceExample1::EExactSol enumfunc = TLaplaceExample1::StringToExactSol(exactname);
        if(enumfunc != TLaplaceExample1::ENone) {
            gExact1[exactname].fExact = enumfunc;
        } else {
            DebugStop();
        }
    } else if (problemtype == "Elastic") {
        if(gExact2.find(exactname) != gExact2.end()) return;
        TElasticity2DAnalytic::EDefState enumfunc = TElasticity2DAnalytic::StringToExactSol(exactname);
        if(enumfunc != TElasticity2DAnalytic::ENone) {
            gExact2[exactname].fProblemType = enumfunc;
        } else {
            DebugStop();
        }
    }

}



int main() {


#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif


    // TPZFNMatrix<9,REAL> coord(3,3,0.);
    // coord(0,0) = 1.;
    // coord(1,1) = 1.;
    // coord(0,2) = M_SQRT1_2;
    // coord(1,2) = M_SQRT1_2;

    // pzgeom::TPZArc3D arc3d(coord);

    // Initializing uniform refinements for reference elements
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    
  	using json = nlohmann::json;
    std::cout << "Reading input file\n";
//    std::string filename = "steklovScalar.json";
    std::string filename = "CrackDef.json";
//    std::string filename = "SqrtSingular.json";
//    std::string filename = "elasticsmooth.json";
//    std::string filename = "SqrtSingularElastic.json";
#ifdef MACOSX22
    filename = "../" + filename;
#endif
	// Read file
	std::ifstream file(filename);
	if (!file){
		std::cout << "\nCouldn't find file \"" << filename << "\""<< std::endl;
		DebugStop();
	}

	// Parse json
	json input;
	// file >> input;
	input = json::parse(file,nullptr,true,true); // to ignore comments in json file

    TPZAutoPointer<TPZGeoMesh> gmesh = CreateGeometry(input);

    {
        std::ofstream out("gmesh.txt");
        gmesh->Print(out);

        std::ofstream out2("gmesh.vtk"); 
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out2);
    }

    TPZBuildSBFem builder(gmesh);
    TPZCompMesh *cmeshH1ptr = nullptr;
    cmeshH1ptr = CreateSBFemSpace(builder, input);

    ComputeSingularMode(*cmeshH1ptr);
    {
        std::ofstream out("cmeshH1.txt");
        cmeshH1ptr->Print(out);
    }
    std::list<int64_t> singularmodes = GetSingularModes(*cmeshH1ptr);
    int refstep = 0;
    for(auto modeindex : singularmodes) {
        LoadEigenMode(*cmeshH1ptr, modeindex);
        PlotSolution("EigenMode", *cmeshH1ptr, refstep, input);
        refstep++;
    }
    {
        std::ofstream out("cmeshH1.txt");
        cmeshH1ptr->Print(out);
    }

    return 0;
}

#include <list>
#include <tuple>
#include "TPZGeoElement.h"
#include "TPZRefLinear.h"
#include "pzgeotriangle.h"
#include "tpzgeoblend.h"

TPZGeoMesh *CreateGeometry(json &input) {

    if(input.find("geometry") == input.end()){
		DebugStop();
    }
    auto geom = input["geometry"];
    REAL radius = geom["radius"];
    std::list<std::tuple<REAL,REAL,int>> sectors;
    REAL prevangle = 0.;
    for(auto &item : geom["angles"]) {
        sectors.push_back({prevangle,item["angle"],item["material"]});
        prevangle = item["angle"];
    }
    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(2);
    // we assume the centernode is at 0,0
    int64_t centerindex = gmesh->NodeVec().AllocateNewElement();
    TPZManVector x(3,0.);
    gmesh->NodeVec()[centerindex].Initialize(x,*gmesh);
    int pointmatid = 0;
    if(input.find("point_singularity") != input.end()) {
        pointmatid = input["point_singularity"];
    } else {
        DebugStop();
    }
    {
        TPZManVector<int64_t> nodeindex = {centerindex};
        int64_t index;
        gmesh->CreateGeoElement(EPoint,nodeindex,pointmatid,index);
    }


    // generate the arc nodes and arc elements
    int arcmatid = input["arcmatid"];
    int64_t ncirclenodes = sectors.size()+1;
    x[0] = radius;
    TPZManVector<int64_t> nodeindexes(ncirclenodes,-1);
    nodeindexes[0] = gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[nodeindexes[0]].Initialize(x,*gmesh);
    int count = 1;
    for(auto &it : sectors) {
        nodeindexes[count] = gmesh->NodeVec().AllocateNewElement();
        REAL prevangle = std::get<0>(it);
        REAL angle = std::get<1>(it);
        x[0] = radius*std::cos(angle);
        x[1] = radius*std::sin(angle);
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
            ax1[0] /= radius;
            ax1[1] /= radius;
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
            gel->Geom().SetAttributes(gmesh, center, rotationmatrix, radius, openingangle);
        }
        // create the triangular element
        TPZManVector<int64_t> triindexes(lineindexes);
        triindexes[2] = centerindex;
        int64_t index;
        gmesh->CreateGeoBlendElement(ETriangle,triindexes,matid,index);
        prevangle = angle;
        count++;
    }

    gmesh->BuildConnectivity();
    

    return gmesh;
}

/// @brief Create a computational mesh based on the SBFem hybrid H1 space
/// @param gmesh autopointer to the geometric mesh
/// @param input json object with the input data
/// @return pointer to the computational mesh
TPZCompMesh *CreateSBFemSpace(TPZBuildSBFem &builder, json &input) {
    TPZAutoPointer<TPZGeoMesh> gmesh = builder.GetGMesh();
    int dim = gmesh->Dimension();
    ConfigureSBFemBuilder(builder, input);
    int64_t singular_partition = -1;
    int64_t singular_element = -1;
    if(input.find("point_singularity") != input.end()) {
        // two steps :
        // 1 create standard elements for all elements that do not touch the singularity
        int singular_matid = input["point_singularity"];
        TPZStack<int64_t> volels;
        std::map<int,int> matidtranslation = builder.GetMatIdTranslation();
        
        int64_t singular_node = -1;
        TPZGeoEl *singular_element = 0;
        int64_t nel = gmesh->NElements();
        std::set<int64_t> singular_group;
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel || gel->MaterialId() != singular_matid) continue;
            singular_node = gel->NodeIndex(0);
            singular_element = gel;
            break;
        }
        if(singular_node == -1) DebugStop();
        TPZGeoElSide singularside (singular_element);
        std::set<int64_t> to_invert;
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
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel || gel->Dimension() != dim) continue;
            if(singular_group.find(el) != singular_group.end()) continue;
            volels.Push(el);
        }
        builder.CreateElementCenterNodes(volels);
        singular_partition = builder.AddPartition(singular_group, singular_node);
        builder.AddSkeletonElements();
        // 2 create collapsed elements towards the singular point
    } else {
        builder.StandardConfiguration();
    }
    // builder.SetNRefSkeleton(0);
    int nref = builder.GetNRefSkeleton();
    builder.DivideSkeleton(nref);
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
    
    // std::set<int> matidsskel = builder.GetBoundaryMatIds();
    // matidsskel.insert(builder.GetSkeletonMatid());
    // cmesh.SetDefaultOrder(builder.GetSkeletonPOrder());
    // cmesh.SetAllCreateFunctionsContinuous();
    // cmesh.AutoBuild(matidsskel);
    // builder.CreateVolumetricElements(cmesh);
    // builder.CreateElementGroups(cmesh);

    return cmeshptr;
}

/// @brief configure the SBFem builder according to the json input
/// @param builder reference to the SBFem builder
/// @param input json object with the input data
void ConfigureSBFemBuilder (TPZBuildSBFem &builder, json &input) {
    int skeletonmatid = input["skeletonMatid"];
    builder.SetSkeletonMatid(skeletonmatid);


    std::map<int,int> matidtranslation;
    std::set<int> sbfemids;
    for (const auto &item : input["matid_translation"]) {
        int from = item["from"];
        int to = item["to"];
        matidtranslation[from] = to;
        sbfemids.insert(to);
    }
    builder.SetMatIdTranslation(matidtranslation);
    std::set<int> boundmatids;
    for (const auto &item : input["boundary_conditions"]) {
        int matid = item["matid"];
        if(sbfemids.find(matid) != sbfemids.end()) continue;
        boundmatids.insert(matid);
    }
    builder.SetBoundaryMatIds(boundmatids);
    int skeletonporder = input["skeleton porder"];
    builder.SetSkeletonPOrder(skeletonporder);
    int nref = input["nref_skeleton"];
    builder.SetNRefSkeleton(nref);
}

void InsertMaterialObjectsH1(TPZCompMesh &cmesh, json &input)
{
    int dim = cmesh.Dimension();
    TPZMaterialT<STATE> *mat = nullptr;
    int nstate = 0;
    std::string problemtype = input["problem_type"];

    for (const auto &item : input["materials"]) {
        int matid = item["matid"];
        if(problemtype == "Scalar") {
            nstate = 1;
            REAL perm = item["permeability"];
            auto darcy = new TPZDarcyFlow(matid, dim);
            //        auto darcy = new TPZDarcyFlow(matid, dim);
            darcy->SetConstantPermeability(perm);
            if(item.find("exactsolution") != item.end()) {
                std::string name = item["exactsolution"];
                AddExact(problemtype,name);
                darcy->SetForcingFunction(gExact1[name].ForceFunc(),5);
                darcy->SetExactSol(gExact1[name].ExactSolution(),5);
            }
            cmesh.InsertMaterialObject(darcy);
            if(!mat) mat = darcy;
        } else if (problemtype == "Elastic") {
            nstate = 2;
            REAL young = item["young"];
            REAL poisson = item["poisson"];
            REAL fx(0.), fy(0.);
            int planestress = 1;
            TPZElasticity2D *matelast = new TPZElasticity2D(matid,young,poisson,fx,fy,planestress);
            if(item.find("exactsolution") != item.end()) {
                std::string name = item["exactsolution"];
                AddExact(problemtype,name);
                matelast->SetForcingFunction(gExact2[name].ForceFunc(),5);
                matelast->SetExactSol(gExact2[name].ExactSolution(),5);
                planestress = gExact2[name].fPlaneStress;
                if(planestress) matelast->SetPlaneStress();
                else matelast->SetPlaneStrain();
                gExact2[name].gE = young;
                gExact2[name].gPoisson = poisson;
            }
            if(!mat) mat = matelast;
            cmesh.InsertMaterialObject(matelast);
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
        if(item.find("exactsolution") != item.end()) {
            std::string name = item["exactsolution"];
            AddExact(problemtype,name);
            if(nstate == 1) {
                bc->SetForcingFunctionBC(gExact1[name].ExactSolution(),5);
            } else if(nstate == 2) {
                bc->SetForcingFunctionBC(gExact2[name].ExactSolution(),5);
            }
        }
        cmesh.InsertMaterialObject(bc);
    }

    int skeletonMatId = input["skeletonMatid"];
    TPZNullMaterial<STATE> *skelmat = new TPZNullMaterial<STATE>(skeletonMatId, dim,nstate);
    cmesh.InsertMaterialObject(skelmat);
}

#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZLinearAnalysis.h"

/// @brief Compute the stiffness matrix and load the singular mode as the mesh solution
void ComputeSingularMode(TPZCompMesh &cmesh) {
    TPZFStructMatrix<STATE> strmat(&cmesh);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    TPZLinearAnalysis an(&cmesh,RenumType::ENone);
    an.SetStructuralMatrix(strmat);
    an.SetSolver(step);
    TPZFMatrix<STATE> rhs;
    an.Assemble();
}

#include "TPZSBFemElementGroup.h"
/// @brief load the singular mode as the mesh solution
std::list<int64_t> GetSingularModes(TPZCompMesh &cmesh) {

    std::list<int64_t> singularmodes;
    REAL range[] = {0.01,0.9};
    /// find the SBFEmElementGroup
    TPZSBFemElementGroup *elgr = 0;
    int64_t nel = cmesh.NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh.Element(el);
        if(!cel) continue;
        elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if(elgr) break;
    }
    if(!elgr) {
        DebugStop();
    }
    // find an eigenvalue larger than zero and smaller than 1
    TPZManVector<std::complex<REAL>,10> eigenvalues;
    eigenvalues = elgr->EigenValues();
    int64_t neigen = eigenvalues.size();
    int64_t selected = -1;
    for (int64_t ieigen = 0; ieigen<neigen; ieigen++) {
        REAL realpart = eigenvalues[ieigen].real();
        if(abs(realpart) > range[0] && abs(realpart) < range[1]) {
            singularmodes.push_back(ieigen);
        }
    }
    std::cout << "Found " << singularmodes.size() << " singular modes in the range (" << range[0] << "," << range[1] << ")\n";
    std::cout << "Selected Eigenvalues:\n";
    for (auto modeindex : singularmodes) {
        std::cout << "  " << modeindex << " : " << eigenvalues[modeindex] << "\n";
    }
    return singularmodes;
}

#include "TPZVTKGenerator.h"
/// @brief Plot the solution of the computational mesh
/// @param cmesh reference to the computational mesh
/// @param output json object with the output data
void PlotSolution(const std::string &rootname, TPZCompMesh &cmesh, int refstep, json &output) {
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
void LoadEigenMode(TPZCompMesh &cmesh, int64_t modeindex) {
    /// find the SBFEmElementGroup
    TPZSBFemElementGroup *elgr = 0;
    int64_t nel = cmesh.NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh.Element(el);
        if(!cel) continue;
        elgr = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if(elgr) break;
    }
    if(!elgr) {
        DebugStop();
    }
    elgr->LoadEigenVector(modeindex);
}

