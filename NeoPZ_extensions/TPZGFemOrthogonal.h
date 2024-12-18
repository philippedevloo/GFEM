#ifndef TPZGFemOrthogonal_H
#define TPZGFemOrthogonal_H

#include "TPZMultiphysicsCompMesh.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzmatred.h"

class TPZGFemOrthogonal {
public:
    // Constructor
    TPZGFemOrthogonal(TPZMultiphysicsCompMesh *multiphysicsCompMesh, TPZMatrix<STATE> *globalMatrix) : fMultiphysicsCompMesh(multiphysicsCompMesh), fGlobalMatrix(globalMatrix){
        fMultiphysicsCompMesh->ComputeNodElCon();
    }

    // Destructor
    ~TPZGFemOrthogonal() {
        // Clean up resources here
    }

    // Data structure linking a connect index to H1 connect indexes
    struct TNodePatch {
        int64_t fGFemConnectIndex;

        // connect indexes of the H1 meshes. First the mesh 0 connects, then the mesh 1 connects
        TPZVec<std::pair<int64_t, int64_t>> fH1ConnectIndexes;

        // smallest eigenvalue of the patch before orthogonalization
        REAL fOriginalSmallestEigenvalue;

        // smallest eigenvalue of the patch after orthogonalization
        REAL fOrthogonalizedSmallestEigenvalue;
    };

public:

    // Orthogonalize Connects
    void OrthogonalizeConnects();

    // Verify is the equations are orthogonalized
    void VerifyOrthogonality();

    // return the connect index with largest eigenvalue ratio
    int64_t LargestEigenvalueRatio();

    // represent the orthogonalization graphically
    void DrawOrthogonalization(int gfemconnectindex, const std::string &filename);
private:
    // Initialize the data structure fNodePatches
    void BuildNodePatches();

    // build the correspondence between the mesh 0 connects and mesh 1 connects
    void BuildConnectCorrespondence(std::map<int64_t, int64_t> &mesh0tomesh1);

    /// build the node connectivity of a GFem node
    void BuildConnectConnectivity(int64_t *firstel, int64_t *lastel, TPZVec<int64_t> &connects);

    /// Compute the element graph as a function of the connect indexes
    void ComputeElementGraph(TPZStack<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex);

    // Compute the orthogonalization matrix for a given connect index
    void AnaliseAndReduceConnect(TNodePatch &patch, TPZMatRed<STATE,TPZFMatrix<STATE>> &matred);

    // Compute the stiffness matrix of a patch
    void ComputePatchMatrix(TNodePatch &patch, TPZFMatrix<STATE> &patchmatrix);

    // For each connect index of the GFem mesh, store the connect indexes of the H1 meshes
    TPZManVector<TNodePatch> fNodePatches;

    // Pointer to a TPZMultiphysicsCompMesh object
    TPZMultiphysicsCompMesh *fMultiphysicsCompMesh;

    // Pointer to an assembled global matrix object
    TPZMatrix<STATE> *fGlobalMatrix;

};

#endif // TPZGFemOrthogonal_H