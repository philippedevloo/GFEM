#include "TPZGFemOrthogonal.h"

#include "TPZRenumbering.h"
#include "TPZGFemCompMesh.h"



void TPZGFemOrthogonal::BuildNodePatches()
{
    // BuildNodePatches implementation
    int64_t nconnects = fMultiphysicsCompMesh->NConnects();
    int64_t nelem = fMultiphysicsCompMesh->NElements();
    TPZVec<TPZCompMesh *> &meshvec = fMultiphysicsCompMesh->MeshVector();
    TPZGFemCompMesh *gfemmesh = dynamic_cast<TPZGFemCompMesh *>(meshvec[2]);
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> nodtoelgraph;
    TPZVec<int64_t> elgraphindex, nodtoelgraphindex;
    ComputeElementGraph(elgraph, elgraphindex);
    TPZRenumbering renumber(nelem,nconnects);
    renumber.NodeToElGraph(elgraph, elgraphindex, nodtoelgraph, nodtoelgraphindex);
    std::map<int64_t, int64_t> mesh0tomesh1;
    BuildConnectCorrespondence(mesh0tomesh1);
    int64_t firstgfemconnect = meshvec[0]->NConnects()+meshvec[1]->NConnects();
    int64_t nconnectsgfem = meshvec[2]->NConnects();
    fNodePatches.resize(gfemmesh->fShapeFunctionMap.size());
    int64_t firstconnect1 = meshvec[0]->NConnects();
    int64_t firstconnect2 = meshvec[0]->NConnects()+meshvec[1]->NConnects();
    if(0)
    {
        int64_t lastcon = fMultiphysicsCompMesh->NConnects()-1;
        TPZConnect &c = fMultiphysicsCompMesh->ConnectVec()[lastcon];
        c.Print(*fMultiphysicsCompMesh,std::cout);
        std::cout << "nodtoelgraphindex " << nodtoelgraphindex[lastcon] << " " << nodtoelgraphindex[lastcon+1] << "\n";
        std::cout << "elements connected to last connect\n";
        for(int64_t el = nodtoelgraphindex[lastcon]; el < nodtoelgraphindex[lastcon+1]; el++) {
            TPZCompEl *cel = fMultiphysicsCompMesh->Element(nodtoelgraph[el]);
            if(!cel) DebugStop();
            cel->Print(std::cout);
        }
    }
    int count = 0;
    for (int64_t ic = 0; ic < nconnectsgfem; ic++) {
        if(gfemmesh->fShapeFunctionMap.find(ic) == gfemmesh->fShapeFunctionMap.end()) {
            continue;
        }
        fNodePatches[count].fGFemConnectIndex = firstgfemconnect+ic;
        int64_t firstel = nodtoelgraphindex[ic+firstgfemconnect];
        int64_t lastel = nodtoelgraphindex[firstgfemconnect+ic+1];
        TPZVec<int64_t> connects;
        {
            std::ofstream out("elgraph.txt");
            for(int64_t el = firstel; el < lastel; el++) {
                TPZCompEl *cel = fMultiphysicsCompMesh->Element(nodtoelgraph[el]);
                if(!cel) DebugStop();
                cel->Print(out);
            }
        }
        BuildConnectConnectivity(&nodtoelgraph[firstel], &nodtoelgraph[lastel], connects);
        int nconnect0 = 0;
        int nconnect1 = 0;
        int nconnect2 = 0;
        for (auto it : connects) {
            if(it < firstconnect1) nconnect0++;
            else if(it < firstconnect2) nconnect1++;
            else nconnect2++;
        }
        if(nconnect0 != nconnect1 || nconnect2 != nconnect1) {
            DebugStop();
        }
        fNodePatches[count].fH1ConnectIndexes.resize(nconnect0);
        nconnect0 = 0;
        for (auto it : connects) {
            if(it < firstconnect1) {
                if(mesh0tomesh1.find(it) == mesh0tomesh1.end()) {
                    DebugStop();
                }
                int64_t secondconnect = mesh0tomesh1[it];
                fNodePatches[count].fH1ConnectIndexes[nconnect0].first = it;
                fNodePatches[count].fH1ConnectIndexes[nconnect0].second = secondconnect;
                nconnect0++;
            }
        }
        count++;
    }

}

    // build the correspondence between the mesh 0 connects and mesh 1 connects
void TPZGFemOrthogonal::BuildConnectCorrespondence(std::map<int64_t, int64_t> &mesh0tomesh1)
{
    // BuildConnectCorrespondence implementation
    mesh0tomesh1.clear();
    int64_t nelem = fMultiphysicsCompMesh->NElements();
    for (int64_t el = 0; el < nelem; el++) {
        TPZCompEl *cel = fMultiphysicsCompMesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mphys) {
            DebugStop();
        }
        TPZCompEl *cel2 = mphys->Element(2);
        if (!cel2) {
            continue;
        }
        TPZCompEl *cel0 = mphys->Element(0);
        TPZCompEl *cel1 = mphys->Element(1);
        if (!cel0 || !cel1) {
            DebugStop();
        }
        int64_t nconnects0 = cel0->NConnects();
        int64_t nconnects1 = cel1->NConnects();
        int64_t nconnects = cel->NConnects();
        for (int64_t ic = 0; ic < nconnects0; ic++) {
            int64_t c0index = cel->ConnectIndex(ic);
            int64_t c1index = cel->ConnectIndex(ic+nconnects0);
            if(mesh0tomesh1.find(c0index) != mesh0tomesh1.end()) {
                if(mesh0tomesh1[c0index] != c1index) {
                    DebugStop();
                }
                continue;
            }
            mesh0tomesh1[c0index] = c1index;
        }
    }
}

    /// build the node connectivity of a GFem node
void TPZGFemOrthogonal::BuildConnectConnectivity(int64_t *firstel, int64_t *lastel, TPZVec<int64_t> &connects)
{
    // BuildNodeConnectivity implementation
    // number of elements connected to each connect
    std::map<int64_t, int> nodelcon;
    std::map<int64_t, int64_t> mesh0tomesh1;
    std::map<int64_t, int64_t> mesh1tomesh2;
    int64_t nel = lastel-firstel;
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fMultiphysicsCompMesh->Element(*(firstel+el));
        if (!cel) {
            continue;
        }
        int nnodes = cel->NConnects();
        if(nnodes%3) DebugStop();
        for (int ic = 0; ic < nnodes/3; ic++) {
            int64_t cindex = cel->ConnectIndex(ic);
            mesh0tomesh1[cindex] = cel->ConnectIndex(ic+nnodes/3);
            mesh1tomesh2[cindex] = cel->ConnectIndex(ic+2*nnodes/3);
            if(nodelcon.find(cindex) == nodelcon.end()) {
                nodelcon[cindex] = 1;
            } else {
                nodelcon[cindex]++;
            }
        }
    }
    int numconnects = 0;
    for(auto it : nodelcon) {
        int64_t cindex = it.first;
        int ncon = it.second;
        if(ncon == fMultiphysicsCompMesh->ConnectVec()[cindex].NElConnected()) numconnects++;
    }
    connects.resize(numconnects*3);
    numconnects = 0;
    for(auto it : nodelcon) {
        int64_t cindex = it.first;
        int ncon = it.second;
        if(ncon == fMultiphysicsCompMesh->ConnectVec()[cindex].NElConnected()) {
            connects[numconnects++] = cindex;
            connects[numconnects++] = mesh0tomesh1[cindex];
            connects[numconnects++] = mesh1tomesh2[cindex];
        }
    }
}

    /// Compute the element graph as a function of the connect indexes
void TPZGFemOrthogonal::ComputeElementGraph(TPZStack<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex)
{

    TPZStack<int64_t> connectstack;
    int64_t nelem = fMultiphysicsCompMesh->NElements();
    elgraphindex.Resize(nelem+1);
    elgraphindex[0] = 0;
    elgraph.Resize(0);
    int64_t curel=0;
    int64_t nindep = fMultiphysicsCompMesh->NConnects();
    for(int64_t el=0; el<nelem; el++) {
        TPZCompEl *cel = fMultiphysicsCompMesh->Element(el);
        if(!cel){
            elgraphindex[curel+1]=elgraph.NElements();
            curel++;
            continue;
        }
        

        
        connectstack.Resize(0);
        cel->BuildConnectList(connectstack);

        int64_t ncon = connectstack.NElements();
        for(int64_t in=0; in<ncon; in++) {
            int ic = connectstack[in];
            elgraph.Push(ic);
        }
        elgraphindex[curel+1]=elgraph.NElements();
        curel++;
    }
    elgraphindex.Resize(curel+1);

}

#include "pzmatred.h"
#include "pzstepsolver.h"

    // Compute the stiffness matrix of a patch
void TPZGFemOrthogonal::ComputePatchMatrix(TNodePatch &nodepatch, TPZFMatrix<STATE> &patchmatrix)
{
    int64_t neq_total = fMultiphysicsCompMesh->NEquations();

    int64_t gfemconnectindex = nodepatch.fGFemConnectIndex;
    TPZConnect &cgfem = fMultiphysicsCompMesh->ConnectVec()[gfemconnectindex];
    int64_t dim1 = cgfem.NShape()*cgfem.NState();
    int64_t dim0 = 0;
    // loop over the H1 connects
    for (auto &h1connectpair : nodepatch.fH1ConnectIndexes) {
        // orthogonalize the connect
        int64_t h1connectindex = h1connectpair.first;
        TPZConnect &ch1 = fMultiphysicsCompMesh->ConnectVec()[h1connectindex];
        dim0 += ch1.NShape()*ch1.NState();
    }
    patchmatrix.Redim(dim0+dim1,dim0+dim1);
    TPZManVector<int64_t,20> globeq(dim0+dim1,0);
    int64_t count = 0;
    for (auto &h1connectpair : nodepatch.fH1ConnectIndexes) {
        int64_t h1connectindex = h1connectpair.first;
        TPZConnect &ch1 = fMultiphysicsCompMesh->ConnectVec()[h1connectindex];
        int64_t nshape = ch1.NShape();
        int64_t nstate = ch1.NState();
        int blsize = nshape*nstate;
        int64_t seqnum = ch1.SequenceNumber();
        int64_t firsteq = fMultiphysicsCompMesh->Block().Position(seqnum);
        if(ch1.NShape()*ch1.NState() != fMultiphysicsCompMesh->Block().Size(seqnum)) DebugStop();
        bool active = true;
        if(ch1.HasDependency() || ch1.IsCondensed()) active = false;
        if(!active) {
            DebugStop();
        }
        if(firsteq+blsize > neq_total) DebugStop();
        for (int64_t is=0; is<blsize; is++) {
            globeq[count++] = firsteq+is;
        }
    }
    {
        int64_t nshape = cgfem.NShape();
        int64_t nstate = cgfem.NState();
        int64_t seqnum = cgfem.SequenceNumber();
        int blsize = fMultiphysicsCompMesh->Block().Size(seqnum);
        if(blsize != nshape*nstate) DebugStop();
        int64_t pos = fMultiphysicsCompMesh->Block().Position(seqnum);
        for (int64_t is=0; is<nshape*nstate; is++) {
            globeq[count++] = pos+is;
        }
    }
    if(count != dim0+dim1) DebugStop();
    for(int i = 0; i < dim0+dim1; i++) {
        int64_t eqi = globeq[i];
        for(int j = 0; j < dim0+dim1; j++) {
            int64_t eqj = globeq[j];
            patchmatrix.PutVal(i,j, fGlobalMatrix->GetVal(eqi,eqj));
        }
    }
}

// Orthogonalize Connects
void TPZGFemOrthogonal::OrthogonalizeConnects()
{
    int64_t neq_total = fMultiphysicsCompMesh->NEquations();
    // loop over the node patches
    for (auto &nodepatch : fNodePatches) {
        int64_t gfemconnectindex = nodepatch.fGFemConnectIndex;
        TPZConnect &cgfem = fMultiphysicsCompMesh->ConnectVec()[gfemconnectindex];
        int64_t dim1 = cgfem.NShape()*cgfem.NState();
        int64_t dim0 = 0;
        // loop over the H1 connects
        for (auto &h1connectpair : nodepatch.fH1ConnectIndexes) {
            // orthogonalize the connect
            int64_t h1connectindex = h1connectpair.first;
            TPZConnect &ch1 = fMultiphysicsCompMesh->ConnectVec()[h1connectindex];
            dim0 += ch1.NShape()*ch1.NState();
        }

        TPZMatRed<STATE, TPZFMatrix<STATE> > matred(dim0+dim1, dim0);
        TPZAutoPointer<TPZMatrix<STATE> > k00auto = new TPZFMatrix<STATE>(dim0,dim0,0.);
        for(int i =0; i<dim0; i++) for(int j=0; j<dim0; j++) k00auto->PutVal(i,j,i+j);
        TPZStepSolver<STATE> *step = new TPZStepSolver<STATE>(k00auto);
        step->SetDirect(ECholesky);
        TPZAutoPointer<TPZMatrixSolver<STATE> > autostep = step;
        matred.SetSolver(autostep);
        TPZFMatrix<STATE> patchmatrix(dim0+dim1,dim0+dim1,0.);
        ComputePatchMatrix(nodepatch, patchmatrix);
        for(int i = 0; i<dim0+dim1; i++) {
            for(int j = 0; j<dim0+dim1; j++) {
                matred(i,j) = patchmatrix.GetVal(i,j);
                if(std::isnan(patchmatrix.GetVal(i,j))) DebugStop();
            }
        }
        AnaliseAndReduceConnect(nodepatch, matred);
    }
}

void NormalizeDiagonal(TPZMatrix<STATE> &mat)
{
    int64_t dim = mat.Rows();
    TPZVec<REAL> diag(dim,0.);
    for(int64_t i = 0; i < dim; i++) {
        diag[i] = std::sqrt(mat.GetVal(i,i));
    }
    for(int64_t i = 0; i < dim; i++) {
        for(int64_t j = 0; j < dim; j++) {
            mat(i,j) /= diag[i]*diag[j];
        }
    }
}

#include <complex>
// Compute the orthogonalization matrix for a given connect index
void TPZGFemOrthogonal::AnaliseAndReduceConnect(TNodePatch &patch, TPZMatRed<STATE, TPZFMatrix<STATE>> &matred)
{
    int64_t dim = matred.Rows();
    int64_t dim0 = matred.Dim0();
    int64_t dim1 = matred.Dim1();
    // NormalizeDiagonal(matred);
    TPZFNMatrix<50,STATE> matlocal(matred);
    NormalizeDiagonal(matlocal);
    // matlocal.Print("matred",std::cout,EMathematicaInput);
    if(0)
    {
        TPZManVector<CSTATE,20> eigenvalues(dim);
        TPZFMatrix<CSTATE> eigenvectors(dim,dim);
        matlocal.SolveEigenProblem(eigenvalues,eigenvectors);
        {
            REAL small = std::norm(eigenvalues[0]);
            for(auto it : eigenvalues) {
                if(std::norm(it) < small) small = std::norm(it);
            }
            patch.fOriginalSmallestEigenvalue = small;
        }
    }
    {
        TPZFNMatrix<10,STATE> eigenvector(dim,1,0.);
        for(int i=0; i<dim; i++) eigenvector(i,0) = rand()*1./RAND_MAX;
        TPZFNMatrix<50,STATE> matinv;
        REAL norm = Norm(eigenvector);
        // matlocal.Print(std::cout);
        eigenvector = eigenvector * (1./norm);
        matlocal.Inverse(matinv,ECholesky);
        // matinv.Print(std::cout);
        for(int it = 0; it<100; it++) {
            eigenvector = matinv*eigenvector;
            norm = Norm(eigenvector);
            eigenvector = eigenvector * (1./norm);
            // std::cout << "norm = " << norm << std::endl;
        }
        patch.fOriginalSmallestEigenvalue = 1./norm;

    }
    matlocal = matred;
    TPZFMatrix<STATE> k11(dim1,dim1,0.),f1(dim1,1,0.);
    matred.K11Reduced(k11,f1);
    // matred.K01().Print("K01Orig =",std::cout,EMathematicaInput);
    int count = 0;
    int64_t gfemconnectindex = patch.fGFemConnectIndex;
    TPZConnect &cgfem = fMultiphysicsCompMesh->ConnectVec()[gfemconnectindex];
    int cgfemsize = cgfem.NShape()*cgfem.NState();
    for(auto it : patch.fH1ConnectIndexes) {
        int64_t h1connectindex = it.second;
        TPZConnect &ch1 = fMultiphysicsCompMesh->ConnectVec()[h1connectindex];
        int64_t ch1size = ch1.NShape()*ch1.NState();
        if(!ch1size) continue;
        TPZFMatrix<STATE> k01 = matred.K01();
        k01 = k01 * (-1.);
        ch1.SetCondensed(false);
        ch1.AddDependency(h1connectindex, gfemconnectindex, k01, count, 0, ch1size, cgfemsize);
        count += ch1size;
    }
    for(int i=dim0; i<dim; i++) {
        for(int j=0; j<dim0; j++) {
            matlocal(i,j) = 0.;
            matlocal(j,i) = 0.;
        }
    }
    for(int i=0; i<dim1; i++) {
        for(int j=0; j<dim1; j++) {
            matlocal(i+dim0,j+dim0) = k11(i,j);
        }
    }
    NormalizeDiagonal(matlocal);
    if(0)
    {
        TPZManVector<CSTATE,20> eigenvalues(dim);
        matlocal.SolveEigenProblem(eigenvalues);
        {
            REAL small = std::norm(eigenvalues[0]);
            for(auto it : eigenvalues) {
                if(std::norm(it) < small) small = std::norm(it);
            }
            patch.fOrthogonalizedSmallestEigenvalue = small;
        }
    }
    {
        TPZFNMatrix<10,STATE> eigenvector(dim,1,0.);
        for(int i=0; i<dim; i++) eigenvector(i,0) = rand()*1./RAND_MAX;
        TPZFNMatrix<50,STATE> matinv;
        REAL norm = Norm(eigenvector);
        // matlocal.Print(std::cout);
        eigenvector = eigenvector * (1./norm);
        matlocal.Inverse(matinv,ECholesky);
        // matinv.Print(std::cout);
        for(int it = 0; it<100; it++) {
            eigenvector = matinv*eigenvector;
            norm = Norm(eigenvector);
            eigenvector = eigenvector * (1./norm);
            // std::cout << "norm = " << norm << std::endl;
        }
        patch.fOrthogonalizedSmallestEigenvalue = 1./norm;

    }

}

    // Verify is the equations are orthogonalized
void TPZGFemOrthogonal::VerifyOrthogonality()
{
    int64_t neq_total = fMultiphysicsCompMesh->NEquations();
    // loop over the node patches
    for (auto &nodepatch : fNodePatches) {
        int64_t gfemconnectindex = nodepatch.fGFemConnectIndex;
        TPZConnect &cgfem = fMultiphysicsCompMesh->ConnectVec()[gfemconnectindex];
        int64_t dim1 = cgfem.NShape()*cgfem.NState();
        int64_t dim0 = 0;
        // loop over the H1 connects
        for (auto &h1connectpair : nodepatch.fH1ConnectIndexes) {
            // orthogonalize the connect
            int64_t h1connectindex = h1connectpair.first;
            TPZConnect &ch1 = fMultiphysicsCompMesh->ConnectVec()[h1connectindex];
            dim0 += ch1.NShape()*ch1.NState();
        }
        TPZFMatrix<STATE> patchmatrix(dim0+dim1,dim0+dim1,0.);
        ComputePatchMatrix(nodepatch, patchmatrix);
        TPZFMatrix<STATE> k01(dim0,dim1,0.);
        for(int i = 0; i < dim0; i++) {
            for(int j = 0; j < dim1; j++) {
                k01(i,j) = patchmatrix.GetVal(i,j+dim0);
            }
        }
        REAL norm = Norm(k01);
        if(norm > 1.e-6) 
        {
            std::cout << "Connect " << gfemconnectindex << " is not orthogonalized\n";
            k01.Print("k01 =",std::cout,EMathematicaInput);
            
        }
    }
}
