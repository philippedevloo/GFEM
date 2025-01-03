/**
 * @file
 * @brief Contains the implementation of the TPZSSpMatRedStructMatrix methods. 
 */

#include "TPZSSpMatRedStructMatrix.h"
#include "pzcmesh.h"
#include "TPZYSMPMatrix.h"
#include "TPZRenumbering.h"
#include "TPZGuiInterface.h"
#include "pzcmesh.h"
#include "TPZSSparseMatRed.h"
#include "TPZTimer.h"
#include "pzstepsolver.h"

#ifdef USING_EIGEN
#include "TPZEigenSparseMatrix.h"
#endif

#ifdef USING_MKL
#include "TPZYSMPPardiso.h"
#endif


#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.StrMatrix");
#endif

using namespace std;

template<class TVar, class TPar>
TPZSSpMatRedStructMatrix<TVar,TPar>::TPZSSpMatRedStructMatrix(TPZCompMesh *mesh) : TPZStructMatrixT<TVar>(mesh) {
    fDim = mesh->NEquations();
    // no equations are condensed
    fDim0 = 0;
}

template<class TVar, class TPar>
TPZSSpMatRedStructMatrix<TVar,TPar>::TPZSSpMatRedStructMatrix(TPZCompMesh *mesh, const std::set<int> &LagrangeLevels) : TPZStructMatrixT<TVar>(mesh), fLagrangeLevels(LagrangeLevels) {
    fDim = mesh->NEquations();
    fDim0 = 0;
    Resequence();
}

/// @brief Create a structural matrix that will condense the first dim0 equations
/// @param mesh computational mesh already reordered by connects
/// @param dim0 number of equations that will be resolved by a direct solver
template<class TVar, class TPar>
TPZSSpMatRedStructMatrix<TVar,TPar>::TPZSSpMatRedStructMatrix(TPZCompMesh *mesh, int64_t dim0): TPZStructMatrixT<TVar>(mesh)
{
    fDim = mesh->NEquations();
    fDim0 = dim0;
}



template<class TVar, class TPar>
TPZStructMatrix * TPZSSpMatRedStructMatrix<TVar,TPar>::Clone(){
    return new TPZSSpMatRedStructMatrix(*this);
}

    /*! Creates solver matrix and assembles it alongside global rhs.
    Avoid overriding it unless there are no other options*/
template<class TVar, class TPar>
TPZBaseMatrix *
TPZSSpMatRedStructMatrix<TVar,TPar>::CreateAssemble(TPZBaseMatrix &rhs,
                   TPZAutoPointer<TPZGuiInterface> guiInterface) {
    auto result = TPar::CreateAssemble(rhs,guiInterface);
     auto spMat = dynamic_cast<TPZSparseMatRed<TVar> *>(result);
    if(!spMat) DebugStop();
    TPZFMatrix<TVar> *rhsF = dynamic_cast<TPZFMatrix<TVar> *>(&rhs);
    if(!(rhsF)) DebugStop();
    spMat->SetF(*rhsF);
    return result;
}

template<class TVar, class TPar>
void TPZSSpMatRedStructMatrix<TVar,TPar>::EndCreateAssemble(TPZBaseMatrix * mat){
    auto spMat = dynamic_cast<TPZSparseMatRed<TVar> *>(mat);
    if(!spMat) DebugStop();
    spMat->SimetrizeMatRed();
    spMat->SetReduced();
}

template<class TVar, class TPar>
TPZMatrix<TVar> * TPZSSpMatRedStructMatrix<TVar,TPar>::Create(){
    int64_t neq = this->fMesh->NEquations();

	
    /**
     *Longhin implementation
	 */
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> elgraphindex;
    //    int nnodes = 0;
    this->fMesh->ComputeElGraph(elgraph,elgraphindex);
    // this will also compute the number of reduced equations
    TPZMatrix<TVar> * mat = SetupMatrixData(elgraph,elgraphindex);

    TPZSparseMatRed<TVar> *sparsemat = new TPZSparseMatRed<TVar>(fDimReduced, fDim0Reduced);
    sparsemat->AllocateSubMatrices(*mat);
    sparsemat->InitializeSolver(ELDLt);
    delete mat;
    return sparsemat;
}

/// @brief return a reference to the K11 matrix
template<class TVar, class TPar>
TPZAutoPointer<TPZMatrix<TVar>>  TPZSSpMatRedStructMatrix<TVar,TPar>::CloneK11(TPZAutoPointer<TPZMatrix<TVar>> matrix) {
    auto *sparsemat = dynamic_cast<TPZSparseMatRed<TVar> *>(matrix.operator->());
    if (!sparsemat) {
        DebugStop();
    }
    TPZSYsmpMatrix<TVar> * matK11 = dynamic_cast<TPZSYsmpMatrix<TVar> *>(&sparsemat->K11());
    if (!matK11) {
        DebugStop();
    }
#ifdef USING_MKL
    TPZFYsmpMatrixPardiso<TVar> * mat = new TPZFYsmpMatrixPardiso<TVar>(*matK11);
#elif USING_EIGEN
    TPZEigenSparseMatrix<TVar> * mat = new TPZEigenSparseMatrix<TVar>(*matK11);
#else
    TPZFYsmpMatrix<TVar> *mat = new TPZFYsmpMatrix<TVar>(*matK11);
    DebugStop();
#endif

    TPZAutoPointer<TPZMatrix<TVar>> K11(mat);
    return K11;

}

/// @brief extract the reduced right hand side
template<class TVar, class TPar>
void TPZSSpMatRedStructMatrix<TVar,TPar>::ExtractF1(TPZAutoPointer<TPZMatrix<TVar>> matrix, TPZFMatrix<TVar> &F1){
    auto *sparsemat = dynamic_cast<TPZSparseMatRed<TVar> *>(matrix.operator->());
    if (!sparsemat) {
        DebugStop();
    }
    sparsemat->F1Red(F1);
}



template<class TVar, class TPar>
TPZMatrix<TVar> *
TPZSSpMatRedStructMatrix<TVar,TPar>::SetupMatrixData(TPZStack<int64_t> & elgraph,
                                              TPZVec<int64_t> &elgraphindex){
    
    const int64_t neq = this->fEquationFilter.NActiveEquations();
#ifdef USING_MKL
    TPZFYsmpMatrixPardiso<TVar> * mat = new TPZFYsmpMatrixPardiso<TVar>(neq,neq);
#elif USING_EIGEN
    TPZEigenSparseMatrix<TVar> * mat = new TPZEigenSparseMatrix<TVar>(neq,neq);
#else
    TPZFYsmpMatrix<TVar> *mat = new TPZFYsmpMatrix<TVar>(neq,neq);
    DebugStop();
#endif
    if(this->fEquationFilter.IsActive()) {
        DebugStop();
    } else {
        fDimReduced = fDim;
        fDim0Reduced = fDim0;
    }
    /**Creates a element graph*/
    TPZRenumbering graph;
    graph.SetElementsNodes(elgraphindex.NElements() -1 ,
                           this->fMesh->NIndependentConnects());
	
    TPZManVector<int64_t> nodegraph;
    TPZManVector<int64_t> nodegraphindex;
    /**
     *converts an element graph structure into a node graph structure
     *those vectors have size ZERO !!!
     */
    graph.ConvertGraph(elgraph,elgraphindex,nodegraph,nodegraphindex);
    /**vector sizes*/
    const int64_t nblock = nodegraphindex.NElements()-1;
    // number of values in the sparse matrix
    int64_t totalvar = 0;
    // number of equations
    int64_t totaleq = 0;
    for(int64_t i=0;i<nblock;i++){
		const int64_t iblsize = this->fMesh->Block().Size(i);
		const int64_t iblpos = this->fMesh->Block().Position(i);
        const int64_t numactive = this->fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
		totaleq += iblsize;
		const int64_t icfirst = nodegraphindex[i];
		const int64_t iclast = nodegraphindex[i+1];
		totalvar+=iblsize*iblsize;
		for(auto j=icfirst;j<iclast;j++) {
			const int64_t col = nodegraph[j];
			const int64_t colsize = this->fMesh->Block().Size(col);
			const int64_t colpos = this->fMesh->Block().Position(col);
            const int64_t numactive = this->fEquationFilter.NumActive(colpos, colpos+colsize);
            if (!numactive) {
                continue;
            }
			totalvar += iblsize*colsize;
		}
    }
	
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Number of equations " << totaleq << " number of nonzero s " << totalvar;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    int64_t ieq = 0;
    // pos is the position where we will put the column value
    int64_t pos = 0;
	
    TPZManVector<int64_t,400> Eq(totaleq+1);
    TPZVec<int64_t> EqCol(totalvar);
    TPZVec<TVar> EqValue(totalvar);
    //lambda for avoid repeating code
    auto AddColEqs =
        [this,&EqCol,&EqValue,&pos](const int colsize, const int colpos)
        {
            TPZManVector<int64_t> destindices(colsize);
            for (int64_t i=0; i<colsize; i++) {
                destindices[i] = colpos+i;
            }
            this->fEquationFilter.Filter(destindices);
            for(auto jbleq=0; jbleq<destindices.size(); jbleq++) {
                EqCol[pos] = destindices[jbleq];
                EqValue[pos] = 0.;
                pos++;
            }
        };

    for(auto i=0;i<nblock;i++){
		const int64_t iblsize = this->fMesh->Block().Size(i);
		const int64_t iblpos = this->fMesh->Block().Position(i);
        const int64_t numactive =
            this->fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
        TPZManVector<int64_t> rowdestindices(iblsize);
        for (int64_t ij=0; ij<iblsize; ij++) {
            rowdestindices[ij] = iblpos+ij;
        }
        this->fEquationFilter.Filter(rowdestindices);
        // working equation by equation
        // rowdestindices contains the equation number of each element in the block number "i"
		for(auto ibleq=0; ibleq<rowdestindices.size(); ibleq++) {
            if (rowdestindices[ibleq] != ieq) {
                DebugStop();
            }
			Eq[ieq] = pos;
            //we need to add the diagonal block
    //since nodegraph does not contain the pointer to itself
            AddColEqs(this->fMesh->Block().Size(i),
                      this->fMesh->Block().Position(i));
			const int64_t icfirst = nodegraphindex[i];
			const int64_t iclast = nodegraphindex[i+1];
			for(auto j=icfirst;j<iclast;j++)
            {
                // col is the block linked to block "i"
				const int64_t col = nodegraph[j];
				const auto colsize = this->fMesh->Block().Size(col);
				const auto colpos = this->fMesh->Block().Position(col);
                const int64_t numactive =
                    this->fEquationFilter.NumActive(colpos, colpos+colsize);
                if (!numactive) {
                    continue;
                }
                AddColEqs(colsize,colpos);
			}
            //we need to sort the entries of ja otherwise pardiso will complain
            const auto firstentry = Eq[ieq];
            const auto nentries = pos-Eq[ieq];

            std::stable_sort(&EqCol[firstentry],
                             &EqCol[firstentry]+nentries);
            
			ieq++;
		}
    }
    Eq[ieq] = pos;
    if(pos != totalvar)
    {
        PZError<<__PRETTY_FUNCTION__<<'\n';
        PZError<<"nvars: "<<totalvar<<'\n';
        PZError<<"l pos: "<<pos<<'\n';
        DebugStop();
    }
    mat->SetData(std::move(Eq),std::move(EqCol),std::move(EqValue));
    return mat;
}

/// resequence the equation as a function of the lagrange levels
/// all connects with lagrange level included in the set will be ordered first
template<class TVar, class TPar>
void TPZSSpMatRedStructMatrix<TVar,TPar>::Resequence() {
    TPZCompMesh *cmesh = this->fMesh;
    cmesh->LoadReferences();
    auto neqsTotal = cmesh->NEquations();
  
    int64_t dim00 = 0;
    int64_t dim11 = 0;
  
    // loop over the connects and create two std::sets one to the "pressure" ones, representing the
    // degrees of freedom to be condensed. The second set contains the "flux" connects, which will not be condensed
    
    int64_t seqNumP = 0;
    cmesh->Block().Resequence();

    // that is soo strange. Putting the sequence number according to connect index??
    for (int icon = 0; icon < cmesh->NConnects(); icon++){
        
        TPZConnect &con = cmesh->ConnectVec()[icon];
        int64_t seqNum = con.SequenceNumber();
        if (con.IsCondensed()) continue;
        if (seqNum < 0) continue;
        int laglevel = con.LagrangeMultiplier();

        if(fLagrangeLevels.find(laglevel) != fLagrangeLevels.end()) {
            con.SetSequenceNumber(seqNumP);
            //Em cada caso precisa atualizar o tamanho do bloco
            int neq = con.NShape()*con.NState();
            seqNumP++;
            dim00 += neq;
            seqNum=con.SequenceNumber();
            cmesh->Block().Set(seqNum,neq);
        }
    }
    
    std::set<int64_t> auxConnects00, auxConnects11;
    int64_t seqNumF = seqNumP;
    
    //Second - Set the sequence number to the flux variables - which will not be condensed
    for (int icon = 0; icon < cmesh->NConnects(); icon++){
        
        TPZConnect &con = cmesh->ConnectVec()[icon];
        int64_t seqNum = con.SequenceNumber();
        if (con.IsCondensed()) continue;
        if (seqNum < 0) continue;
        
        int laglevel = con.LagrangeMultiplier();

        if(fLagrangeLevels.find(laglevel) == fLagrangeLevels.end()) {

            con.SetSequenceNumber(seqNumF);
            
            //Em cada caso precisa atualizar o tamanho do bloco
            int neq = con.NShape()*con.NState();
            seqNumF++;
            dim11 += neq;
            seqNum=con.SequenceNumber();
            cmesh->Block().Set(seqNum,neq);
        }
    }
    cmesh->Block().Resequence();
    cmesh->ExpandSolution();
    
    
    fDim = dim00+dim11;
    fDim0 = dim00;
}


template<class TVar, class TPar>
int TPZSSpMatRedStructMatrix<TVar,TPar>::ClassId() const{
    return Hash("TPZSSpMatRedStructMatrix") ^
        TPZStructMatrixT<TVar>::ClassId() << 1 ^
        TPar::ClassId() << 2;
}

template<class TVar, class TPar>
void TPZSSpMatRedStructMatrix<TVar,TPar>::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPar::Read(buf,context);
}

template<class TVar, class TPar>
void TPZSSpMatRedStructMatrix<TVar,TPar>::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPar::Write(buf,withclassid);
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"

template class TPZSSpMatRedStructMatrix<STATE,TPZStructMatrixOR<STATE>>;
template class TPZSSpMatRedStructMatrix<STATE,TPZStructMatrixOT<STATE>>;
template class TPZSSpMatRedStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>;

#ifndef USING_EIGEN
template class TPZSSpMatRedStructMatrix<CSTATE,TPZStructMatrixOR<CSTATE>>;
template class TPZSSpMatRedStructMatrix<CSTATE,TPZStructMatrixOT<CSTATE>>;
template class TPZSSpMatRedStructMatrix<CSTATE,TPZStructMatrixTBBFlow<CSTATE>>;
#endif
