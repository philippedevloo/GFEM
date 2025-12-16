/**
 * @file
 * @brief Contains the TPZSSpStructMatrix class which implements sparse structural matrices.
 */

#ifndef TPZSSpMatRedStructMatrix_H
#define TPZSSpMatRedStructMatrix_H

#include "TPZSSpStructMatrix.h"
/**
 * @brief Implements a sparse symmetric structural matrix using TPZSYsmpMatrix as a storage format.
 * @ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZSSpMatRedStructMatrix : public TPZStructMatrixT<TVar>,
                                  public TPar {
    
public:
	using TPZStructMatrixT<TVar>::TPZStructMatrixT;

    TPZSSpMatRedStructMatrix() = default;

    TPZSSpMatRedStructMatrix(TPZCompMesh *mesh);

    TPZSSpMatRedStructMatrix(const TPZSSpMatRedStructMatrix &copy) : TPZStructMatrixT<TVar>(copy), fDim(copy.fDim), fDim0(copy.fDim0), fDim0Reduced(copy.fDim0Reduced), fDimReduced(copy.fDimReduced), fLagrangeLevels(copy.fLagrangeLevels) {
    }

    TPZSSpMatRedStructMatrix(TPZCompMesh *mesh, const std::set<int> &LagrangeLevels);

    /// @brief Create a structural matrix that will condense the first dim0 equations
    /// @param mesh computational mesh already reordered by connects
    /// @param dim0 number of equations that will be resolved by a direct solver
    TPZSSpMatRedStructMatrix(TPZCompMesh *mesh, int64_t dim0);

    TPZMatrix<TVar> * Create() override;
	TPZStructMatrix * Clone() override;

        /*! Creates solver matrix and assembles it alongside global rhs.
     Avoid overriding it unless there are no other options*/
    virtual TPZBaseMatrix *
    CreateAssemble(TPZBaseMatrix &rhs) override;

    void EndCreateAssemble(TPZBaseMatrix *mat) override;

    /// @brief return a reference to the K11 matrix
    static TPZAutoPointer<TPZMatrix<TVar>>  CloneK11(TPZAutoPointer<TPZMatrix<TVar>> matrix);

    /// @brief extract the reduced right hand side
    static void ExtractF1(TPZAutoPointer<TPZMatrix<TVar>> matrix, TPZFMatrix<TVar> &F1);

    /// resequence the equation as a function of the lagrange levels
    /// all connects with lagrange level included in the set will be ordered first
    void Resequence();
    //@{
    int ClassId() const override;

    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;
    //@}
protected:
    virtual TPZMatrix<TVar> * SetupMatrixData(TPZStack<int64_t> & elgraph, TPZVec<int64_t> &elgraphindex);
private :
    /// @brief Size of the system of equations before condensing
    int64_t fDim = -1;

    /// @brief Number of equations that will be condensed
    int64_t fDim0 = -1;

    /// @brief Size of the system of equations after applying equation filter
    int64_t fDimReduced = -1;

    /// @brief Number of condensed equations after applying equation filter
    int64_t fDim0Reduced = -1;

    /// @brief Lagrange levels which determine the condensed and non condensed equations
    std::set<int> fLagrangeLevels;
    
    friend TPZPersistenceManager;
};

#endif //TPZSSpStructMatrix_H
