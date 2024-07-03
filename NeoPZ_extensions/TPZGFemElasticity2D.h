/**
 * @file
 * @brief Contains the TPZElasticity2D class which implements a two dimensional elastic material in plane stress or strain.
 */

#ifndef TPZGFEMELASTICITY2D_H
#define TPZGFEMELASTICITY2D_H


#include "Elasticity/TPZElasticity2D.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"


/**
 * @ingroup material
 * @brief This class implements a two dimensional elastic material in plane stress or strain
 */
class TPZGFemElasticity2D : public TPZMatCombinedSpacesT<STATE>, TPZMatErrorCombinedSpaces<STATE>, public TPZElasticity2D
{
   
public :
	/**
	 * @brief Creates an elastic material with:
	 * @param id material id
	 * @param E elasticity modulus
	 * @param nu poisson coefficient
	 * @param fx forcing function \f$ -x = fx \f$ 
	 * @param fy forcing function \f$ -y = fy \f$
	 * @param plainstress \f$ plainstress = 1 \f$ indicates use of plainstress
	 */
	TPZGFemElasticity2D(int id, STATE E, STATE nu, STATE fx, STATE fy, int planestress = 1) :
    TPZElasticity2D(id,E,nu,fx,fy,planestress)
    {
        
    }
    
    TPZGFemElasticity2D(int id) : TPZElasticity2D(id) {}
    
    TPZGFemElasticity2D() : TPZElasticity2D(){}

    TPZGFemElasticity2D(const TPZGFemElasticity2D &copy) : TPZElasticity2D(copy), TPZMatCombinedSpacesT<STATE>(copy),
    TPZMatErrorCombinedSpaces<STATE>(copy){

    }
	/** @brief Returns the material name*/
	std::string Name()  const override { return "TPZGFemElasticity2D"; }
		
	/** @name Contribute methods */
	/** @{ */
	
	/** @brief Calculates the element stiffness matrix */
	void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &data, STATE weight,
                    TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;

    /** @brief Applies the element boundary conditions */
	void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &data,STATE weight,
                      TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;

    // Compute the D matrix for the element contribution
    void ComputeDMatrix(TPZMaterialData &data, TPZFMatrix<STATE> &D) const;

    /*
     * @brief Fill requirements for volumetric contribute
     */
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    /*
     * @brief Fill requirements for boundary contribute
     */
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

	
	/** @} */

    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &data, int var,
                  TPZVec<STATE> &Solout) override;
    
	/** @} */

    /** @brief Creates a new material from the current object   ??*/
	TPZMaterial * NewMaterial()  const override;
    /**
       @name ReadWrite
       @{*/
    //!Read and Write methods
    int ClassId() const override;

	void Read(TPZStream &buf, void *context) override;
	
	void Write(TPZStream &buf, int withclassid) const override;
    /**@}*/
	/** 
	 * @brief Computes the error due to the difference between the interpolated flux \n
	 * and the flux computed based on the derivative of the solution
	 */

    /** @name Errors */
	/** @{ */
	void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data,
                TPZVec<STATE> &values) override;
    /** @} */
    /** @brief Creates an associated boundary condition.
     @param[in] reference The volumetric material associated with the BC.
     @param[in] id Boundary condition identifier.
     @param[in] type Type of the boundary condition.
     @param[in] val1 Value to be set at the element matrix.
     @param[in] val2 Value to be set at the rhs vector.
    */
    virtual TPZBndCondT<STATE>* CreateBC(TPZMaterial *reference,
                                        int id, int type,
                                        const TPZFMatrix<STATE> &val1,
                                        const TPZVec<STATE> &val2) override
    {
        return new  TPZBndCondBase<STATE,TPZMatCombinedSpacesBC<STATE>, TPZMatErrorCombinedSpacesBC<STATE>,TPZMatLoadCasesBC<STATE> >
        (reference,id, type,val1,val2);
    }

    /// transform a scalar shape structure to a vectorial shape structure
    // this method will expand phi and dphix to a vector representation
    static void AtomicToVec(int nstate, TPZMaterialData &data);

    /// transform a vectorial shape structure to a scalar shape structure
    /// this method will only resize the data structure
    static void VecToAtomic(int nstate, TPZMaterialData &data);

    
    // transform derivatives in the direction of axes to derivatives in xy
    static void VecAxesToXYZ(int nstate, TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &dudaxes, TPZFMatrix<REAL> &dudx);

    // transform derivatives in the direction of xy to derivatives in axes
    static void VecXYZToAxes(int nstate, TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &dudx, TPZFMatrix<REAL> &dudaxes);
};

#endif
