//
// Created by Gustavo Batistela on 5/13/21.
//

#ifndef TPZGFEMDARCYFLOW_H
#define TPZGFEMDARCYFLOW_H

#include "DarcyFlow/TPZDarcyFlow.h"

/**
 * @ingroup material
 * @brief This class implements an H1-conforming approximation for the Darcy flow equation for isotropic materials.
 *
 * The Darcy flow equation is given by: \f[- \nabla \cdot K \nabla u = f,\f] where \f$u\f$ is the pressure field
 * to be solved, \f$K\f$ is the permeability tensor and \f$f\f$ is the source term.
 *
 * @see TPZMixedDarcyFlow For an approximation using the mixed method.
 * @see TPZIsotropicPermeability For setting the permeability field.
 */

#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"


class TPZGFemDarcyFlow : public TPZDarcyFlow, public TPZMatCombinedSpacesT<STATE>,public TPZMatErrorCombinedSpaces<STATE> {


public:
    /**
     * @brief Default constructor
     */
    TPZGFemDarcyFlow();

    /**
	 * @brief Class constructor
	 * @param [in] id material id
	 * @param [in] dim problem dimension
	 */
    TPZGFemDarcyFlow(int id, int dim);

    TPZGFemDarcyFlow(const TPZGFemDarcyFlow &copy) : TPZDarcyFlow(copy), TPZMatCombinedSpacesT<STATE>(copy),
    TPZMatErrorCombinedSpaces<STATE>(copy){

    }
    TPZGFemDarcyFlow &operator=(const TPZGFemDarcyFlow &copy){
        TPZDarcyFlow::operator=(copy);
        TPZMatCombinedSpacesT<STATE>::operator=(copy);
        TPZMatErrorCombinedSpaces<STATE>::operator=(copy);
        return *this;
    }
    /**
	 * @brief Returns a 'std::string' with the name of the material
	 */
    [[nodiscard]] std::string Name() const override { return "TPZGFemDarcyFlow"; }


    /** @name Contribute */
    /** @{ */
    /**
     * @brief It computes a contribution to the stiffness matrix
     * and load vector at one integration point.
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                            REAL weight,TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override;
    /**@}*/

    /** @name ContributeBC
        @ingroup Contribute*/
    /**@{*/
    /**
     * @brief It computes a contribution to the stiffness matrix
     * and load vector at one BC integration point.
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                              REAL weight, TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef,
                              TPZBndCondT<STATE> &bc) override;

    /**@}*/
    /** @brief Returns the solution associated with a given index
        based on the finite element approximation.
        @param[in] datavec Stores all the input data.
        @param[in] var Index of the queried solution
        @param[out] sol FEM Solution at the integration point
    */
    virtual void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                          int var, TPZVec<STATE> &sol) override {
        TPZDarcyFlow::Solution(datavec[0],var,sol);
        for(int i = 1; i < datavec.size(); i++){
            TPZManVector<STATE> soltemp(sol.size(),0.);
            TPZDarcyFlow::Solution(datavec[i],var,soltemp);
            for(int j = 0; j < sol.size(); j++){
                sol[j] += soltemp[j];
            }
        }
    }

    //! @name Error
    /** @{*/
    /*!
      \brief Calculates the error at a given point x.
      \param[in] datavec input data
      \param[out] errors The calculated errors.
     */
    virtual void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data,
                        TPZVec<REAL> &errors) override;

    /**
     * @brief Returns an unique class identifier
     */
    [[nodiscard]] int ClassId() const override;
    
    /** @brief Writes this object to the TPZStream buffer. Include the classid if `withclassid = true` */
    virtual void Write(TPZStream &buf, int withclassid) const override
    {
        
    }
    
    /** @brief Reads an objects from the TPZStream buffer. */
    virtual void Read(TPZStream &buf, void *context) override
    {
        
    }


    /**
     * @brief Creates another material of the same type
     */
    [[nodiscard]] TPZMaterial *NewMaterial() const override;

    /**
     * @brief Prints data associated with the material.
     */
    void Print(std::ostream & out) const override;

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
        return new  TPZBndCondBase<STATE,TPZMatCombinedSpacesBC<STATE>, TPZMatErrorCombinedSpacesBC<STATE> >
        (reference,id, type,val1,val2);
    }

};


#endif //TPZDARCYFLOW_H
